#ifndef TATAMI_HDF5_SPARSE_MATRIX_HPP
#define TATAMI_HDF5_SPARSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../../base/Matrix.hpp"
#include "../../base/CompressedSparseMatrix.hpp"
#include "utils.hpp"

/**
 * @file HDF5CompressedSparseMatrix.hpp
 *
 * @brief Defines a class for a HDF5-backed compressed sparse matrix.
 */

namespace tatami {

/**
 * @brief Compressed sparse matrix in a HDF5 file.
 *
 * This class retrieves sparse data from the HDF5 file on demand rather than loading it all in at the start.
 * This allows us to handle very large datasets in limited memory at the cost of speed.
 *
 * We manually handle the chunk caching to speed up access for consecutive rows or columns (for compressed sparse row and column matrices, respectively).
 * The policy is to minimize the number of calls to the HDF5 library by requesting large contiguous slices where possible.
 * The size of the slice is determined by the cache limit in the constructor.
 *
 * Callers should follow the `prefer_rows()` suggestion when extracting data,
 * as this tries to minimize the number of chunks that need to be read per access request.
 * This recommendation is even stronger than for the `HDF5DenseMatrix`,
 * as the access pattern on disk for the non-preferred dimension is very suboptimal.
 *
 * As the HDF5 library is not generally thread-safe, the HDF5-related operations should only be run in a single thread.
 * For OpenMP, this is handled automatically by putting all HDF5 operations in a critical region.
 * For other parallelization schemes, callers should define the `TATAMI_HDF5_PARALLEL_LOCK` macro;
 * this should be a function that accepts and executes a no-argument lambda within an appropriate serial region (e.g., based on a global mutex).
 *
 * @tparam ROW Whether the matrix is stored in compressed sparse row format.
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 */
template<bool ROW, typename T, typename IDX>
class HDF5CompressedSparseMatrix : public tatami::Matrix<T, IDX> {
    size_t nrows, ncols;
    std::string file_name;
    std::string data_name, index_name;
    std::vector<hsize_t> pointers;

    std::vector<size_t> primary_cache_id;
    std::vector<std::pair<size_t, size_t> > primary_cache_limits;

public:
    /**
     * @param nr Number of rows in the matrix.
     * @param nc Number of columns in the matrix.
     * @param file Path to the file.
     * @param vals Name of the 1D dataset inside `file` containing the non-zero elements.
     * @param idx Name of the 1D dataset inside `file` containing the indices of the non-zero elements.
     * If `ROW = true`, this should contain column indices sorted within each row, otherwise it should contain row indices sorted within each column.
     * @param ptr Name of the 1D dataset inside `file` containing the index pointers for the start and end of each row (if `ROW = true`) or column (otherwise).
     * This should have length equal to the number of rows (if `ROW = true`) or columns (otherwise).
     * @param cache_limit Limit to the size of the cache, in bytes.
     *
     * The cache is created by extracting multiple columns (for CSC matrices) or rows (CSR) on every call to the HDF5 library.
     * These are held in memory in the `Workspace` while the relevant column/row is returned to the user by `row()` or `column()`.
     * The aim is to minimize the number of calls to the HDF5 library - and thus expensive file reads - for consecutive accesses.
     */
    HDF5CompressedSparseMatrix(size_t nr, size_t nc, std::string file, std::string vals, std::string idx, std::string ptr, size_t cache_limit = 100000000) :
        nrows(nr),
        ncols(nc),
        file_name(file),
        data_name(std::move(vals)),
        index_name(std::move(idx)),
        pointers(ROW ? nr + 1 : nc + 1)
    {
        size_t total_element_size;

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        H5::H5File file_handle(file_name, H5F_ACC_RDONLY);

        auto dhandle = HDF5::open_and_check_dataset<false>(file_handle, data_name);
        const size_t nonzeros = HDF5::get_array_dimensions<1>(dhandle, "vals")[0];

        auto ihandle = HDF5::open_and_check_dataset<true>(file_handle, index_name);
        if (HDF5::get_array_dimensions<1>(ihandle, "idx")[0] != nonzeros) {
            throw std::runtime_error("number of non-zero elements is not consistent between 'data' and 'idx'");
        }

        auto phandle = HDF5::open_and_check_dataset<true>(file_handle, ptr);
        const size_t ptr_size = HDF5::get_array_dimensions<1>(phandle, "ptr")[0];
        if (ptr_size != pointers.size()) {
            throw std::runtime_error("'ptr' dataset should have length equal to the number of " + (ROW ? std::string("rows") : std::string("columns")) + " plus 1");
        }

        // Checking the contents of the index pointers.
        phandle.read(pointers.data(), H5::PredType::NATIVE_HSIZE);
        if (pointers[0] != 0) {
            throw std::runtime_error("first index pointer should be zero");
        }
        if (pointers.back() != nonzeros) {
            throw std::runtime_error("last index pointer should be equal to the number of non-zero elements");
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif

        // Setting up the cache parameters.
        size_t primary_dim = pointers.size() - 1;
        primary_cache_id.resize(primary_dim);
        if (primary_dim) {
            primary_cache_limits.emplace_back(0, pointers[1]);
        }

        size_t effective_cache_limit = cache_limit / (sizeof(T) + sizeof(IDX));
        size_t counter = 0, start = 0;
        for (size_t i = 1; i < primary_dim; ++i) {
            size_t end = pointers[i + 1];
            if (end - start <= effective_cache_limit) {
                primary_cache_limits.back().second = end;
            } else {
                ++counter;
                start = pointers[i];
                primary_cache_limits.emplace_back(start, end);
            }
            primary_cache_id[i] = counter;
        }

        return;
    }

public:
    size_t nrow() const {
        return nrows;
    }

    size_t ncol() const {
        return ncols;
    }

    /**
     * @return `true`.
     */
    bool sparse() const {
        return true;
    }

    /**
     * @return `true` if this is in compressed sparse row format.
     */
    bool prefer_rows() const {
        return ROW;
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::sparse_row;

    using Matrix<T, IDX>::sparse_column;

private:
    struct SparseWorkspaceBase {
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        std::vector<IDX> index_cache;
    };

    void fill_base(SparseWorkspaceBase& base) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        // TODO: set more suitable chunk cache values here, to avoid re-reading
        // chunks on the boundaries of the primary cache.
        base.file.openFile(file_name, H5F_ACC_RDONLY);

        base.data = base.file.openDataSet(data_name);
        base.index = base.file.openDataSet(index_name);
        base.dataspace = base.data.getSpace();

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif
    }

    struct PrimarySparseWorkspaceBase : public SparseWorkspaceBase {
        std::vector<T> data_cache;
        size_t current_cache_id = 0;
        bool init = false;

        std::vector<size_t> starts;
    };

    void fill_base(PrimarySparseWorkspaceBase& base, size_t cache_size) const {
        fill_base(base);
        base.starts.resize(cache_size, -1);
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct HDF5PrimarySparseWorkspace : public Workspace<WORKROW> {
        PrimarySparseWorkspaceBase base;
    };

    template<bool WORKROW>
    struct HDF5SecondarySparseWorkspace : public Workspace<WORKROW> {
        SparseWorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace(bool = false) const {
        if constexpr(ROW) {
            auto ptr = new HDF5PrimarySparseWorkspace<true>;
            std::shared_ptr<RowWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5SecondarySparseWorkspace<true>;
            std::shared_ptr<RowWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace(bool = false) const {
        if constexpr(ROW) {
            auto ptr = new HDF5SecondarySparseWorkspace<false>;
            std::shared_ptr<ColumnWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5PrimarySparseWorkspace<false>;
            std::shared_ptr<ColumnWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(ROW) {
            return extract_primary(r, buffer, 0, ncols, static_cast<HDF5PrimarySparseWorkspace<true>*>(work)->base);
        } else {
            return extract_secondary(r, buffer, 0, ncols, static_cast<HDF5SecondarySparseWorkspace<true>*>(work)->base);
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(ROW) {
            return extract_secondary(c, buffer, 0, nrows, static_cast<HDF5SecondarySparseWorkspace<false>*>(work)->base);
        } else {
            return extract_primary(c, buffer, 0, nrows, static_cast<HDF5PrimarySparseWorkspace<false>*>(work)->base);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, RowWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_primary(r, dbuffer, ibuffer, 0, ncols, static_cast<HDF5PrimarySparseWorkspace<true>*>(work)->base);
        } else {
            return extract_secondary(r, dbuffer, ibuffer, 0, ncols, static_cast<HDF5SecondarySparseWorkspace<true>*>(work)->base);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_secondary(c, dbuffer, ibuffer, 0, nrows, static_cast<HDF5SecondarySparseWorkspace<false>*>(work)->base);
        } else {
            return extract_primary(c, dbuffer, ibuffer, 0, nrows, static_cast<HDF5PrimarySparseWorkspace<false>*>(work)->base);
        }
    }

private:
    void populate_primary_cache(size_t i, PrimarySparseWorkspaceBase& work) const {
        if (primary_cache_id[i] == work.current_cache_id && work.init) {
            return;
        }

        work.init = true;
        work.current_cache_id = primary_cache_id[i];

        // Pulling out the entire chunk of primary dimensions containing
        // 'i'. We have to do this for all indices, regardless of the
        // slicing/indexing. 
        const auto& limits = primary_cache_limits[work.current_cache_id];
        hsize_t offset = limits.first;
        hsize_t count = limits.second - limits.first;

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        work.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        work.memspace.setExtentSimple(1, &count);
        work.memspace.selectAll();

        work.index_cache.resize(count);
        work.index.read(work.index_cache.data(), HDF5::define_mem_type<IDX>(), work.memspace, work.dataspace);

        // In theory, we could avoid extracting data for the entire column when
        // populating the cache for sliced or indexed queries. In practice,
        // each chunk is often larger than the number of non-zeroes in a
        // column, so we end up having to pull out the entire column anyway.
        // Also, I want to get out of this critical region ASAP and do my
        // indexing elsewhere. 

        work.data_cache.resize(count);
        work.data.read(work.data_cache.data(), HDF5::define_mem_type<T>(), work.memspace, work.dataspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif
    }

    template<class Function>
    void extract_primary(size_t i, Function fill, size_t start, size_t length, PrimarySparseWorkspaceBase& work) const {
        if (length == 0) {
            return;
        }

        size_t primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, work);

        const auto& limits = primary_cache_limits[work.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = work.index_cache.begin() + offset;
        auto dstart = work.data_cache.begin() + offset;

        size_t request_start = 0;
        if (start > *istart) { // 'istart' guaranteed to be valid from length check above.
            bool do_cache = !work.starts.empty();
            if (do_cache && work.starts[i] != -1) {
                request_start = work.starts[i];
            } else {
                request_start = std::lower_bound(istart, istart + primary_length, start) - istart;
                if (do_cache) {
                    work.starts[i] = request_start;
                }
            }
        }

        istart += request_start;
        dstart += request_start;

        size_t end = start + length;
        for (size_t i = request_start; i < primary_length && *istart < end; ++i, ++istart, ++dstart) {
            fill(*istart, *dstart);
        }

        return;
    }

    const T* extract_primary(size_t i, T* dbuffer, size_t start, size_t length, PrimarySparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + length, 0);
        extract_primary(i, [&](IDX pos, T value) -> void {
            dbuffer[pos - start] = value;
        }, start, length, work);
        return dbuffer;
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, size_t start, size_t length, PrimarySparseWorkspaceBase& work) const {
        size_t counter = 0;

        extract_primary(i, [&](IDX pos, T value) -> void {
            dbuffer[counter] = value;
            ibuffer[counter] = pos;
            ++counter;
        }, start, length, work);

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

private:
    // This could be improved by extracting multiple rows at any given call and
    // caching them for subsequent requests. However, even then, we'd require
    // multiple re-reads from file when we exceed the cache. So, any caching
    // would be just turning an extremely bad access pattern into a very bad
    // pattern, when users shouldn't even be calling this at all... 

    template<class Function>
    bool extract_secondary_raw(IDX primary, size_t secondary, Function& fill, SparseWorkspaceBase& work) const {
        size_t left = pointers[primary], right = pointers[primary + 1];
        work.index_cache.resize(right - left);

        // Serial locks should be applied by the callers.
        hsize_t offset = left;
        hsize_t count = work.index_cache.size();
        work.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        work.memspace.setExtentSimple(1, &count);
        work.memspace.selectAll();
        work.index.read(work.index_cache.data(), HDF5::define_mem_type<IDX>(), work.memspace, work.dataspace);

        auto it = std::lower_bound(work.index_cache.begin(), work.index_cache.end(), secondary);
        if (it != work.index_cache.end() && *it == secondary) {
            offset = left + (it - work.index_cache.begin());
            count = 1;
            work.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
            work.memspace.setExtentSimple(1, &count);
            work.memspace.selectAll();

            T dest;
            work.data.read(&dest, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
            fill(primary, dest);
            return true;
        } else {
            return false;
        }
    }

    template<class Function>
    void extract_secondary(size_t i, Function fill, size_t start, size_t length, SparseWorkspaceBase& work) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        size_t end = start + length;
        for (size_t j = start; j < end; ++j) {
            extract_secondary_raw(j, i, fill, work);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif
    }

    const T* extract_secondary(size_t i, T* dbuffer, size_t start, size_t length, SparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + length, 0);
        extract_secondary(i, [&](IDX pos, T value) -> void {
            dbuffer[pos - start] = value;
        }, start, length, work);
        return dbuffer;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, size_t start, size_t length, SparseWorkspaceBase& work) const {
        size_t counter = 0;

        extract_secondary(i, [&](IDX pos, T value) -> void {
            dbuffer[counter] = value;
            ibuffer[counter] = pos;
            ++counter;
        }, start, length, work);

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct HDF5PrimarySparseBlockWorkspace : public BlockWorkspace<WORKROW> {
        HDF5PrimarySparseBlockWorkspace(size_t s, size_t l) : details(s, l) {}
        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }
        PrimarySparseWorkspaceBase base;
    };

    template<bool WORKROW>
    struct HDF5SecondarySparseBlockWorkspace : public BlockWorkspace<WORKROW> {
        HDF5SecondarySparseBlockWorkspace(size_t s, size_t l) : details(s, l) {}
        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }
        SparseWorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t s, size_t l, bool cache = false) const {
        if constexpr(ROW) {
            auto ptr = new HDF5PrimarySparseBlockWorkspace<true>(s, l);
            std::shared_ptr<RowBlockWorkspace> output(ptr);
            fill_base(ptr->base, cache ? nrows : 0);
            return output;
        } else {
            auto ptr = new HDF5SecondarySparseBlockWorkspace<true>(s, l);
            std::shared_ptr<RowBlockWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t s, size_t l, bool cache = false) const {
        if constexpr(ROW) {
            auto ptr = new HDF5SecondarySparseBlockWorkspace<false>(s, l);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5PrimarySparseBlockWorkspace<false>(s, l);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);
            fill_base(ptr->base, cache ? ncols : 0);
            return output;
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<true>*>(work);
            auto details = wptr->details;
            return extract_primary(r, buffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<true>*>(work);
            auto details = wptr->details;
            return extract_secondary(r, buffer, details.first, details.second, wptr->base);
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<false>*>(work);
            auto details = wptr->details;
            return extract_secondary(c, buffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<false>*>(work);
            auto details = wptr->details;
            return extract_primary(c, buffer, details.first, details.second, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<true>*>(work);
            auto details = wptr->details;
            return extract_primary(r, dbuffer, ibuffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<true>*>(work);
            auto details = wptr->details;
            return extract_secondary(r, dbuffer, ibuffer, details.first, details.second, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<false>*>(work);
            auto details = wptr->details;
            return extract_secondary(c, dbuffer, ibuffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<false>*>(work);
            auto details = wptr->details;
            return extract_primary(c, dbuffer, ibuffer, details.first, details.second, wptr->base);
        }
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct HDF5PrimarySparseIndexWorkspace : public IndexWorkspace<IDX, WORKROW> {
        HDF5PrimarySparseIndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}
        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }
        PrimarySparseWorkspaceBase base;
    };

    template<bool WORKROW>
    struct HDF5SecondarySparseIndexWorkspace : public IndexWorkspace<IDX, WORKROW> {
        HDF5SecondarySparseIndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}
        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }
        SparseWorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> i, bool cache = false) const {
        if constexpr(ROW) {
            auto ptr = new HDF5PrimarySparseIndexWorkspace<true>(std::move(i));
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base, cache ? nrows : 0);
            return output;
        } else {
            auto ptr = new HDF5SecondarySparseIndexWorkspace<true>(std::move(i));
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> i, bool cache = false) const {
        if constexpr(ROW) {
            auto ptr = new HDF5SecondarySparseIndexWorkspace<false>(std::move(i));
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5PrimarySparseIndexWorkspace<false>(std::move(i));
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base, cache ? ncols : 0);
            return output;
        }
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<true>*>(work);
            return extract_primary(r, buffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<true>*>(work);
            return extract_secondary(r, buffer, wptr->indices_, wptr->base);
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<false>*>(work);
            return extract_secondary(c, buffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<false>*>(work);
            return extract_primary(c, buffer, wptr->indices_, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<true>*>(work);
            return extract_primary(r, dbuffer, ibuffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<true>*>(work);
            return extract_secondary(r, dbuffer, ibuffer, wptr->indices_, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<false>*>(work);
            return extract_secondary(c, dbuffer, ibuffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<false>*>(work);
            return extract_primary(c, dbuffer, ibuffer, wptr->indices_, wptr->base);
        }
    }

private:
    template<class Function, class Skip>
    void extract_primary(size_t i, Function fill, Skip skip, const std::vector<IDX>& indices, PrimarySparseWorkspaceBase& work) const {
        if (indices.empty()) {
            return;
        }

        size_t primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, work);

        const auto& limits = primary_cache_limits[work.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = work.index_cache.begin() + offset;
        auto dstart = work.data_cache.begin() + offset;
        auto iend = istart + primary_length;

        size_t quick_shift = 0;
        if (indices[0] > *istart) { // Both are guaranteed to valid, from length check above.
            bool do_cache = !work.starts.empty();
            if (do_cache && work.starts[i] != -1) {
                quick_shift = work.starts[i];
            } else {
                quick_shift = std::lower_bound(istart, istart + primary_length, indices[0]) - istart;
                if (do_cache) {
                    work.starts[i] = quick_shift;
                }
            }
        }

        istart += quick_shift;
        dstart += quick_shift;

        for (auto idx : indices) {
            while (istart != iend && *istart < idx) {
                ++istart;
                ++dstart;
            }
            if (istart == iend) {
                break;
            }
            if (*istart == idx) {
                fill(idx, *dstart);
            } else {
                skip();
            }
        }
    }

    const T* extract_primary(size_t i, T* dbuffer, const std::vector<IDX>& indices, PrimarySparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + indices.size(), 0);
        auto original = dbuffer;
        extract_primary(i, 
            [&](IDX, T value) -> void {
                *dbuffer = value;
                ++dbuffer;
            }, 
            [&]() -> void{
                ++dbuffer;
            },
            indices, 
            work
        );
        return original;
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, const std::vector<IDX>& indices, PrimarySparseWorkspaceBase& work) const {
        size_t counter = 0;

        extract_primary(i, 
            [&](IDX pos, T value) -> void {
                dbuffer[counter] = value;
                ibuffer[counter] = pos;
                ++counter;
            }, 
            []() -> void {},
            indices, 
            work
        );

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

private:
    template<class Function, class Skip>
    size_t extract_secondary(size_t i, Function fill, Skip skip, const std::vector<IDX>& indices, SparseWorkspaceBase& work) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        for (auto j : indices) {
            if (!extract_secondary_raw(j, i, fill, work)) {
                skip();
            }
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif
    }

    const T* extract_secondary(size_t i, T* dbuffer, const std::vector<IDX>& indices, SparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + indices.size(), 0);
        auto original = dbuffer;
        extract_secondary(i, 
            [&](IDX pos, T value) -> void {
                *dbuffer = value;
                ++dbuffer;
            }, 
            [&]() -> void {
                ++dbuffer;
            },
            indices, 
            work
        );
        return original;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, const std::vector<IDX>& indices, SparseWorkspaceBase& work) const {
        size_t counter = 0;

        extract_secondary(i, 
            [&](IDX pos, T value) -> void {
                dbuffer[counter] = value;
                ibuffer[counter] = pos;
                ++counter;
            }, 
            []() -> void {},
            indices, 
            work
        );

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }
};

}

#endif
