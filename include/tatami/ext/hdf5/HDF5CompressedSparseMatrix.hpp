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
        primary_cache_id.resize(pointers.size() - 1);
        if (cache_id.size()) {
            primary_cache_limits.emplace_back(0, pointers[1]);
        }

        size_t effective_cache_limit = cache_limit / (sizeof(T) + sizeof(IDX));
        size_t counter = 0, start = 0;
        for (size_t i = 1; i < cache_id.size(); ++i) {
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

        return output;
    }

    struct PrimarySparseWorkspaceBase : public SparseWorkspaceBase {
        std::vector<T> data_cache;
        std::vector<IDX> index_cache;
        size_t current_cache_id;
        bool init = false;
    };

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct HDF5PrimarySparseWorkspace : public Workspace<ROW> {
        PrimarySparseWorkspaceBase base;
    };

    template<bool ROW>
    struct HDF5SecondarySparseWorkspace : public Workspace<ROW> {
        SparseWorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        if constexpr(ROW) {
            auto ptr = new HDF5PrimarySparseWorkspace<ROW>;
            std::shared_ptr<RowWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5SecondarySparseWorkspace<ROW>;
            std::shared_ptr<RowWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    std::shared_ptr<ColumnWorkspace> new_row_workspace() const {
        if constexpr(ROW) {
            auto ptr = new HDF5SecondarySparseWorkspace<ROW>;
            std::shared_ptr<ColumnWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5PrimarySparseWorkspace<ROW>;
            std::shared_ptr<ColumnWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(ROW) {
            return extract_primary(r, buffer, 0, ncols, static_cast<HDF5PrimarySparseWorkspace<ROW>*>(work)->base);
        } else {
            return extract_secondary(r, buffer, 0, ncols, static_cast<HDF5SecondarySparseWorkspace<ROW>*>(work)->base);
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(ROW) {
            return extract_secondary(c, buffer, 0, nrows, static_cast<HDF5SecondarySparseWorkspace<ROW>*>(work)->base);
        } else {
            return extract_primary(c, buffer, 0, nrows, static_cast<HDF5PrimarySparseWorkspace<ROW>*>(work)->base);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, RowWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_primary(r, dbuffer, ibuffer, 0, ncols, static_cast<HDF5PrimarySparseWorkspace<ROW>*>(work)->base);
        } else {
            return extract_secondary(r, dbuffer, ibuffer, 0, ncols, static_cast<HDF5SecondarySparseWorkspace<ROW>*>(work)->base);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_secondary(c, dbuffer, ibuffer, 0, nrows, static_cast<HDF5PrimarySparseWorkspace<ROW>*>(work)->base);
        } else {
            return extract_primary(c, dbuffer, ibuffer, 0, nrows, static_cast<HDF5SecondarySparseWorkspace<ROW>*>(work)->base);
        }
    }

private:
    void populate_primary_cache(size_t i, PrimarySparseWorkspaceBase& work) const {
        if (cache_id[i] == work.current_cache_id && work.init) {
            return;
        }

        work.init = true;

        // Pulling out the entire chunk of primary dimensions containing
        // 'i'. We have to do this for all indices, regardless of the
        // slicing/indexing. 
        work.current_cache_id = cache_id[i];
        const auto& limits = cache_limits[work.current_cache_id];
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

        // In theory, we could avoid extracting entire extent when populating
        // the cache for sliced or indexed queries. In practice, each chunk is
        // often larger than the number of non-zeroes in a column, so we end up
        // having to pull out the entire column anyway. Also, I want to get out
        // of this critical region ASAP and do any indexing elsewhere. 

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

        const auto& limits = cache_limits[work.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = work.index_cache.begin() + offset;
        auto dstart = work.data_cache.begin() + offset;

        size_t request_start = (start > *istart ? std::lower_bound(istart, istart + primary_len, start) - istart : 0);
        istart += request_start;
        dstart += request_start;

        size_t end = start + length;
        for (size_t i = request_start; i < primary_length && *istart < end; ++i, ++istart, ++dstart) {
            fill(*istart, *dstart);
        }

        return;
    }

    const T* extract_primary(size_t i, T* dbuffer, const ExtractType& start, size_t length, PrimarySparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + length, 0);
        extract_primary(i, [&](IDX pos, T value) -> void {
            dbuffer[pos] = value;
        }, start, length, work);
        return dbuffer;
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, const ExtractType& start, size_t length, PrimarySparseWorkspaceBase& work) const {
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
    size_t extract_secondary_raw(IDX primary, size_t secondary, Function& fill, SparseWorkspaceBase& work) const {
        size_t left = pointers[primary], right = pointers[primary + 1];
        index_cache.resize(right - left);

        // Serial locks should be applied by the callers.
        hsize_t offset = left;
        hsize_t count = index_cache.size();
        work.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        work.memspace.setExtentSimple(1, &count);
        work.memspace.selectAll();
        work.index.read(index_cache.data(), HDF5::define_mem_type<IDX>(), work.memspace, work.dataspace);

        auto it = std::lower_bound(index_cache.begin(), index_cache.end(), secondary);
        if (it != index_cache.end() && *it == secondary) {
            offset = left + (it - index_cache.begin());
            count = 1;
            work.dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
            work.memspace.setExtentSimple(1, &count);
            work.memspace.selectAll();

            T dest;
            work.data.read(&dest, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
            fill(primary, dest);
        }
    }

    template<class Function, class ExtractType>
    size_t extract_secondary(size_t i, Function fill, size_t start, size_t length, SparseWorkspaceBase& work) const {
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

        return counter;
    }

    const T* extract_secondary(size_t i, T* dbuffer, const ExtractType& start, size_t length, SparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + length, 0);
        extract_secondary(i, [&](IDX pos, T value) -> void {
            dbuffer[pos] = value;
        }, start, length, work);
        return dbuffer;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, const ExtractType& start, size_t length, SparseWorkspaceBase& work) const {
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
    template<bool ROW>
    struct HDF5PrimarySparseBlockWorkspace : public BlockWorkspace<ROW> {
        HDF5PrimarySparseBlockWorkspace(size_t s, size_t l) : details(s, l) {}
        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }
        PrimarySparseWorkspaceBase base;
    };

    template<bool ROW>
    struct HDF5SecondarySparseBlockWorkspace : public BlockWorkspace<ROW> {
        HDF5SecondarySparseBlockWorkspace(size_t s, size_t l) : details(s, l) {}
        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }
        SparseWorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t s, size_t l) const {
        if constexpr(ROW) {
            auto ptr = new HDF5PrimarySparseBlockWorkspace<ROW>(s, l);
            std::shared_ptr<RowBlockWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5SecondarySparseBlockWorkspace<ROW>(s, l);
            std::shared_ptr<RowBlockWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_row_workspace() const {
        if constexpr(ROW) {
            auto ptr = new HDF5SecondarySparseBlockWorkspace<ROW>(s, l);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5PrimarySparseBlockWorkspace<ROW>(s, l);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_primary(r, buffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_secondary(r, buffer, details.first, details.second, wptr->base);
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_secondary(c, buffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_primary(c, buffer, details.first, details.second, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_primary(r, dbuffer, ibuffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_secondary(r, dbuffer, ibuffer, details.first, details.second, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_secondary(c, dbuffer, ibuffer, details.first, details.second, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseBlockWorkspace<ROW>*>(work);
            auto details = wptr->details;
            return extract_primary(c, dbuffer, ibuffer, details.first, details.second, wptr->base);
        }
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct HDF5PrimarySparseIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        HDF5PrimarySparseIndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}
        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }
        PrimarySparseWorkspaceBase base;
    };

    template<bool ROW>
    struct HDF5SecondarySparseIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        HDF5SecondarySparseIndexWorkspace(std::vector<IDX> i) : indices_(std::move(i)) {}
        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }
        SparseWorkspaceBase base;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(size_t s, size_t l) const {
        if constexpr(ROW) {
            auto ptr = new HDF5PrimarySparseIndexWorkspace<ROW>(s, l);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5SecondarySparseIndexWorkspace<ROW>(s, l);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_row_workspace() const {
        if constexpr(ROW) {
            auto ptr = new HDF5SecondarySparseIndexWorkspace<ROW>(s, l);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base);
            return output;
        } else {
            auto ptr = new HDF5PrimarySparseIndexWorkspace<ROW>(s, l);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);
            fill_base(ptr->base);
            return output;
        }
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<ROW>*>(work);
            return extract_primary(r, buffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<ROW>*>(work);
            return extract_secondary(r, buffer, wptr->indices_, wptr->base);
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<ROW>*>(work);
            return extract_secondary(c, buffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<ROW>*>(work);
            return extract_primary(c, buffer, wptr->indices_, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<ROW>*>(work);
            return extract_primary(r, dbuffer, ibuffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<ROW>*>(work);
            return extract_secondary(r, dbuffer, ibuffer, wptr->indices_, wptr->base);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(ROW) {
            auto wptr = static_cast<HDF5SecondarySparseIndexWorkspace<ROW>*>(work);
            return extract_secondary(c, dbuffer, ibuffer, wptr->indices_, wptr->base);
        } else {
            auto wptr = static_cast<HDF5PrimarySparseIndexWorkspace<ROW>*>(work);
            return extract_primary(c, dbuffer, ibuffer, wptr->indices_, wptr->base);
        }
    }

private:
    template<class Function>
    void extract_primary(size_t i, Function fill, const std::vector<IDX>& indices, PrimarySparseWorkspaceBase& work) const {
        if (indices.empty()) {
            return;
        }

        size_t primary_length = pointers[i + 1] - pointers[i];
        if (primary_length ==0) {
            return;
        }

        populate_primary_cache(i, work);

        const auto& limits = cache_limits[work.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = work.index_cache.begin() + offset;
        auto dstart = work.data_cache.begin() + offset;
        auto iend = istart + primary_length;

        // Guaranteed to be of non-zero length from check above.
        size_t quick_shift = (indices[0] > *istart ? std::lower_bound(istart, istart + primary_length, indices[0]) - istart : 0);
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
            }
        }
    }

    const T* extract_primary(size_t i, T* dbuffer, const std::vector<IDX>& indices, PrimarySparseWorkspace& work) const {
        std::fill(dbuffer, dbuffer + indices.size(), 0);
        extract_primary(i, [&](IDX pos, T value) -> void {
            dbuffer[pos] = value;
        }, indices, work);
        return dbuffer;
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, const std::vector<IDX>& indices, PrimarySparseWorkspace& work) const {
        size_t counter = 0;

        extract_primary(i, [&](IDX pos, T value) -> void {
            dbuffer[counter] = value;
            ibuffer[counter] = pos;
            ++counter;
        }, indices, work);

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }

private:
    template<class Function, class ExtractType>
    size_t extract_secondary(size_t i, Function fill, const std::vector<IDX>& indices, SparseWorkspaceBase& work) const {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        for (auto j : indices) {
            extract_secondary_raw(j, i, fill, work);
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif

        return counter;
    }

    const T* extract_secondary(size_t i, T* dbuffer, const ExtractType& start, size_t length, SparseWorkspaceBase& work) const {
        std::fill(dbuffer, dbuffer + length, 0);
        extract_secondary(i, [&](IDX pos, T value) -> void {
            dbuffer[pos] = value;
        }, start, length, work);
        return dbuffer;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, const ExtractType& start, size_t length, SparseWorkspaceBase& work) const {
        size_t counter = 0;

        extract_secondary(i, [&](IDX pos, T value) -> void {
            dbuffer[counter] = value;
            ibuffer[counter] = pos;
            ++counter;
        }, start, length, work);

        return SparseRange<T, IDX>(counter, dbuffer, ibuffer);
    }
};

}

#endif
