#ifndef TATAMI_HDF5_SPARSE_MATRIX_HPP
#define TATAMI_HDF5_SPARSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../base/Matrix.hpp"
#include "HDF5DenseMatrix.hpp"

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

    std::vector<size_t> cache_id;
    std::vector<std::pair<size_t, size_t> > cache_limits;

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
     * @param cache_limit Limit to the size of the chunk cache, in bytes.
     *
     * The cache is created by extracting multiple columns (for CSC matrices) or rows (CSR) on every call to the HDF5 library.
     * These are held in memory in the workspace created by `new_workspace()`, while the relevant column/row is returned to the user by `row()` or `column()`.
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

        auto check_size = [&](const H5::DataSet& handle, const std::string& name) -> size_t {
            auto type = handle.getTypeClass();
            if (type != H5T_INTEGER && type != H5T_FLOAT) { 
                throw std::runtime_error(std::string("expected numeric values in the '") + name + "' dataset");
            }
            auto space = handle.getSpace();
            
            int ndim = space.getSimpleExtentNdims();
            if (ndim != 1) {
                throw std::runtime_error(std::string("'") + name + "' should be a one-dimensional array");
            }

            hsize_t dims_out[1];
            space.getSimpleExtentDims(dims_out, NULL);
            return dims_out[0];
        };

        auto dhandle = file_handle.openDataSet(data_name);
        const size_t nonzeros = check_size(dhandle, "vals");

        auto ihandle = file_handle.openDataSet(index_name);
        if (check_size(ihandle, "idx") != nonzeros) {
            throw std::runtime_error("number of non-zero elements is not consistent between 'data' and 'idx'");
        }

        total_element_size = dhandle.getDataType().getSize() + ihandle.getDataType().getSize();

        auto phandle = file_handle.openDataSet(ptr);
        const size_t ptr_size = check_size(phandle, "ptr");
        if (ptr_size != pointers.size()) {
            throw std::runtime_error("'ptr' is not of the appropriate length");
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

        size_t maxgap = 0;
        for (size_t i = 1; i < pointers.size(); ++i) {
            if (pointers[i] < pointers[i-1]) {
                throw std::runtime_error("index pointers should be sorted");
            }
            size_t delta = pointers[i] - pointers[i - 1];
            if (delta > maxgap) {
                maxgap = delta;
            }
        }

        // Setting up the cache parameters.
        cache_id.resize(pointers.size() - 1);
        if (cache_id.size()) {
            cache_limits.emplace_back(0, pointers[1]);
        }

        size_t effective_cache_limit = cache_limit / total_element_size;

        size_t counter = 0, start = 0;
        for (size_t i = 1; i < cache_id.size(); ++i) {
            size_t end = pointers[i + 1];
            if (end - start <= effective_cache_limit) {
                cache_limits.back().second = end;
            } else {
                ++counter;
                start = pointers[i];
                cache_limits.emplace_back(start, end);
            }
            cache_id[i] = counter;
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

public:
    /**
     * @cond
     */
    struct HDF5SparseWorkspace : public Workspace {
        H5::H5File file;
        H5::DataSet data, index;
        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        std::vector<T> data_cache;
        std::vector<IDX> index_cache;
        size_t current_cache_id = -1;
    };
    /**
     * @endcond
     */

    /**
     * @param row Should a workspace be created for row-wise extraction?
     * @return A shared pointer to a `Workspace` object is returned.
     */
    std::shared_ptr<Workspace> new_workspace(bool row) const {
        std::shared_ptr<Workspace> output;

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        auto ptr = new HDF5SparseWorkspace;
        output.reset(ptr);

        H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
        fapl.setCache(0, 0, 0, 0);
        ptr->file.openFile(file_name, H5F_ACC_RDONLY, fapl);

        ptr->data = ptr->file.openDataSet(data_name);
        ptr->index = ptr->file.openDataSet(index_name);
        ptr->dataspace = ptr->data.getSpace();

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif

        return output;
    }

private:
    template<typename IndexIt, typename DataIt, typename Thing> 
    size_t copy_primary_to_buffer(IndexIt istart, DataIt dstart, size_t len, size_t first, size_t last, T* dbuffer, Thing thing) const {
        size_t request_start = (first > *istart ? std::lower_bound(istart, istart + len, first) - istart : 0);
        size_t request_end = (last < *(istart + len - 1) ? std::lower_bound(istart + request_start, istart + len, last) - istart : len);

        if constexpr(std::is_same<typename std::remove_reference<Thing>::type, IDX*>::value) {
            std::copy(istart + request_start, istart + request_end, thing);
            std::copy(dstart + request_start, dstart + request_end, dbuffer);
        } else {
            std::fill(dbuffer, dbuffer + (last - first), 0);
            for (size_t i = request_start; i < request_end; ++i) {
                dbuffer[*(istart + i) - first] = *(dstart + i);
            }
        }

        return request_end - request_start; // number of non-zeros.
    }

    template<typename Thing>
    size_t extract_primary(size_t i, T* dbuffer, Thing thing, size_t first, size_t last, HDF5SparseWorkspace& work) const {
        if (cache_id[i] != work.current_cache_id) {
            // Pulling out the entire chunk containing 'i'.
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

            work.data_cache.resize(count);
            work.data.read(work.data_cache.data(), HDF5::define_mem_type<T>(), work.memspace, work.dataspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
            }
#else
            });
#endif
        }

        size_t len = pointers[i + 1] - pointers[i];
        if (len == 0) {
            return 0;
        }

        const auto& limits = cache_limits[work.current_cache_id];
        size_t offset = pointers[i] - limits.first;
        auto istart = work.index_cache.begin() + offset;
        auto dstart = work.data_cache.begin() + offset;
        return copy_primary_to_buffer(istart, dstart, len, first, last, dbuffer, thing);
    }

    template<typename Thing>
    size_t extract_primary(size_t i, T* dbuffer, Thing thing, size_t first, size_t last) const {
        // Ignore all caches.
        hsize_t offset = pointers[i];
        hsize_t count = pointers[i + 1] - pointers[i];
        if (count == 0) {
            return 0;
        }

        std::vector<IDX> index_cache(count);
        std::vector<T> data_cache(count);

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
        fapl.setCache(0, 0, 0, 0);
        H5::H5File file(file_name, H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT, fapl);

        auto data = file.openDataSet(data_name);
        auto index = file.openDataSet(index_name);
        auto dataspace = data.getSpace();

        dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        H5::DataSpace memspace(1, &count);
        memspace.selectAll();

        index.read(index_cache.data(), HDF5::define_mem_type<IDX>(), memspace, dataspace);
        data.read(data_cache.data(), HDF5::define_mem_type<T>(), memspace, dataspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif

        return copy_primary_to_buffer(index_cache.begin(), data_cache.begin(), count, first, last, dbuffer, thing);
    }

    SparseRange<T, IDX> extract_primary(size_t i, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            return SparseRange<T, IDX>(extract_primary(i, dbuffer, ibuffer, first, last, *wptr), dbuffer, ibuffer);
        } else {
            return SparseRange<T, IDX>(extract_primary(i, dbuffer, ibuffer, first, last), dbuffer, ibuffer);
        }
    }

private:
    template<typename Thing>
    size_t extract_secondary(size_t i, T* dbuffer, Thing thing, size_t first, size_t last, 
        H5::DataSet& data, H5::DataSet& index, 
        H5::DataSpace& dataspace, H5::DataSpace& memspace, 
        std::vector<IDX>& index_cache) const
    {
        // Looping over all secondary dimensions.
        size_t counter = 0;
        for (size_t j = first; j < last; ++j) {
            size_t left = pointers[j], right = pointers[j + 1];
            index_cache.resize(right - left);

            // Serial locks should be applied by the callers.
            hsize_t offset = left;
            hsize_t count = index_cache.size();
            dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
            memspace.setExtentSimple(1, &count);
            memspace.selectAll();
            index.read(index_cache.data(), HDF5::define_mem_type<IDX>(), memspace, dataspace);

            auto it = std::lower_bound(index_cache.begin(), index_cache.end(), i);
            if (it != index_cache.end() && *it == i) {
                offset = left + (it - index_cache.begin());
                count = 1;
                dataspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
                memspace.setExtentSimple(1, &count);
                memspace.selectAll();

                auto dest = dbuffer;
                if constexpr(std::is_same<Thing, IDX*>::value) {
                    dest += counter;
                } else {
                    dest += j - first;
                }
                data.read(dest, HDF5::define_mem_type<T>(), memspace, dataspace);

                if constexpr(std::is_same<typename std::remove_reference<Thing>::type, IDX*>::value) {
                    thing[counter] = j;
                    ++counter;
                }
            }
        }

        return counter;
    }

    template<typename Thing>
    size_t extract_secondary(size_t i, T* dbuffer, Thing thing, size_t first, size_t last, HDF5SparseWorkspace& work) const {
        size_t n;

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        // Reusing the index cache inside the workspace as a holding ground for the indices extracted for each column.
        n = extract_secondary(i, dbuffer, thing, first, last, work.data, work.index, work.dataspace, work.memspace, work.index_cache);

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif

        return n;
    }

    template<typename Thing>
    size_t extract_secondary(size_t i, T* dbuffer, Thing thing, size_t first, size_t last) const {
        size_t n;
        
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif
        
        // Ignore all caches, not that it really makes a difference here anyway.
        H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
        fapl.setCache(0, 0, 0, 0);
        H5::H5File file(file_name, H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT, fapl);

        auto data = file.openDataSet(data_name);
        auto index = file.openDataSet(index_name);
        auto dataspace = data.getSpace();

        H5::DataSpace memspace;
        std::vector<IDX> index_cache;
        n = extract_secondary(i, dbuffer, thing, first, last, data, index, dataspace, memspace, index_cache);

#ifndef TATAMI_HDF5_PARALLEL_LOCK 
        }
#else
        });
#endif

        return n;
    }

    SparseRange<T, IDX> extract_secondary(size_t i, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            return SparseRange<T, IDX>(extract_secondary(i, dbuffer, ibuffer, first, last, *wptr), dbuffer, ibuffer);
        } else {
            return SparseRange<T, IDX>(extract_secondary(i, dbuffer, ibuffer, first, last), dbuffer, ibuffer);
        }
    }

public:
    SparseRange<T, IDX> sparse_row(size_t r, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_primary(r, dbuffer, ibuffer, first, last, work);
        } else {
            return extract_secondary(r, dbuffer, ibuffer, first, last, work);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* dbuffer, IDX* ibuffer, size_t first, size_t last, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(ROW) {
            return extract_secondary(c, dbuffer, ibuffer, first, last, work);
        } else {
            return extract_primary(c, dbuffer, ibuffer, first, last, work);
        }
    }

    using Matrix<T, IDX>::sparse_row;

    using Matrix<T, IDX>::sparse_column;

private:
    const T* expand_primary(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            extract_primary(i, buffer, false, first, last, *wptr);
        } else {
            extract_primary(i, buffer, false, first, last);
        }
        return buffer;
    }

    const T* expand_secondary(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5SparseWorkspace*>(work);
            extract_secondary(i, buffer, false, first, last, *wptr);
        } else {
            extract_secondary(i, buffer, false, first, last);
        }
        return buffer;
    }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            return expand_primary(r, buffer, first, last, work);
        } else {
            return expand_secondary(r, buffer, first, last, work);
        }
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        if constexpr(ROW) {
            return expand_secondary(c, buffer, first, last, work);
        } else {
            return expand_primary(c, buffer, first, last, work);
        }
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;
};

}

#endif
