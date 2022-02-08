#ifndef TATAMI_HDF5_DENSE_MATRIX_HPP
#define TATAMI_HDF5_DENSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../base/Matrix.hpp"

/**
 * @file HDF5DenseMatrix.hpp
 *
 * @brief Defines a class for a HDF5-backed dense matrix.
 */

namespace tatami {

/**
 * @cond
 */

namespace HDF5 {

template<typename T>
const H5::PredType& define_mem_type() {
    if constexpr(std::is_same<int, T>::value) {
        return H5::PredType::NATIVE_INT;
    } else if (std::is_same<unsigned int, T>::value) {
        return H5::PredType::NATIVE_UINT;
    } else if (std::is_same<long, T>::value) {
        return H5::PredType::NATIVE_LONG;
    } else if (std::is_same<unsigned long, T>::value) {
        return H5::PredType::NATIVE_ULONG;
    } else if (std::is_same<long long, T>::value) {
        return H5::PredType::NATIVE_LLONG;
    } else if (std::is_same<unsigned long long, T>::value) {
        return H5::PredType::NATIVE_ULLONG;
    } else if (std::is_same<short, T>::value) {
        return H5::PredType::NATIVE_SHORT;
    } else if (std::is_same<unsigned short, T>::value) {
        return H5::PredType::NATIVE_USHORT;
    } else if (std::is_same<char, T>::value) {
        return H5::PredType::NATIVE_CHAR;
    } else if (std::is_same<unsigned char, T>::value) {
        return H5::PredType::NATIVE_UCHAR;
    } else if (std::is_same<double, T>::value) {
        return H5::PredType::NATIVE_DOUBLE;
    } else if (std::is_same<float, T>::value) {
        return H5::PredType::NATIVE_FLOAT;
    } else if (std::is_same<long double, T>::value) {
        return H5::PredType::NATIVE_LDOUBLE;
    } else if (std::is_same<uint8_t, T>::value) {
        return H5::PredType::NATIVE_UINT8;
    } else if (std::is_same<int8_t, T>::value) {
        return H5::PredType::NATIVE_INT8;
    } else if (std::is_same<uint16_t, T>::value) {
        return H5::PredType::NATIVE_UINT16;
    } else if (std::is_same<int16_t, T>::value) {
        return H5::PredType::NATIVE_INT16;
    } else if (std::is_same<uint32_t, T>::value) {
        return H5::PredType::NATIVE_UINT32;
    } else if (std::is_same<int32_t, T>::value) {
        return H5::PredType::NATIVE_INT32;
    } else if (std::is_same<uint64_t, T>::value) {
        return H5::PredType::NATIVE_UINT64;
    } else if (std::is_same<int64_t, T>::value) {
        return H5::PredType::NATIVE_INT64;
    }
    static_assert("unsupported HDF5 type for template parameter 'T'");
}

}

/**
 * @endcond
 */

/**
 * @brief Dense matrix backed by a DataSet in a HDF5 file.
 *
 * This class retrieves data from the HDF5 file on demand rather than loading it all in at the start.
 * This allows us to handle very large datasets in limited memory at the cost of some speed.
 *
 * We manually handle the chunk caching to speed up access for consecutive rows and columns.
 * The policy is to minimize the number of calls to the HDF5 library by requesting large contiguous slices where possible.
 * The size of the slice is determined by the cache limit in the constructor.
 *
 * Callers should follow the `prefer_rows()` suggestion when extracting data,
 * as this tries to minimize the number of chunks that need to be read per access request.
 * If they do not, the access pattern on disk may be slightly to highly suboptimal, depending on the chunk dimensions.
 *
 * As the HDF5 library is not generally thread-safe, the HDF5-related operations should only be run in a single thread.
 * For OpenMP, this is handled automatically by putting all HDF5 operations in a critical region.
 * For other parallelization schemes, callers should define the `TATAMI_HDF5_PARALLEL_LOCK` macro;
 * this should be a function that accepts and executes a no-argument lambda within an appropriate serial region (e.g., based on a global mutex).
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 * @tparam transpose Whether the dataset is transposed in its storage order, i.e., rows in HDF5 are columns in this matrix.
 */
template<typename T, typename IDX, bool transpose = false>
class HDF5DenseMatrix : public tatami::Matrix<T, IDX> {
    size_t firstdim, seconddim;
    std::string file_name, dataset_name;

    hsize_t chunk_firstdim, chunk_seconddim;
    hsize_t cache_firstdim, cache_seconddim;
    bool prefer_firstdim;

public:
    /**
     * @param file Path to the file.
     * @param name Path to the dataset inside the file.
     * @param cache_limit Limit to the size of the chunk cache, in bytes.
     *
     * The cache size should be large enough to fit all chunks spanned by a row or column, for (near-)consecutive row and column access respectively.
     * Otherwise, performance will degrade as the same chunks may need to be repeatedly read back into memory.
     */
    HDF5DenseMatrix(std::string file, std::string name, size_t cache_limit = 100000000) : 
        file_name(std::move(file)), 
        dataset_name(std::move(name))
    {
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        #pragma omp critical
        {
#else
        TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

        H5::H5File fhandle(file_name, H5F_ACC_RDONLY);
        auto dhandle = fhandle.openDataSet(dataset_name);
        auto space = dhandle.getSpace();

        int ndim = space.getSimpleExtentNdims();
        if (ndim != 2) {
           throw std::runtime_error("data in HDF5 dataset is not a two-dimensional array");
        }

        hsize_t dims_out[2];
        space.getSimpleExtentDims(dims_out, NULL);
        firstdim = dims_out[0];
        seconddim = dims_out[1];

        auto curtype = dhandle.getTypeClass();
        if (curtype != H5T_INTEGER && curtype != H5T_FLOAT) { 
            throw std::runtime_error("expected numeric data in the HDF5 dataset");
        }

        auto dparms = dhandle.getCreatePlist();
        if (dparms.getLayout() != H5D_CHUNKED) {
            // If contiguous, each firstdim is treated as a chunk.
            chunk_firstdim = 1;
            chunk_seconddim = seconddim;
        } else {
            hsize_t chunk_dims[2];
            dparms.getChunk(2, chunk_dims);
            chunk_firstdim = chunk_dims[0];
            chunk_seconddim = chunk_dims[1];
        }

        auto limit_cache = [](size_t chunk_dim, size_t mydim, size_t otherdim, size_t type_size, size_t limit) -> size_t {
            size_t nelements_along_chunkdim = std::min(mydim, limit / (type_size * otherdim));
            size_t nchunks = nelements_along_chunkdim / chunk_dim;
            return std::min(nelements_along_chunkdim, nchunks * chunk_dim); 
        };

        size_t data_size = dhandle.getDataType().getSize();
        cache_firstdim = limit_cache(chunk_firstdim, firstdim, seconddim, data_size, cache_limit);
        cache_seconddim = limit_cache(chunk_seconddim, seconddim, firstdim, data_size, cache_limit);

        // Favoring extraction along the first dimension if it's not crippled
        // by the chunking. This favors pulling of contiguous sections along
        // the last (i.e., second) dimension, which is how HDF5 stores and
        // returns 2D arrays, so we avoid needing some manual transposition.
        if (cache_firstdim >= chunk_firstdim) {
            prefer_firstdim = true;
        } else if (cache_seconddim >= chunk_seconddim) {
            // Of course, if the first dimension is chunk-compromised, then we
            // would rather do the transpositions than repeat a read from file.
            prefer_firstdim = false;
        } else {
            // If both are crippled, then we might as well extract along the first dimension.
            prefer_firstdim = true;
        }

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif

        return;
    }

public:
    size_t nrow() const {
        if constexpr(transpose) {
            return seconddim;
        } else {
            return firstdim;
        }
    }

    size_t ncol() const {
        if constexpr(transpose) {
            return firstdim;
        } else {
            return seconddim;
        }
    }

    /**
     * @return Boolean indicating whether to prefer row extraction.
     *
     * We favor extraction on the first dimension (rows by default, columns when `transpose = true`) as this matches the HDF5 storage order.
     * However, for some chunking scheme and `cache_limit`, this might require repeated reads from file;
     * in such cases, we switch to extraction on the second dimension.
     */
    bool prefer_rows() const {
        if constexpr(transpose) {
            return !prefer_firstdim;
        } else {
            return prefer_firstdim;
        }
    }

public:
    /**
     * @cond
     */
    struct HDF5DenseWorkspace : public Workspace {
        H5::H5File file;
        H5::DataSet dataset;

        H5::DataSpace dataspace;
        H5::DataSpace memspace;

        std::vector<T> cache, buffer;
        size_t cached_chunk = -1;
        hsize_t cached_first = 0, cached_last = 0;
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

        auto ptr = new HDF5DenseWorkspace;
        output.reset(ptr);

        // Turn off HDF5's caching, as we'll be handling that.
        H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
        fapl.setCache(0, 0, 0, 0);

        ptr->file.openFile(file_name, H5F_ACC_RDONLY, fapl);
        ptr->dataset = ptr->file.openDataSet(dataset_name);
        ptr->dataspace = ptr->dataset.getSpace();

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
        }
#else
        });
#endif

        return output;
    }

private:
    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, H5::DataSet& dataset, H5::DataSpace& dataspace, H5::DataSpace& memspace) const {
        hsize_t offset[2];
        hsize_t count[2];

        constexpr int x = (row != transpose);
        offset[1-x] = i;
        offset[x] = first;
        count[1-x] = 1;
        count[x] = last - first;

        // Serial locks are applied in callers.
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
        memspace.setExtentSimple(2, count);
        memspace.selectAll();
        dataset.read(buffer, HDF5::define_mem_type<T>(), memspace, dataspace);

        return buffer;
    }

    /* We manually handle the chunk caching, in particular so that we only
     * pay the cost of rearranging the chunked data into contiguous arrays once.
     */
    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, HDF5DenseWorkspace& work) const {
        // Figuring out which chunk the request belongs to.
        hsize_t cache_mydim, chunk_otherdim, mydim, otherdim;
        if constexpr(row != transpose) {
            cache_mydim = cache_firstdim;
            chunk_otherdim = chunk_seconddim;
            mydim = firstdim;
            otherdim = seconddim;
        } else {
            cache_mydim = cache_seconddim;
            chunk_otherdim = chunk_firstdim;
            mydim = seconddim;
            otherdim = firstdim;
        }

        // No caching can be done here, so we just extract directly.
        if (cache_mydim == 0) {
            const T* out;
#ifndef TATAMI_HDF5_PARALLEL_LOCK        
            #pragma omp critical
            {
#else
            TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

            out = extract<row>(i, buffer, first, last, work.dataset, work.dataspace, work.memspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
            }
#else
            });
#endif
            return out;
        }

        size_t chunk = i / cache_mydim;
        size_t index = i % cache_mydim;

        // If we're in the same chunk for 'i' and the 'first'/'last' are within
        // range, we re-use the cache.  In theory, we could minimize
        // extractions if the 'first' and 'last' are overlapping with their
        // cached counterparts, but it's pretty rare that 'first' and 'last'
        // will change, so we'll just keep things simple.
        if (chunk != work.cached_chunk || first < work.cached_first || last > work.cached_last) {
            size_t first_chunk = first / chunk_otherdim;

            // Need the chunk _after_ the chunk containing the _before-last_ element.
            size_t last_chunk = last ? (last - 1) / chunk_otherdim + 1 : 0; 

            hsize_t new_cached_first = first_chunk * chunk_otherdim;
            hsize_t new_cached_last = std::min(otherdim, last_chunk * chunk_otherdim);
            size_t new_span = new_cached_last - new_cached_first;

            hsize_t cache_mydim_start = chunk * cache_mydim;
            hsize_t cache_mydim_end = std::min(mydim, cache_mydim_start + cache_mydim);
            hsize_t cache_mydim_actual = cache_mydim_end - cache_mydim_start;
            size_t new_cache_size = new_span * cache_mydim;

            T* destination;
            work.cache.resize(new_cache_size);
            if constexpr(row != transpose) {
                destination = work.cache.data();
            } else {
                work.buffer.resize(new_cache_size);
                destination = work.buffer.data();
            }

            // For chunks not present in the cache, we need to loop through and read them in.
            hsize_t offset[2];
            hsize_t count[2];
            {
                constexpr int x = (row != transpose);
                offset[1-x] = cache_mydim_start;
                offset[x] = new_cached_first;
                count[1-x] = cache_mydim_actual;
                count[x] = new_span;
            }

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
            #pragma omp critical
            {
#else
            TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif

            work.dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

            // HDF5 is a lot faster when the memspace and dataspace match in dimensionality.
            // Presumably there is some shuffling that happens inside when dimensions don't match.
            work.memspace.setExtentSimple(2, count);
            work.memspace.selectAll();

            work.dataset.read(destination, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
            }
#else
            });
#endif

            if constexpr(row == transpose) {
                auto output = work.cache.begin();
                for (hsize_t x = 0; x < count[1]; ++x, output += new_span) {
                    auto in = work.buffer.begin() + x;
                    for (hsize_t y = 0; y < count[0]; ++y, in += count[1]) {
                        *(output + y) = *in;
                    }
                }
            }

            work.cached_first = new_cached_first;
            work.cached_last = new_cached_last;
            work.cached_chunk = chunk;
        }

        auto start = work.cache.begin() + index * (work.cached_last - work.cached_first);
        std::copy(
            start + (first - work.cached_first), 
            start + (last - work.cached_first), 
            buffer
        );
        return buffer;
    }

    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5DenseWorkspace*>(work);
            return extract<row>(i, buffer, first, last, *wptr);
        } else {
            const T* out;

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
            #pragma omp critical
            {
#else
            TATAMI_HDF5_PARALLEL_LOCK([&]() -> void {
#endif
            // Bypass all caching, manual and HDF5.
            H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
            fapl.setCache(10000, 0, 0, 0);

            H5::H5File file(file_name, H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT, fapl);
            auto dataset = file.openDataSet(dataset_name);
            auto dataspace = dataset.getSpace();

            H5::DataSpace memspace;
            out = extract<row>(i, buffer, first, last, dataset, dataspace, memspace);

#ifndef TATAMI_HDF5_PARALLEL_LOCK        
            }
#else
            });
#endif
            return out;
        }
    }

public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return extract<true>(r, buffer, first, last, work);
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return extract<false>(c, buffer, first, last, work);
    }

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::column;
};

}

#endif
