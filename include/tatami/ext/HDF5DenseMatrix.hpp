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
 * We use HDF5's chunk caching to speed up access for consecutive rows and columns.
 * However, cache parameters seem to be ignored if any existing HDF5 file handle is open;
 * callers should close other handles before creating a `HDF5DenseMatrix`.
 *
 * Callers should follow the `prefer_rows()` suggestion when extracting data,
 * as this tries to minimize the number of chunks that need to be read per access request.
 * If they do not, the access pattern on disk may be highly suboptimal, depending on the chunk dimensions.
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
    bool prefer_firstdim = false;

public:
    /**
     * @param file Path to the file.
     * @param name Path to the dataset inside the file.
     * @param cache_limit Limit to the size of the chunk cache, in bytes.
     *
     * The actual size of the chunk cache is automatically chosen to optimize for row or column extraction during a call to `new_workspace()`.
     * As long as this automatic choice is less than `cache_limit`, the former is used.
     * The latter is only used if the automatic choice would be too large and we need to constrain memory usage (at the likely cost of some speed).
     */
    HDF5DenseMatrix(std::string file, std::string name, size_t cache_limit = 100000000) : 
        file_name(std::move(file)), 
        dataset_name(std::move(name))
    {
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
        size_t num_chunks_per_firstdim, num_chunks_per_seconddim;

        if (dparms.getLayout() != H5D_CHUNKED) {
            // If contiguous, each firstdim is treated as a chunk.
            // So the number of chunks along the seconddim
            num_chunks_per_firstdim = 1;
            num_chunks_per_seconddim = firstdim;
        } else {
            hsize_t chunk_dims[2];
            dparms.getChunk(2, chunk_dims);
            chunk_firstdim = chunk_dims[0];
            chunk_seconddim = chunk_dims[1];

            // Need to use the other dimension to figure out
            // the number of chunks along one dimension.
            num_chunks_per_firstdim = std::ceil(static_cast<double>(seconddim)/chunk_seconddim); 
            num_chunks_per_seconddim = std::ceil(static_cast<double>(firstdim)/chunk_firstdim);
        }

        prefer_firstdim = num_chunks_per_firstdim < num_chunks_per_seconddim;
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
     * @return `true` if the row extraction requires loading of fewer chunks compared to column extraction.
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
        hsize_t offset[2];
        hsize_t count[2];

        hsize_t memsize;
        H5::DataSpace memspace;

        std::vector<T> cache, buffer, workspace;
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
        auto ptr = new HDF5DenseWorkspace;
        std::shared_ptr<Workspace> output(ptr);

        ptr->file.openFile(file_name, H5F_ACC_RDONLY);
        ptr->dataset = ptr->file.openDataSet(data_name);
        ptr->dataspace = ptr->dataset.getSpace();
        ptr->buffer.resize(chunk_firstdim * chunk_seconddim);
        return output;
    }

private:
    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, HDF5DenseWorkspace& work) const {
        // Figuring out which chunk the request belongs to.
        hsize_t chunk_mydim, chunk_otherdim;
        if constexpr(row != transpose) {
            chunk_mydim = chunk_firstdim;
            chunk_otherdim = chunk_seconddim;
        } else {
            chunk_mydim = chunk_seconddim;
            chunk_otherdim = chunk_firstdim;
        }

        size_t chunk = i / chunk_mydim;
        size_t index = i % chunk_mydim;

        // If we're in the same chunk for 'i' and the first/last are within range, we re-use the cache.
        bool complete_refresh = false, partial_refresh = false;
        if (chunk != work.cached_chunk) {
            complete_refresh = true;
        } else if (first < work.cached_first || last > work.cached_last) {
            partial_refresh = true;
        }

        if (complete_refresh || partial_refresh) {
            size_t first_chunk = first / chunk_otherdim;
            size_t last_chunk = last ? (last - 1) / chunk_otherdim + 1 : 0; // need the chunk after the chunk containing the last included line.

            hsize_t new_cached_first = first_chunk * chunk_otherdim;
            hsize_t new_cached_last = last_chunk * chunk_otherdim;
            if (partial_refresh) {
                // Always expanding the cache for the current 'chunk'. This
                // ensures that we are (eventually) performant for multiple
                // queries to the same 'i' but different 'first' and 'last'.
                new_cached_first = std::min(work.cached_first, new_cached_first);
                new_cached_last = std::max(work.cached_last, new_cached_last);
            }
            size_t new_span = new_cached_last - new_cached_last;

            // Copying the pieces we already have in memory. This uses a 
            // separate workspace as I don't think I can modify the cache in-place.
            T* destination;
            if (partial_refresh) {
                work.workspace.resize(chunk_mydim * new_span);
                size_t src_span = work.cached_last - work.cached_start;
                auto dest = work.workspace.begin() + (new_cached_first - work.cached_start);
                for (size_t c = 0; c < chunk_mydim; ++c, src += src_span, dest += new_span) {
                    std::copy(src, src + src_span, dest);
                }
                destination = work.workspace.data();
            } else {
                work.cache.resize(chunk_mydim * new_span);
                destination = work.cache.data();
            }

            for (size_t c = first_chunk; c < last_chunk; ++c) {
                if (partial_refresh) {
                    size_t position = c * chunk_otherdim;
                    if (position >= work.cached_first && position < work.cached_last) {
                        continue;
                    }
                }

                // For chunks not present in the cache, we need to loop through and read them in.
                if constexpr(row != transpose) {
                    work.offset[0] = chunk * chunk_firstdim;
                    work.offset[1] = c * chunk_seconddim;
                } else {
                    work.offset[0] = c * chunk_firstdim;
                    work.offset[1] = chunk * chunk_seconddim;
                }

                work.count[0] = std::min(work.offset[0] + chunk_firstdim, static_cast<hsize_t>(firstdim)) - work.offset[0];
                work.count[1] = std::min(work.offset[1] + chunk_seconddim, static_cast<hsize_t>(seconddim)) - work.offset[1];
                work.dataspace.selectHyperslab(H5S_SELECT_SET, work.count, work.offset);

                work.memsize = work.count[0] * work.count[1];
                work.memspace.setExtentSimple(1, &work.memsize);
                work.memspace.selectAll();

                work.dataset.read(work.buffer.data(), HDF5::define_mem_type<T>(), work.memspace, work.dataspace);

                // Transferring it to our internal buffer, with any relevant transposition.
                if constexpr(row != transpose) {
                    auto in = work.buffer.begin();
                    for (hsize_t x = 0; x < work.count[0]; ++x, in += work.count[1]) {
                        std::copy(in, in + work.count[1], destination + span * x);
                    }
                } else {
                    for (hsize_t x = 0; x < work.count[1]; ++x) {
                        auto in = work.buffer.begin() + x;
                        auto out = destination + x * span;
                        for (hsize_t y = 0; y < work.count[0]; ++y, in += work.count[1]) {
                            *(out + y) = *in;
                        }
                    }
                }
            }

            if (partial_refresh) {
                work.cache.swap(work.cache.workspace);
            }
            work.cached_first = new_cached_first;
            work.cached_last = new_cached_last;
            work.cached_chunk = chunk;
        }

        auto start = work.cache.begin() + index * (work.cache_last - work.cache_first);
        std::copy(
            start + (first - work.cache_first), 
            start + (last - work.cache_first), 
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
            auto full = new_workspace(row);
            auto wptr = dynamic_cast<HDF5DenseWorkspace*>(full.get());
            return extract<row>(i, buffer, first, last, *wptr);
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
