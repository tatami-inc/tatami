#ifndef TATAMI_HDF5_DENSE_MATRIX_HPP
#define TATAMI_HDF5_DENSE_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../base/Matrix.hpp"

namespace tatami {

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

template<typename T, typename IDX, bool transpose = false>
class HDF5DenseMatrix : public tatami::Matrix<T, IDX> {
    size_t firstdim, seconddim;
    std::string file_name, dataset_name;

    size_t nslots = 0;
    size_t cache_size_firstdim = 0;
    size_t cache_size_seconddim = 0;
    bool prefer_firstdim = false;

public:
    HDF5DenseMatrix(std::string file, std::string path, size_t cache_limit = 100000000) : 
        file_name(std::move(file)), 
        dataset_name(std::move(path))
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
            const size_t chunk_firstdim = chunk_dims[0];
            const size_t chunk_seconddim = chunk_dims[1];

            // Need to use the other dimension to figure out
            // the number of chunks along one dimension.
            num_chunks_per_firstdim = std::ceil(static_cast<double>(seconddim)/chunk_seconddim); 
            num_chunks_per_seconddim = std::ceil(static_cast<double>(firstdim)/chunk_firstdim);

            size_t chunk_size = dhandle.getDataType().getSize() * chunk_firstdim * chunk_seconddim;
            cache_size_firstdim = std::min(cache_limit, num_chunks_per_firstdim * chunk_size);
            cache_size_seconddim = std::min(cache_limit, num_chunks_per_seconddim * chunk_size);

            // We used to be able to compute the nslots exactly, but then we got:
            // https://forum.hdfgroup.org/t/unintended-behaviour-for-hash-values-during-chunk-caching/4869/5
            // ... so we'll just set the number of slots to the total number of chunks and forget about it.
            nslots = num_chunks_per_firstdim * num_chunks_per_seconddim; 
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

    bool prefer_rows() const {
        if constexpr(transpose) {
            return !prefer_firstdim;
        } else {
            return prefer_firstdim;
        }
    }

public:
    struct HDF5DenseWorkspace : public Workspace {
        HDF5DenseWorkspace(const std::string& fpath, const std::string& dpath, size_t nslots, size_t cache_size, size_t dim) : 
            current_size(dim),
            memspace(1, &current_size)
        { 
            // Configuring the chunk cache.
            if (nslots) {
                H5::FileAccPropList plist(H5::FileAccPropList::DEFAULT.getId());

                /* The first argument is ignored, according to https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html.
                 * Setting w0 to 0 to evict the last used chunk; no need to worry about full vs partial reads here.
                 */
                plist.setCache(10000, nslots, cache_size, 0);

                file.openFile(fpath, H5F_ACC_RDONLY, plist);
            } else {
                file.openFile(fpath, H5F_ACC_RDONLY);
            }

            dataset = file.openDataSet(dpath);
            dataspace = dataset.getSpace();
            dataspace.selectNone();
            memspace.selectAll();
            return;
        }

        H5::H5File file;
        H5::DataSet dataset;

        H5::DataSpace dataspace;
        hsize_t offset[2];
        hsize_t count[2];

        hsize_t current_size;
        H5::DataSpace memspace;
    };

    std::shared_ptr<Workspace> new_workspace(bool row) const {
        return std::shared_ptr<Workspace>(
            new HDF5DenseWorkspace(
                file_name, 
                dataset_name, 
                nslots, 
                row != transpose ? cache_size_firstdim : cache_size_seconddim,
                row ? ncol() : nrow()
            )
        );
    }

private:
    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, HDF5DenseWorkspace& work) const {
        if constexpr(row != transpose) {
            work.offset[0] = i;
            work.offset[1] = first;
            work.count[0] = 1;
            work.count[1] = last - first;
        } else {
            work.offset[0] = first;
            work.offset[1] = i;
            work.count[0] = last - first;
            work.count[1] = 1;
        }

        work.dataspace.selectHyperslab(H5S_SELECT_SET, work.count, work.offset);

        if (static_cast<size_t>(work.current_size) != last - first) {
            work.current_size = last - first;
            work.memspace.setExtentSimple(1, &(work.current_size));
            work.memspace.selectAll();
        }

        work.dataset.read(buffer, HDF5::define_mem_type<T>(), work.memspace, work.dataspace);
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
