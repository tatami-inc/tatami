#ifndef TATAMI_HDF5_DENSE_MATRIX_HPP
#define TATAMI_HDF5_DENSE_MATRIX_HPP

#include "H5Cpp.h"
#include <string>
#include <cstdint>
#include <type_traits>
#include "../base/Matrix.hpp"

namespace tatami {

namespace HDF5 {

template<typename T>
constexpr H5::PredType& define_type() {
    if constexpr(std::is_same<int, T>::value) {
        return H5::NATIVE_INT;
    } else if (std::is_same<unsigned int, T>::value) {
        return H5::NATIVE_UINT;
    } else if (std::is_same<long, T>::value) {
        return H5::NATIVE_LONG;
    } else if (std::is_same<unsigned long, T>::value) {
        return H5::NATIVE_ULONG;
    } else if (std::is_same<long long, T>::value) {
        return H5::NATIVE_LLONG;
    } else if (std::is_same<unsigned long long, T>::value) {
        return H5::NATIVE_ULLONG;
    } else if (std::is_same<short, T>::value) {
        return H5::NATIVE_SHORT;
    } else if (std::is_same<unsigned short, T>::value) {
        return H5::NATIVE_USHORT;
    } else if (std::is_same<char, T>::value) {
        return H5::NATIVE_CHAR;
    } else if (std::is_same<unsigned char, T>::value) {
        return H5::NATIVE_UCHAR;
    } else if (std::is_same<double, T>::value) {
        return H5::NATIVE_DOUBLE;
    } else if (std::is_same<float, T>::value) {
        return H5::NATIVE_FLOAT;
    } else if (std::is_same<long double, T>::value) {
        return H5::NATIVE_LDOUBLE;
    } else if (std::is_same<uint8_t, T>::value) {
        return H5::NATIVE_UINT8;
    } else if (std::is_same<int8_t, T>::value) {
        return H5::NATIVE_INT8;
    } else if (std::is_same<uint16_t, T>::value) {
        return H5::NATIVE_UINT16;
    } else if (std::is_same<int16_t, T>::value) {
        return H5::NATIVE_INT16;
    } else if (std::is_same<uint32_t, T>::value) {
        return H5::NATIVE_UINT32;
    } else if (std::is_same<int32_t, T>::value) {
        return H5::NATIVE_INT32;
    } else if (std::is_same<uint64_t, T>::value) {
        return H5::NATIVE_UINT64;
    } else if (std::is_same<int64_t, T>::value) {
        return H5::NATIVE_INT64;
    }
    static_assert("unsupported HDF5 type for template parameter 'T'");
}

}

template<typename T, typename IDX, bool transpose = false>
class HDF5DenseMatrix : public tatami::Matrix<T, IDX> {
public:
    struct CacheDetails {
        CacheDetails() : CacheDetails(0, 0) {}
        CacheDetails(size_t n, size_t c) : n_slots(n), cache_size(c) {}
        size_t n_slots;
        size_t cache_size;
    };

    static CacheDetails compute_best_cache_parameters(const H5::DataSet& dataset, bool row, bool col) {
        auto dparms = dataset.getCreatePlist();

        hsize_t dims[2];
        space.getSimpleExtentDims(dims, NULL);
        const size_t nrows = (transpose ? dims[1] : dims[0]);
        const size_t ncols = (transpose ? dims[0] : dims[1]);

        hsize_t chunk_dims[2];
        dparms.getChunk(2, chunk_dims);
        const size_t chunk_nrows = (transpose ? chunk_dims[1] : chunk_dims[0]);
        const size_t chunk_ncols = (transpose ? chunk_dims[0] : chunk_dims[1]);

        const size_t num_chunks_per_row = std::ceil(static_cast<double>(ncol)/chunk_ncols); // per row needs to divide by column dimensions.
        const size_t num_chunks_per_col = std::ceil(static_cast<double>(nrow)/chunk_nrows);

        // We used to be able to compute the nslots exactly, but then we got:
        // https://forum.hdfgroup.org/t/unintended-behaviour-for-hash-values-during-chunk-caching/4869/5
        // ... so we'll just set the number of slots to the total number of chunks and forget about it.
        const size_t num_slots = num_chunks_per_row * num_chunks_per_col; 
    
        auto chunk_size = dataset.getType().getSize() * chunk_nrows * chunk_ncols;
        const size_t cache_size_row = num_chunks_per_row * chunk_size;
        const size_t cache_size_col = num_chunks_per_col * chunk_size;

        if (row && col) {
            return CacheDetails(num_slots, std::max(cache_size_row, cache_size_col)); 
        } else if (row) {
            return CacheDetails(num_slots, cache_size_row); 
        } else if (column) {
            return CacheDetails(num_slots, cache_size_col); 
        } else {
            return CacheDetails();
        }
    }

public:
    HDF5DenseMatrix(std::string file, std::string path, size_t cache_limit) {
        file.openFile(path, H5F_ACC_RDONLY);
        dhandle = file.openDataSet(path);
        auto space = dhandle.getSpace();

        int ndim = space.getSimpleExtentNdims();
        if (ndim != 2) {
           throw std::runtime_error("data in HDF5 dataset is not a two-dimensional array");
        }

        hsize_t dims_out[2];
        space.getSimpleExtentDims(dims_out, NULL);
        firstdim = dims_out[0];
        seconddim = dims_out[1];

        auto curtype=hdata.getTypeClass();
        if (curtype != H5T_INTEGER || curtype != H5T_FLOAT) { 
            throw std::runtime_error("expected numeric data in the HDF5 dataset");
        }

        return;
    }

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

public:
    struct HDF5DenseWorkspace {
        HDF5DenseWorkspace(H5::DataSpace ds, size_t n) : dataspace(std::move(ds)), current_n(n), memspace(1, &current_n) {
            dataspace.selectNone();
            memspace.selectAll();
            return;
        }

        H5::DataSpace dataspace;
        hsize_t offset[2];
        hsize_t count[2];

        hsize_t current_n;
        H5::DataSpace memspace;
    }

    std::shared_ptr<Workspace> new_workspace(bool row) const {
        return std::shared_ptr<Workspace>(new HDF5DenseWorkspace(dataset.getSpace(), row ? ncol() : nrow()));
    }

private:
    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, H5::DataSpace& dataspace, hsize_t* offset, hsize_t* count, hsize_t& n, H5::DataSpace& memspace) const {
        if constexpr(row != transpose) {
            offset[0] = i;
            offset[1] = first;
            count[0] = 1;
            count[1] = last - first;
        } else {
            offset[1] = i;
            offset[0] = first;
            count[1] = 1;
            count[0] = last - first;
        }

        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset):

        if (static_cast<size_t>(n) != last - first) {
            n = last - first;
            memspace.setExtentSimple(1, &n);
            memspace.selectAll();
        }

        dataset.read(buffer, define_type<T>(), memspace, dataspace);
        return buffer;
    }

    template<bool row>
    const T* extract(size_t i, T* buffer, size_t first, size_t last, Workspace* work) const {
        if (work) {
            auto wptr = dynamic_cast<HDF5DenseWorkspace*>(work);
            return extract(i, buffer, first, last, wptr->dataspace, wptr->offset, wptr->count, wptr->current_n, wptr->memspace);
        } else {
            H5::DataSpace dspace = dataset.getSpace();
            hsize_t offset[2];
            hsize_t count[2];
            hsize_t current_n = last - first;
            H5::DataSpace memspace(1, &current_n);
            return extract(i, buffer, first, last, wptr->dataspace, offset, count, current_n, memspace);
        }
    }


public:
    const T* row(size_t r, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return extract<true>(r, buffer, first, last, work);
    }

    const T* column(size_t c, T* buffer, size_t first, size_t last, Workspace* work=nullptr) const {
        return extract<false>(r, buffer, first, last, work);
    }

private:
    size_t firstdim, seconddim;
    H5::H5File file;
    H5::DataSet dataset;
}

}

