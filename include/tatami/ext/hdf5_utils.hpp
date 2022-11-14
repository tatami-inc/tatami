#ifndef TATAMI_HDF5_UTILS_HPP
#define TATAMI_HDF5_UTILS_HPP

#include "H5Cpp.h"
#include <cstdint>

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

template<bool integer_only>
size_t get_1d_array_length(const H5::DataSet& handle, const std::string& name) {
    auto type = handle.getTypeClass();

    if constexpr(integer_only) {
        if (type != H5T_INTEGER) {
            throw std::runtime_error(std::string("expected integer values in the '") + name + "' dataset");
        }
    } else {
        if (type != H5T_INTEGER && type != H5T_FLOAT) { 
            throw std::runtime_error(std::string("expected numeric values in the '") + name + "' dataset");
        }
    }

    auto space = handle.getSpace();
    int ndim = space.getSimpleExtentNdims();
    if (ndim != 1) {
        throw std::runtime_error(std::string("'") + name + "' should be a one-dimensional array");
    }

    hsize_t dims_out[1];
    space.getSimpleExtentDims(dims_out, NULL);
    return dims_out[0];
}

inline std::pair<size_t, size_t> get_2d_array_dims(const H5::DataSet& handle, const std::string& name) {
    auto space = handle.getSpace();

    int ndim = space.getSimpleExtentNdims();
    if (ndim != 2) {
       throw std::runtime_error("'" + name + "' dataset is not a two-dimensional array");
    }

    hsize_t dims_out[2];
    space.getSimpleExtentDims(dims_out, NULL);

    auto curtype = handle.getTypeClass();
    if (curtype != H5T_INTEGER && curtype != H5T_FLOAT) { 
        throw std::runtime_error("expected numeric data in the '" + name + "' dataset");
    }

    return std::pair<size_t, size_t>(dims_out[0], dims_out[1]);
}

}
/**
 * @endcond
 */

}

#endif
