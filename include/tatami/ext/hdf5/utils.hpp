#ifndef TATAMI_HDF5_UTILS_HPP
#define TATAMI_HDF5_UTILS_HPP

#include "H5Cpp.h"
#include <cstdint>
#include <array>
#include <string>
#include <type_traits>
#include <stdexcept>

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

template<bool integer_only, class GroupLike>
H5::DataSet open_and_check_dataset(const GroupLike& handle, const std::string& name) {
    // Avoid throwing H5 exceptions.
    if (!H5Lexists(handle.getId(), name.c_str(), H5P_DEFAULT) || handle.childObjType(name) != H5O_TYPE_DATASET) {
        throw std::runtime_error("no child dataset named '" + name + "'");
    }

    auto dhandle = handle.openDataSet(name);
    auto type = dhandle.getTypeClass();
    if constexpr(integer_only) {
        if (type != H5T_INTEGER) {
            throw std::runtime_error(std::string("expected integer values in the '") + name + "' dataset");
        }
    } else {
        if (type != H5T_INTEGER && type != H5T_FLOAT) { 
            throw std::runtime_error(std::string("expected numeric values in the '") + name + "' dataset");
        }
    }

    return dhandle;
}

template<int N>
std::array<hsize_t, N> get_array_dimensions(const H5::DataSet& handle, const std::string& name) {
    auto space = handle.getSpace();

    int ndim = space.getSimpleExtentNdims();
    if (ndim != N) {
        throw std::runtime_error(std::string("'") + name + "' should be a " + std::to_string(N) + "-dimensional array");
    }

    std::array<hsize_t, N> dims_out;
    space.getSimpleExtentDims(dims_out.data(), NULL);
    return dims_out;
}

}
/**
 * @endcond
 */

}

#endif
