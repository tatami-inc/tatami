#ifndef TATAMI_WRITE_SPARSE_MATRIX_TO_HDF5_HPP
#define TATAMI_WRITE_SPARSE_MATRIX_TO_HDF5_HPP

#include "../../base/Matrix.hpp"
#include "utils.hpp"

#include "H5Cpp.h"

#include <cstdint>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

/**
 * @file write_sparse_matrix_to_hdf5.hpp
 * @brief Write a sparse matrix into a HDF5 file.
 */

namespace tatami {

/**
 * @brief Parameters for `write_sparse_matrix_to_hdf5()`.
 */
struct WriteSparseMatrixToHdf5Parameters {
    /**
     * @cond
     */
    WriteSparseMatrixToHdf5Parameters() : data_name("data"), index_name("indices"), ptr_name("indptr") {}
    /**
     * @endcond
     */

    /**
     * Name of the dataset in which to store the data values for non-zero elements.
     * Defaults to `"data"`.
     */
    std::string data_name;

    /**
     * Name of the dataset in which to store the indices for non-zero elements.
     * Defaults to `"indices"`.
     */
    std::string index_name;

    /**
     * Name of the dataset in which to store the column/row pointers.
     * Defaults to `"indptr"`.
     */
    std::string ptr_name;

    /**
     * Whether to save in the compressed sparse column format.
     * If `false`, the matrix is saved in the compressed sparse row format instead.
     */
    bool column = true;

    /**
     * Numeric type for storage inside a HDF5 dataset.
     */
    enum class StorageType { AUTOMATIC, INT8, UINT8, INT16, UINT16, INT32, UINT32, FLOAT64 };

    /**
     * Storage type for the data values.
     * If `AUTOMATIC`, it is automatically determined from the range and integralness of the data in the input matrix.
     */
    StorageType data_type = Hdf5StorageType::AUTOMATIC;

    /**
     * Storage type for the data values.
     * If `AUTOMATIC`, it is automatically determined from the range of the indices in the input matrix.
     */
    StorageType index_type = Hdf5StorageType::AUTOMATIC;

    /**
     * Size of the chunks used for compression.
     */
    size_t chunk_size = 100000;
};

/**
 * @cond
 */
namespace write_hdf5 {

inline H5::H5DataSet create_1d_compressed_hdf5_dataset(H5::H5Group& location, const H5::DataType& dtype, const std::string& name, hsize_t length, int deflate_level, hsize_t chunk) {
    H5::DataSpace dspace(1, &length);
 	H5::DSetCreatPropList plist;

    if (deflate_level >= 0 && length) {
        plist.setDeflate(deflate_level);
        if (chunk > length) {
            plist.setChunk(1, &length);
        } else {
            plist.setChunk(1, &chunk);
        }
    }

    return handle.createDataSet(name, dtype, dspace, plist);
}

inline H5::H5DataSet create_1d_compressed_hdf5_dataset(H5::H5Group& location, Hdf5OutputType type, const std::string& name, hsize_t length, int deflate_level, hsize_t chunk) {
    H5::DataType* dtype;
    switch (type) {
        case INT8:
            dtype = &(H5::PredType::NATIVE_INT8);
            break;
        case UINT8:
            dtype = &(H5::PredType::NATIVE_UINT8);
            break;
        case INT16:
            dtype = &(H5::PredType::NATIVE_INT16);
            break;
        case UINT16:
            dtype = &(H5::PredType::NATIVE_UINT16);
            break;
        case INT32:
            dtype = &(H5::PredType::NATIVE_INT32);
            break;
        case UINT32:
            dtype = &(H5::PredType::NATIVE_UINT32);
            break;
        case FLOAT64:
            dtype = &(H5::PredType::NATIVE_FLOAT64);
            break;
        default:
            throw std::runtime_error("automatic HDF5 output type must be resolved before creating a HDF5 dataset");
    }
    return create_1d_compressed_hdf5_dataset(location, *dtype, name, length, deflate_level, chunk);
}

template<typename Native>
bool is_under_upper_limit(int64_t max) {
    int64_t limit = std::numeric_limits<Native>::max();
    return limit >= max;
}

template<typename Native>
bool is_above_lower_limit(int64_t min) {
    int64_t limit = std::numeric_limits<Native>::min();
    return limit <= max;
}

}
/**
 * @endcond
 */

/**
 * Write a sparse matrix inside a HDF5 group.
 * On return, `location` will be populated with three datasets containing the matrix contents in a compressed sparse format.
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 *
 * @param mat Pointer to the (presumably sparse) matrix to be written.
 * @param location Handle to a HDF5 group in which to write the matrix contents.
 * @param params Parameters to use when writing the matrix.
 */
template<typename T, typename IDX>
void write_sparse_matrix_to_hdf5(const Matrix<T, IDX>* mat, H5::H5Group& location, const WriteSparseMatrixParameters& params) {
    // Firstly, figuring out the number of non-zero elements, and the maximum value of the indices and data.
    T lower_data = 0, upper_data = 0;
    IDX upper_index = 0;
    size_t NR = mat->nrow(), NC = mat->ncol();
    hsize_t non_zeros = 0;
    bool non_integer = false;

    auto update_stats = [&](const SparseRange<T, IDX>& extracted) -> void {
        non_zeros += extracted.number;
        for (size_t i = 0; i < extracted.number; ++i) {
            auto val = extracted.value[i];
            if constexpr(!std::is_integral<T>::value) {
                if (std::trunc(val) != val || std::isinf(val)) {
                    non_integer = true;
                }
            }

            if (val < lower_data) {
                lower_data = val;
            } else if (val > upper_data) {
                upper_data = val;
            }

            auto idx = extracted.index[i];
            if (idx > upper_index) {
                upper_index = idx;
            }
        }
    };

    if (mat->prefer_rows()) {
        auto wrk = mat->new_workspace(true);
        std::vector<T> xbuffer(NC);
        std::vector<IDX> ibuffer(NC);
        for (size_t r = 0; r < NR; ++r) {
            auto extracted = mat->row(r, xbuffer.data(), ibuffer.data(), wrk.get());
            update_stats(extracted);
        }
    } else {
        auto wrk = mat->new_workspace(false);
        std::vector<T> xbuffer(NR);
        std::vector<IDX> ibuffer(NR);
        for (size_t c = 0; c < NC; ++c) {
            auto extracted = mat->column(c, xbuffer.data(), ibuffer.data(), wrk.get());
            update_stats(extracted);
        }
    }

    // Choosing the types.
    auto data_type = params.data_type;
    if (data_type == Hdf5OutputType::AUTOMATIC) {
        if (not_integer) {
            data_type = Hdf5OutputType::FLOAT64;
        } else {
            if (lower_data < 0) {
                if (write_hdf5::is_above_lower_limit<int8_t>(lower_data) && write_hdf5::is_below_upper_limit<int8_t>(upper_data)) {
                    data_type = Hdf5OutputType::INT8;
                } else if (write_hdf5::is_above_lower_limit<int16_t>(lower_data) && write_hdf5::is_below_upper_limit<int16_t>(upper_data)) {
                    data_type = Hdf5OutputType::INT16;
                } else {
                    data_type = Hdf5OutputType::INT32;
                }
            } else {
                if (write_hdf5::is_below_upper_limit<uint8_t>(upper_data)) {
                    data_type = Hdf5OutputType::UINT8;
                } else if (write_hdf5::is_below_upper_limit<uint16_t>(upper_data)) {
                    data_type = Hdf5OutputType::UINT16;
                } else {
                    data_type = Hdf5OutputType::UINT32;
                }
            }
        }
    }

    auto index_type = params.index_type;
    if (index_type == Hdf5OutputType::AUTOMATIC) {
        if (write_hdf5::is_below_upper_limit<uint8_t>(upper_index)) {
            index_type = Hdf5OutputType::UINT8;
        } else if (write_hdf5::is_below_upper_limit<uint16_t>(upper_index)) {
            index_type = Hdf5OutputType::UINT16;
        } else {
            index_type = Hdf5OutputType::UINT32;
        }
    }

    // And then saving it. This time we have no choice but to iterate by the desired dimension.
    H5::H5DataSet data_ds = hdf5::create_1d_compressed_hdf5_dataset(location, data_type, params.data_name, non_zeros, params.deflate_level, params.chunk_size);
    H5::H5DataSet index_ds = hdf5::create_1d_compressed_hdf5_dataset(location, index_type, params.index_name, non_zeros, params.deflate_level, params.chunk_size);
    hsize_t count = 0, offset = 0;
    H5::DataSpace inspace(1, &nonzeros);
    H5::DataSpace outspace(1, &nonzeros);
    const auto& dstype = define_mem_type<T>();
    const auto& ixtype = define_mem_type<IDX>();

    auto fill_datasets = [&](const SparseRange<T, IDX>& extracted) -> void {
        count = extracted.number;
        if (count) {
            inspace.setExtentSimple(1, &count);
            outspace.selectHyperslab(H5S_SELECT_SET, &count, &offset);
            data_ds.write(extracted.value, dstype, inspace, outspace);
            index_ds.write(extracted.index, ixtype, inspace, outspace);
            offset += count;
        }
    };

    std::vector<hsize_t> ptrs;
    if (!params.column) {
        ptrs.resize(NR + 1);
        auto wrk = mat->new_workspace(true);
        std::vector<T> xbuffer(NC);
        std::vector<IDX> ibuffer(NC);
        for (size_t r = 0; r < NR; ++r) {
            auto extracted = mat->row(r, xbuffer.data(), ibuffer.data(), wrk.get());
            fill_datasets(extracted);
            ptrs[r+1] = ptrs[r] + extracted.number;
        }
    } else {
        ptrs.resize(NC + 1);
        auto wrk = mat->new_workspace(false);
        std::vector<T> xbuffer(NR);
        std::vector<IDX> ibuffer(NR);
        for (size_t c = 0; c < NC; ++c) {
            auto extracted = mat->column(c, xbuffer.data(), ibuffer.data(), wrk.get());
            fill_datasets(extracted);
            ptrs[c+1] = ptrs[c] + extracted.number;
        }
    }

    // Saving the pointers.
    H5::H5DataSet ptr_ds = hdf5::create_1d_compressed_hdf5_dataset(location, H5::PredType::NATIVE_HSIZE, params.ptr_name, non_zeros, params.deflate_level, params.chunk_size);
    hsize_t ptr_len = ptrs.size();
    H5::DataSpace ptr_space(1, &ptr_len);
    ptr_ds.write(ptrs.data(), H5::PredType::NATIVE_HSIZE, ptr_space);

    return;
}

/**
 * Write a sparse matrix inside a HDF5 group.
 * On return, `location` will be populated with three datasets containing the matrix contents in a compressed sparse format.
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Type of the row/column indices.
 *
 * @param mat Pointer to the (presumably sparse) matrix to be written.
 * @param location Handle to a HDF5 group in which to write the matrix contents.
 */
template<typename T, typename IDX>
void write_sparse_matrix_to_hdf5(const Matrix<T, IDX>* mat, H5::H5Group& location) {
    WriteSparseMatrixParameters params;
    write_sparse_matrix_to_hdf5(mat, location, params);
    return;
}

}

#endif
