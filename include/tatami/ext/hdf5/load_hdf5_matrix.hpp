#ifndef TATAMI_LOAD_HDF5_MATRIX_HPP
#define TATAMI_LOAD_HDF5_MATRIX_HPP

#include "H5Cpp.h"

#include <string>
#include <cstdint>
#include <type_traits>
#include <cmath>

#include "../../base/Matrix.hpp"
#include "utils.hpp"
#include "../../base/sparse/CompressedSparseMatrix.hpp"
#include "../../base/dense/DenseMatrix.hpp"

/**
 * @file load_hdf5_matrix.hpp
 *
 * @brief Load HDF5 data into `Matrix` objects.
 */

namespace tatami {

/**
 * Load a `CompressedSparseMatrix` from a HDF5 file.
 *
 * @tparam ROW Whether the matrix is stored in compressed sparse row format.
 * @tparam T Type of the matrix values in the `Matrix` interface.
 * @tparam IDX Type of the row/column indices.
 * @tparam U Vector type for storing the values of the non-zero elements.
 * Elements of this vector may be of a different type than `T` for more efficient storage.
 * @tparam V Vector type for storing the indices.
 * Elements of this vector may be of a different type than `IDX` for more efficient storage.
 * @tparam W Vector type for storing the index pointers.
 *
 * @param nr Number of rows in the matrix.
 * @param nc Number of columns in the matrix.
 * @param file Path to the file.
 * @param vals Name of the 1D dataset inside `file` containing the non-zero elements.
 * @param idx Name of the 1D dataset inside `file` containing the indices of the non-zero elements.
 * If `ROW = true`, this should contain column indices sorted within each row, otherwise it should contain row indices sorted within each column.
 * @param ptr Name of the 1D dataset inside `file` containing the index pointers for the start and end of each row (if `ROW = true`) or column (otherwise).
 * This should have length equal to the number of rows (if `ROW = true`) or columns (otherwise).
 *
 * @return A `CompressedSparseMatrix` containing all values and indices in memory.
 * This differs from a `HDF5CompressedSparseMatrix`, where the loading of data is deferred until requested.
 */
template<bool ROW, typename T, typename IDX = int, class U = std::vector<T>, class V = std::vector<IDX>, class W = std::vector<size_t> >
CompressedSparseMatrix<ROW, T, IDX, U, V, W> load_hdf5_compressed_sparse_matrix(
    size_t nr, 
    size_t nc, 
    const std::string& file, 
    const std::string& vals, 
    const std::string& idx, 
    const std::string& ptr) 
{
    H5::H5File file_handle(file, H5F_ACC_RDONLY);

    auto dhandle = HDF5::open_and_check_dataset<false>(file_handle, vals);
    const size_t nonzeros = HDF5::get_array_dimensions<1>(dhandle, "vals")[0];
    U x(nonzeros);
    dhandle.read(x.data(), HDF5::define_mem_type<typename U::value_type>());
    
    auto ihandle = HDF5::open_and_check_dataset<true>(file_handle, idx);
    if (HDF5::get_array_dimensions<1>(ihandle, "idx")[0] != nonzeros) {
        throw std::runtime_error("number of non-zero elements is not consistent between 'data' and 'idx'");
    }
    V i(nonzeros);
    ihandle.read(i.data(), HDF5::define_mem_type<typename V::value_type>());

    auto phandle = HDF5::open_and_check_dataset<true>(file_handle, ptr);
    const size_t ptr_size = HDF5::get_array_dimensions<1>(phandle, "ptr")[0];
    if (ptr_size != (ROW ? nr : nc) + 1) {
        throw std::runtime_error("'ptr' dataset should have length equal to the number of " + (ROW ? std::string("rows") : std::string("columns")) + " plus 1");
    }

    // Because HDF5 doesn't have a native type for size_t.
    W p(ptr_size);
    if constexpr(std::is_same<size_t, typename W::value_type>::value) {
        if constexpr(std::is_same<size_t, hsize_t>::value) {
            phandle.read(p.data(), H5::PredType::NATIVE_HSIZE);
        } else {
            std::vector<hsize_t> p0(ptr_size);
            phandle.read(p0.data(), H5::PredType::NATIVE_HSIZE);
            std::copy(p0.begin(), p0.end(), p.begin());
        }
    } else {
        phandle.read(p.data(), HDF5::define_mem_type<typename W::value_type>());
    }

    return CompressedSparseMatrix<ROW, T, IDX, U, V, W>(nr, nc, std::move(x), std::move(i), std::move(p));
}

/**
 * Load a `DenseMatrix` from a HDF5 DataSet.
 *
 * @tparam T Type of the matrix values in the `Matrix` interface.
 * @tparam IDX Type of the row/column indices.
 * @tparam V Vector type for storing the matrix values.
 * This may be different from `T` for more efficient storage.
 * @tparam transpose Whether the dataset is transposed in its storage order, i.e., rows in HDF5 are columns in the matrix.
 *
 * @param file Path to the HDF5 file.
 * @param name Name of the dataset inside the file.
 * This should refer to a 2-dimensional dataset of integer or floating-point type.
 *
 * @return A `DenseMatrix` where all values are in memory.
 * This differs from a `HDF5DenseMatrix`, where the loading of data is deferred until requested.
 */
template<typename T, typename IDX = int, class V = std::vector<T>, bool transpose = false>
DenseMatrix<!transpose, T, IDX, V> load_hdf5_dense_matrix(const std::string& file, const std::string& name) {
    H5::H5File fhandle(file, H5F_ACC_RDONLY);
    auto dhandle = HDF5::open_and_check_dataset<false>(fhandle, name);

    auto dims = HDF5::get_array_dimensions<2>(dhandle, name);
    V values(dims[0] * dims[1]);
    dhandle.read(values.data(), HDF5::define_mem_type<typename V::value_type>());

    if constexpr(transpose) {
        return DenseMatrix<false, T, IDX, V>(dims[1], dims[0], std::move(values));
    } else {
        return DenseMatrix<true, T, IDX, V>(dims[0], dims[1], std::move(values));
    }
}

}

#endif
