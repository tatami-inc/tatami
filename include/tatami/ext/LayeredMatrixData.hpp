#ifndef TATAMI_LAYERED_MATRIX_DATA_HPP
#define TATAMI_LAYERED_MATRIX_DATA_HPP

#include <memory>
#include <vector>
#include "../base/Matrix.hpp"

/**
 * @file LayeredMatrixData.hpp
 *
 * @brief Defines the layered sparse matrix data structure.
 */

namespace tatami {

/**
 * @brief Pointer and permutations for a layered sparse matrix.
 *
 * In a layered sparse matrix, values are stored in the smallest integer type that will fit all entries on the same row.
 * This is achieved by reordering the rows into separate `tatami::CompressedSparseMatrix` submatrices,
 * where each submatrix contains all rows that can fit into a certain integer type.
 * Submatrices are then combined together using an instance of the `tatami::DelayedBind` class.
 * The idea is to reduce memory usage in common situations involving small integers.
 *
 * Note that the rows of the layered matrix will be in a different order from the input data used to construct it.
 * This is a consequence of the reordering to ensure that rows of the same type can be put into the same submatrix.
 * As a result, the `LayeredMatrixData` object will contain a permutation vector indicating the new position of each original row,
 * i.e., `permutation[i]` will specify the new position of the original row `i`.
 * This can be passed to `tatami::DelayedSubset` to restore the original order, if so desired.
 *
 * Currently, only unsigned integers are supported.
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Integer type for the index.
 */
template<typename T = double, typename IDX = int>
struct LayeredMatrixData {
    /**
     * Pointer to the matrix data.
     * Note that rows will not follow their original order, see `permutation`.
     */
    std::shared_ptr<Matrix<T, IDX> > matrix;

    /**
     * Permutation vector indicating the position of each original row in `matrix`.
     * Specifically, for an original row `r`, its new row index in `matrix` is defined as `permutation[i]`.
     */
    std::vector<size_t> permutation;
};

}

#endif
