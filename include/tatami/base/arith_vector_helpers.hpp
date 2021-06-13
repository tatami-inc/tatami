#ifndef TATMI_ARITH_VECTOR_HELPER_H
#define TATMI_ARITH_VECTOR_HELPER_H

/**
 * @file arith_vector_helpers.hpp
 *
 * Helper functions focusing on arithmetic operations with a vector parallel to the rows or columns,
 * to be used as the `OP` in the `DelayedIsometricOp` class.
 */

#include <vector>

namespace tatami {

/**
 * @brief Add a vector along the rows or columns of a matrix.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after addition.
 * @tparam MARGIN Dimension along which the addition is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and the same value is added to all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam V Class of the vector holding the values to be added.
 */
template<int MARGIN, typename T = double, class V = std::vector<T> >
struct DelayedAddVectorHelper {
    /**
     * @param v Vector of values to be added.
     * This should be of length equal to the number of rows in the matrix if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedAddVectorHelper(V v) : vec(std::move(v)) {}

    /**
     * @param r Row index.
     * @param c Column index.
     * @param val Matrix value to be added.
     *
     * @return `val` plus the vector element at `r` (if `MARGIN = 0`) or at `c` (if `MARGIN = 1`).
     */
    T operator()(size_t r, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            return val + vec[r]; 
        } else {
            return val + vec[c]; 
        }
    }

    /**
     * Addition is always assumed to discard structural sparsity, even when the added value is zero.
     */
    static const bool sparse = false; 
private:
    const V vec;
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 * See the `tatami::DelayedAddVectorHelper` documentation for more details on the arguments.
 */
template<int MARGIN, typename T = double, class V = std::vector<T> >
DelayedAddVectorHelper<MARGIN, T, V> make_DelayedAddVectorHelper(V v) {
    return DelayedAddVectorHelper<MARGIN, T, V>(std::move(v));
}

/**
 * @brief Subtract a vector from the rows/columns of a matrix, or vice versa.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after subtraction.
 * @tparam RIGHT Should the vector be subtracted from the matrix rows/columns?
 * If `false`, the matrix rows/columns are instead subtracted from the vector.
 * @tparam MARGIN Dimension along which the subtraction is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and the same value is subtracted from all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam V Class of the vector holding the values to use in the subtraction.
 */
template<bool RIGHT, int MARGIN, typename T = double,  class V = std::vector<T> >
struct DelayedSubtractVectorHelper {
    /**
     * @param v Vector of values to use for subtraction.
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedSubtractVectorHelper(V v) : vec(std::move(v)) {}

    /**
     * @param r Row index.
     * @param c Column index.
     * @param val Matrix value to use in the subtraction.
     *
     * @return If `RIGHT = true`, `val` minus the vector element at `r` (if `MARGIN = 0`) or at `c` (if `MARGIN = 1`).
     * If `RIGHT = false`, the vector element minus `val` is returned instead.
     */
    T operator()(size_t r, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            if constexpr(RIGHT) {
                return val - vec[r];
            } else {
                return vec[r] - val;
            }
        } else {
            if constexpr(RIGHT) {
                return val - vec[c];
            } else {
                return vec[c] - val;
            }
        }
    }

    /**
     * Subtraction is always assumed to discard structural sparsity, even when the subtracted value is zero.
     */
    static const bool sparse = false; 
private:
    const V vec;
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 * See the `tatami::DelayedSubtractVectorHelper` documentation for more details on the arguments.
 */
template<bool RIGHT, int MARGIN, typename T = double, class V = std::vector<T> >
DelayedSubtractVectorHelper<RIGHT, MARGIN, T, V> make_DelayedSubtractVectorHelper(V v) {
    return DelayedSubtractVectorHelper<RIGHT, MARGIN, T, V>(std::move(v));
}

/**
 * @brief Multiply a vector along the rows or columns of a matrix.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after multiplication.
 * @tparam MARGIN Dimension along which the multipication is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and the same value is multiplied to all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam V Class of the vector holding the values to use for multiplication
 */
template<int MARGIN, typename T = double, class V = std::vector<T> >
struct DelayedMultiplyVectorHelper {
    /**
     * @param v Vector of values to use for multiplication.
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedMultiplyVectorHelper(V v) : vec(std::move(v)) {}

    /**
     * @param r Row index.
     * @param c Column index.
     * @param val Matrix value to use in the multiplication.
     *
     * @return `val` multiplied by the vector element at `r` (if `MARGIN = 0`) or at `c` (if `MARGIN = 1`).
     */
    T operator()(size_t r, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            return val * vec[r]; 
        } else {
            return val * vec[c]; 
        }
    }

    /**
     * Multiplication is always assumed to preserve structural sparsity.
     */
    static const bool sparse = true;
private:
    const V vec;
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 * See the `tatami::DelayedMultiplyVectorHelper` documentation for more details on the arguments.
 */
template<int MARGIN, typename T = double, class V = std::vector<T> >
DelayedMultiplyVectorHelper<MARGIN, T, V> make_DelayedMultiplyVectorHelper(V v) {
    return DelayedMultiplyVectorHelper<MARGIN, T, V>(std::move(v));
}

/**
 * @brief Divide the rows/columns of a matrix by a vector, or vice versa.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after division.
 * @tparam RIGHT Should the matrix rows/columns be divided by the vector?
 * If `false`, the vector is instead divided by the matrix rows/columns.
 * @tparam MARGIN Dimension along which the division is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and all entries in the same row of the matrix are divided by the same vector entry.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam V Class of the vector holding the values to use in the division.
 */
template<bool RIGHT, int MARGIN, typename T = double, class V = std::vector<T> >
struct DelayedDivideVectorHelper {
    /**
     * @param v Vector of values to use for division.
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedDivideVectorHelper(V v) : vec(std::move(v)) {}

    /**
     * @param r Row index.
     * @param c Column index.
     * @param val Matrix value to use in the division.
     *
     * @return If `RIGHT = true`, `val` divided by the vector element at `r` (if `MARGIN = 0`) or at `c` (if `MARGIN = 1`).
     * If `RIGHT = false`, the vector element divided by `val` is returned instead.
     */
    T operator()(size_t r, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            if constexpr(RIGHT) {
                return val / vec[r];
            } else {
                return vec[r] / val;
            }
        } else {
            if constexpr(RIGHT) {
                return val / vec[c];
            } else {
                return vec[c] / val;
            }
        }
    }

    /**
     * Division is always assumed to preserve structural sparsity.
     * Non-finite or zero `scalar` values are not considered here.
     */
    static const bool sparse = true;
private:
    const V vec;
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 * See the `tatami::DelayedDivideVectorHelper` documentation for more details on the arguments.
 */
template<bool RIGHT, int MARGIN, typename T = double, class V = std::vector<T> >
DelayedDivideVectorHelper<RIGHT, MARGIN, T, V> make_DelayedDivideVectorHelper(V v) {
    return DelayedDivideVectorHelper<RIGHT, MARGIN, T, V>(std::move(v));
}


}

#endif
