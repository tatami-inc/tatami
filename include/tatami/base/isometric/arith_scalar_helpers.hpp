#ifndef TATAMI_ARITH_SCALAR_HELPERS_H
#define TATAMI_ARITH_SCALAR_HELPERS_H

#include <limits>

/**
 * @file arith_scalar_helpers.hpp
 *
 * @brief Helper classes for scalar arithmetic operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedIsometricOp` class.
 */

namespace tatami {

/**
 * @brief Add a scalar to all values of a matrix.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after addition.
 */
template<typename T = double>
struct DelayedAddScalarHelper {
    /**
     * @param s Scalar value to be added.
     */
    DelayedAddScalarHelper(T s) : scalar(s) {}

    /**
     * Coordinates are ignored here and are only listed for compatibility purposes.
     *
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value to be added to the scalar.
     *
     * @return `val` plus the scalar.
     */
    T operator()(size_t r, size_t c, T val) const { 
        return val + scalar; 
    }

    /**
     * Addition is always assumed to discard structural sparsity, even when the scalar is zero.
     */
    static const bool sparse_ = false;

    /**
     * This does not require row indices.
     */
    static const bool needs_row_ = false;

    /**
     * This does not require column indices.
     */
    static const bool needs_column_= false;

private:
    const T scalar;
};

/**
 * @brief Multiply a scalar with all values of a matrix.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after multiplication.
 */
template<typename T = double>
struct DelayedMultiplyScalarHelper { 
    /**
     * @param s Scalar value to be multiplied.
     */
    DelayedMultiplyScalarHelper(T s) : scalar(s) {}

    /**
     * Coordinates are ignored here and are only listed for compatibility purposes.
     *
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value to be added to the scalar.
     *
     * @return `val` multiplied by the scalar.
     */
    T operator()(size_t r, size_t c, T val) const { 
        return val * scalar; 
    }

    /**
     * Multiplication is always assumed to preserve structural sparsity.
     * Non-finite `scalar` values are not considered.
     */
    static const bool sparse_ = true;

    /**
     * This does not require row indices.
     */
    static const bool needs_row_ = false;

    /**
     * This does not require column indices.
     */
    static const bool needs_column_= false;

private:
    const T scalar;
};

/**
 * @brief Subtract a scalar from all values of a matrix, or vice versa.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after subtraction.
 * @tparam RIGHT Should the scalar be subtracted from the matrix value?
 * If `false`, the matrix value is subtracted from the scalar.
 */
template<bool RIGHT, typename T = double>
struct DelayedSubtractScalarHelper {
    /**
     * @param s Scalar value to be subtracted from the matrix, or to subtract from the matrix.
     */
    DelayedSubtractScalarHelper(T s) : scalar(s) {}

    /**
     * Coordinates are ignored here and are only listed for compatibility purposes.
     *
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value to use in the subtraction.
     *
     * @return `val` minus the scalar if `RIGHT = true`, otherwise `val` is subtracted from the scalar.
     */
    T operator()(size_t r, size_t c, T val) const { 
        if constexpr(RIGHT) {
            return val - scalar; 
        } else {
            return scalar - val;
        }
    }

    /**
     * Subtraction is always assumed to discard structural sparsity, even when the scalar is zero.
     */
    static const bool sparse_ = false;

    /**
     * This does not require row indices.
     */
    static const bool needs_row_ = false;

    /**
     * This does not require column indices.
     */
    static const bool needs_column_= false;

private:
    const T scalar;
};

/**
 * @brief Divide a scalar from all values of a matrix, or vice versa.
 *
 * This should be used as the `OP` in the `DelayedIsometricOp` class.
 *
 * @tparam T Type to be returned after division.
 * @tparam RIGHT Should the matrix value be divided by the scalar?
 * If `false`, the scalar is divided by the matrix value.
 */
template<bool RIGHT, typename T = double>
struct DelayedDivideScalarHelper { 
    /**
     * @param s Scalar value to use in the division.
     * This should be non-zero.
     */
    DelayedDivideScalarHelper(T s) : scalar(s) {}

    /**
     * Coordinates are ignored here and are only listed for compatibility purposes.
     *
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value to use in the division.
     *
     * @return `val` divided by the scalar if `RIGHT = true`, otherwise the scalar is divided by `val`.
     * If `RIGHT = false` and `val = 0`, an infinite value is returned.
     */
    T operator()(size_t r, size_t c, T val) const { 
        if constexpr(RIGHT) {
            return val / scalar; 
        } else {
            if (val) {
                return scalar / val;
            } else {
                return std::numeric_limits<T>::infinity();
            }
        }
    }

    /**
     * Division on the right is always assumed to preserve structural sparsity.
     * Non-finite or zero `scalar` values are not considered here.
     *
     * Division of the scalar by the matrix value is assumed to discard structural sparsity,
     * as any matrix zeros will yield an infinite value for a non-zero scalar.
     */
    static const bool sparse_ = RIGHT;

    /**
     * This does not require row indices.
     */
    static const bool needs_row_ = false;

    /**
     * This does not require column indices.
     */
    static const bool needs_column_= false;

private:
    const T scalar;
};

}

#endif
