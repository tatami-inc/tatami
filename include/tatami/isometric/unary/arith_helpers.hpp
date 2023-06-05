#ifndef TATAMI_UNARY_ARITH_HELPERS_H
#define TATAMI_UNARY_ARITH_HELPERS_H

#include "../arith_utils.hpp"

/**
 * @file arith_helpers.hpp
 *
 * @brief Helper classes for delayed unary arithmetic operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 */

namespace tatami {

/**
 * @cond
 */
template<DelayedArithOp op_, bool right_, typename Scalar_, typename Value_, typename Index_>
void delayed_arith_run_simple(Scalar_ scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_arith_run<op_, right_>(buffer[i], scalar);
    }
}

template<DelayedArithOp op_, bool right_, typename Scalar_>
bool delayed_arith_actual_sparse(Scalar_ scalar) {
    if constexpr(op_ == DelayedArithOp::ADD || op_ == DelayedArithOp::SUBTRACT) {
        return scalar == 0;
    } else { // DIVIDE only, as MULTIPLY is always_sparse.
        return scalar != 0; 
    }
}
/**
 * @endcond
 */

/**
 * @brief Delayed scalar arithmetic.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the scalar should be on the right hand side of the arithmetic operation.
 * Ignored for commutative operations, e.g., `ADD` and `MULTIPLY`.
 * @tparam Scalar_ Type of the scalar value.
 */
template<DelayedArithOp op_, bool right_, typename Scalar_>
struct DelayedArithScalarHelper {
    /**
     * @param s Scalar value to be added.
     */
    DelayedArithScalarHelper(Scalar_ s) : scalar(s) {}

private:
    const Scalar_ scalar;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool always_dense = (op_ == DelayedArithOp::DIVIDE && !right_);

    static constexpr bool always_sparse = (op_ == DelayedArithOp::MULTIPLY);

    bool actual_sparse() const {
        return delayed_arith_actual_sparse<op_, right_>(scalar);
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_arith_run_simple<op_, right_>(scalar, length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_arith_run_simple<op_, right_>(scalar, number, buffer);
    }

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_arith_run_simple<op_, right_>(scalar, length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector arithmetic.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the vector's values should be on the right hand side of the arithmetic operation.
 * Ignored for some `op_`.
 * @tparam margin_ Matrix dimension along which the operation is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and that value is subtracted from all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam Vector_ Type of the vector.
 */
template<DelayedArithOp op_, bool right_, int margin_, typename Vector_>
struct DelayedArithVectorHelper {
    /**
     * @param v Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedArithVectorHelper(Vector_ v) : vec(std::move(v)) {
        for (auto x : vec) {
            if (!delayed_arith_actual_sparse<op_, right_>(x)) {
                still_sparse = false;
                break;
            }
        }
    }

private:
    const Vector_ vec;
    bool still_sparse = true;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = (margin_ == 0);

    static constexpr bool needs_column = (margin_ == 1);

    static constexpr bool always_dense = (op_ == DelayedArithOp::DIVIDE && !right_);

    static constexpr bool always_sparse = (op_ == DelayedArithOp::MULTIPLY);

    bool actual_sparse() const {
        return still_sparse;
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool accrow_, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_ idx, ExtractType_ start, Index_ length, Value_* buffer) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_arith_run_simple<op_, right_>(vec[idx], length, buffer);

        } else if constexpr(std::is_same<ExtractType_, Index_>::value) {
            for (Index_ i = 0; i < length; ++i) {
                delayed_arith_run<op_, right_>(buffer[i], vec[i + start]);
            }

        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_arith_run<op_, right_>(buffer[i], vec[start[i]]);
            }
        }
    }

    template<bool accrow_, typename Value_, typename Index_>
    void sparse(Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_arith_run_simple<op_, right_>(vec[idx], number, buffer);

        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_arith_run<op_, right_>(buffer[i], vec[indices[i]]);
            }
        }
    }

    template<bool accrow_, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_ idx, ExtractType_&& start, Index_ length, Value_* buffer) const {
        dense<accrow_>(idx, std::forward<ExtractType_>(start), length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be added.
 * @return A helper class for delayed scalar addition.
 */
template<typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::ADD, true, Scalar_> make_DelayedAddScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::ADD, true, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be subtracted.
 * @return A helper class for delayed scalar subtraction.
 */
template<bool right_, typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::SUBTRACT, right_, Scalar_> make_DelayedSubtractScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::SUBTRACT, right_, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be multiplied.
 * @return A helper class for delayed scalar multiplication.
 */
template<typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::MULTIPLY, true, Scalar_> make_DelayedMultiplyScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::MULTIPLY, true, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be divided.
 * @return A helper class for delayed scalar division.
 */
template<bool right_, typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::DIVIDE, right_, Scalar_> make_DelayedDivideScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::DIVIDE, right_, Scalar_>(std::move(s));
}

/**
 * @tparam margin_ Matrix dimension along which the addition is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to be added to the rows/columns.
 * @return A helper class for delayed vector addition.
 */
template<int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::ADD, true, margin_, Vector_> make_DelayedAddVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::ADD, true, margin_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam margin_ Matrix dimension along which the subtraction is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to subtract from (or be subtracted by) the rows/columns.
 * @return A helper class for delayed vector subtraction.
 */
template<bool right_, int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::SUBTRACT, right_, margin_, Vector_> make_DelayedSubtractVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::SUBTRACT, right_, margin_, Vector_>(std::move(v));
}

/**
 * @tparam margin_ Matrix dimension along which the multiplication is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to multiply the rows/columns.
 * @return A helper class for delayed vector multiplication.
 */
template<int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::MULTIPLY, true, margin_, Vector_> make_DelayedMultiplyVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::MULTIPLY, true, margin_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam margin_ Matrix dimension along which the division is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to divide (or be divided by) the rows/columns.
 * @return A helper class for delayed vector division.
 */
template<bool right_, int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::DIVIDE, right_, margin_, Vector_> make_DelayedDivideVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::DIVIDE, right_, margin_, Vector_>(std::move(v));
}

}

#endif
