#ifndef TATAMI_UNARY_ARITH_HELPERS_H
#define TATAMI_UNARY_ARITH_HELPERS_H

#include "../arith_utils.hpp"
#include <vector>
#include <limits>

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

template<DelayedArithOp op_, bool right_, typename Value_, typename Scalar_>
constexpr bool delayed_arith_unsupported_division_by_zero() {
    return !std::numeric_limits<Value_>::is_iec559 && op_ == DelayedArithOp::DIVIDE && !right_;
}

template<DelayedArithOp op_, bool right_, typename Value_, typename Scalar_>
bool delayed_arith_actual_sparse(Scalar_ scalar) {
    if constexpr(delayed_arith_unsupported_division_by_zero<op_, right_, Value_, Scalar_>()) {
        // If we didn't catch this case, the else() condition would be dividing
        // by zero in a Value_ that doesn't support it, and that would be visible
        // at compile time - possibly resulting in compiler warnings. So we
        // declare that this is always non-sparse, and hope that the equivalent
        // zero() method doesn't get called.
        return false;
    } else {
        // Empirically testing this, to accommodate special values (e.g., NaN, Inf) for scalars.
        Value_ output = 0;
        delayed_arith_run<op_, right_>(output, scalar);
        return output == 0;
    }
}

template<DelayedArithOp op_, bool right_, typename Value_, typename Scalar_>
Value_ delayed_arith_zero(Scalar_ scalar) {
    if constexpr(delayed_arith_unsupported_division_by_zero<op_, right_, Value_, Scalar_>()) {
        // Avoid potential problems with division by zero that can be detected
        // at compile time (e.g., resulting in unnecessary compiler warnings).
        throw std::runtime_error("division by zero is not supported with IEEE-754 floats");
    } else {
        Value_ output = 0;
        delayed_arith_run<op_, right_>(output, scalar);
        return output;
    }
}

template<DelayedArithOp op_, bool right_, typename Value_, typename Scalar_>
constexpr bool delayed_arith_always_dense() {
    // If the scalar is operated on by the matrix, return true if zeros in the matrix yield non-zero results.
    if constexpr(!right_) {
        if constexpr(op_ == DelayedArithOp::DIVIDE) {
            return true;
        } else if constexpr(op_ == DelayedArithOp::POWER) {
            return true;
        } else if constexpr(op_ == DelayedArithOp::MODULO) {
            return true;
        } else if constexpr(op_ == DelayedArithOp::INTEGER_DIVIDE) {
            return true;
        }
    }

    return false;    
}

template<DelayedArithOp op_, bool right_, typename Value_, typename Scalar_>
constexpr bool delayed_arith_always_sparse() {
    // Multiplication is always sparse if the Scalar_ type cannot have special values.
    if constexpr(op_ == DelayedArithOp::MULTIPLY && 
        !std::numeric_limits<Scalar_>::has_infinity &&
        !std::numeric_limits<Scalar_>::has_quiet_NaN && 
        !std::numeric_limits<Scalar_>::has_signaling_NaN)
    {
        return true;
    }

    return false;
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
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar value.
 */
template<DelayedArithOp op_, bool right_, typename Value_ = double, typename Scalar_ = Value_>
struct DelayedArithScalarHelper {
    /**
     * @param s Scalar value to be added.
     */
    DelayedArithScalarHelper(Scalar_ s) : scalar(s) {
        still_sparse = delayed_arith_actual_sparse<op_, right_, Value_>(scalar);
    }

private:
    const Scalar_ scalar;
    bool still_sparse;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool always_dense = delayed_arith_always_dense<op_, right_, Value_, Scalar_>();

    static constexpr bool always_sparse = delayed_arith_always_sparse<op_, right_, Value_, Scalar_>();

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
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_arith_run_simple<op_, right_>(scalar, length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_arith_run_simple<op_, right_>(scalar, number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return delayed_arith_zero<op_, right_, Value_>(scalar);
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
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 */
template<DelayedArithOp op_, bool right_, int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
struct DelayedArithVectorHelper {
    /**
     * @param v Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedArithVectorHelper(Vector_ v) : vec(std::move(v)) {
        for (auto x : vec) {
            if (!delayed_arith_actual_sparse<op_, right_, Value_>(x)) {
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

    typedef typename std::remove_reference<decltype(std::declval<Vector_>()[0])>::type Scalar_;

    static constexpr bool always_dense = delayed_arith_always_dense<op_, right_, Value_, Scalar_>();

    static constexpr bool always_sparse = delayed_arith_always_sparse<op_, right_, Value_, Scalar_>();

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
    template<bool accrow_, typename Index_, typename ExtractType_>
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

    template<bool accrow_, typename Index_>
    void sparse(Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_arith_run_simple<op_, right_>(vec[idx], number, buffer);

        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_arith_run<op_, right_>(buffer[i], vec[indices[i]]);
            }
        }
    }

    template<bool, typename Index_>
    Value_ zero(Index_ idx) const {
        return delayed_arith_zero<op_, right_, Value_>(vec[idx]);
    }
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be added.
 * @return A helper class for delayed scalar addition.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::ADD, true, Value_, Scalar_> make_DelayedAddScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::ADD, true, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be subtracted.
 * @return A helper class for delayed scalar subtraction.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::SUBTRACT, right_, Value_, Scalar_> make_DelayedSubtractScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::SUBTRACT, right_, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be multiplied.
 * @return A helper class for delayed scalar multiplication.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::MULTIPLY, true, Value_, Scalar_> make_DelayedMultiplyScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::MULTIPLY, true, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be divided.
 * @return A helper class for delayed scalar division.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::DIVIDE, right_, Value_, Scalar_> make_DelayedDivideScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::DIVIDE, right_, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the power transformation.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be power transformed.
 * @return A helper class for delayed scalar power transformation.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::POWER, right_, Value_, Scalar_> make_DelayedPowerScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::POWER, right_, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the modulus.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be modulo transformed.
 * @return A helper class for delayed scalar modulus.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::MODULO, right_, Value_, Scalar_> make_DelayedModuloScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::MODULO, right_, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the integer division.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be integer divided.
 * @return A helper class for delayed scalar integer division.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedArithScalarHelper<DelayedArithOp::INTEGER_DIVIDE, right_, Value_, Scalar_> make_DelayedIntegerDivideScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::INTEGER_DIVIDE, right_, Value_, Scalar_>(std::move(s));
}

/**
 * @tparam margin_ Matrix dimension along which the addition is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to be added to the rows/columns.
 * @return A helper class for delayed vector addition.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::ADD, true, margin_, Value_, Vector_> make_DelayedAddVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::ADD, true, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam margin_ Matrix dimension along which the subtraction is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to subtract from (or be subtracted by) the rows/columns.
 * @return A helper class for delayed vector subtraction.
 */
template<bool right_, int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::SUBTRACT, right_, margin_, Value_, Vector_> make_DelayedSubtractVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::SUBTRACT, right_, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam margin_ Matrix dimension along which the multiplication is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to multiply the rows/columns.
 * @return A helper class for delayed vector multiplication.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::MULTIPLY, true, margin_, Value_, Vector_> make_DelayedMultiplyVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::MULTIPLY, true, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam margin_ Matrix dimension along which the division is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to divide (or be divided by) the rows/columns.
 * @return A helper class for delayed vector division.
 */
template<bool right_, int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::DIVIDE, right_, margin_, Value_, Vector_> make_DelayedDivideVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::DIVIDE, right_, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the power transformation.
 * @tparam margin_ Matrix dimension along which the power transformation is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to use in the power transformation of the rows/columns.
 * @return A helper class for delayed vector power transformation.
 */
template<bool right_, int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::POWER, right_, margin_, Value_, Vector_> make_DelayedPowerVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::POWER, right_, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the modulus.
 * @tparam margin_ Matrix dimension along which the modulus is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to use in the modulus of the rows/columns.
 * @return A helper class for delayed vector modulus.
 */
template<bool right_, int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::MODULO, right_, margin_, Value_, Vector_> make_DelayedModuloVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::MODULO, right_, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the integer division.
 * @tparam margin_ Matrix dimension along which the integer division is to occur, see `DelayedArithVectorHelper`.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to integer divide (or be integer divided by) the rows/columns.
 * @return A helper class for delayed vector division.
 */
template<bool right_, int margin_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedArithVectorHelper<DelayedArithOp::INTEGER_DIVIDE, right_, margin_, Value_, Vector_> make_DelayedIntegerDivideVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::INTEGER_DIVIDE, right_, margin_, Value_, Vector_>(std::move(v));
}

}

#endif
