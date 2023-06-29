#ifndef TATAMI_ARITH_UTILS_HPP
#define TATAMI_ARITH_UTILS_HPP

#include <cmath>

/**
 * @file arith_utils.hpp
 *
 * @brief Utilities for delayed arithmetic operations.
 */

namespace tatami {

/**
 * Type of the delayed arithmetic operation.
 */
enum class DelayedArithOp : char { 
    ADD, 
    SUBTRACT,
    MULTIPLY,
    DIVIDE,
    POWER,
    MODULO,
    INTEGER_DIVIDE
};

/**
 * @cond
 */
template<DelayedArithOp op_, bool right_, typename Scalar_, typename Value_>
void delayed_arith_run(Value_& val, Scalar_ scalar) {
    if constexpr(op_ == DelayedArithOp::ADD) {
        val += scalar;
    } else if constexpr(op_ == DelayedArithOp::MULTIPLY) {
        val *= scalar;
    } else if constexpr(op_ == DelayedArithOp::SUBTRACT) {
        if constexpr(right_) {
            val -= scalar;
        } else {
            val = scalar - val;
        }
    } else if constexpr(op_ == DelayedArithOp::DIVIDE) {
        // Assume that either Value_ is an IEEE-754 float, or that division by
        // zero is impossible in this context. We don't apply manual checks
        // here to avoid performance degradation; we also don't check that the
        // other operations yield a value that doesn't overflow/underflow, so
        // it would be odd to make an exception for div-by-zero errors.
        if constexpr(right_) {
            val /= scalar;
        } else {
            val = scalar / val;
        }
    } else if constexpr(op_ == DelayedArithOp::POWER) {
        if constexpr(right_) {
            val = std::pow(val, scalar);
        } else {
            val = std::pow(scalar, val);
        }
    } else if constexpr(op_ == DelayedArithOp::MODULO) {
        if constexpr(right_) {
            val = std::fmod(val, scalar);
        } else {
            val = std::fmod(scalar, val);
        }
    } else if constexpr(op_ == DelayedArithOp::INTEGER_DIVIDE) {
        if constexpr(right_) {
            val = std::floor(val / scalar);
        } else {
            val = std::floor(scalar / val);
        }
    }
}
/**
 * @endcond
 */
}

#endif
