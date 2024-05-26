#ifndef TATAMI_ARITHMETIC_UTILS_HPP
#define TATAMI_ARITHMETIC_UTILS_HPP

#include <cmath>

/**
 * @file arithmetic_utils.hpp
 *
 * @brief Utilities for delayed arithmetic operations.
 */

namespace tatami {

/**
 * Type of arithmetic operation.
 */
enum class ArithmeticOperation : char { 
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
template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
Value_ delayed_arithmetic(Value_ val, Scalar_ scalar) {
    if constexpr(op_ == ArithmeticOperation::ADD) {
        return val + scalar;
    } else if constexpr(op_ == ArithmeticOperation::MULTIPLY) {
        return val * scalar;
    } else if constexpr(op_ == ArithmeticOperation::SUBTRACT) {
        if constexpr(right_) {
            return val - scalar;
        } else {
            return scalar - val;
        }
    } else if constexpr(op_ == ArithmeticOperation::DIVIDE) {
        // Assume that either Value_ is an IEEE-754 float, or that division by
        // zero is impossible in this context. We don't apply manual checks
        // here to avoid performance degradation; we also don't check that the
        // other operations yield a value that doesn't overflow/underflow, so
        // it would be odd to make an exception for div-by-zero errors.
        if constexpr(right_) {
            return val / scalar;
        } else {
            return scalar / val;
        }
    } else if constexpr(op_ == ArithmeticOperation::POWER) {
        if constexpr(right_) {
            return std::pow(val, scalar);
        } else {
            return std::pow(scalar, val);
        }
    } else if constexpr(op_ == ArithmeticOperation::MODULO) {
        if constexpr(right_) {
            return std::fmod(val, scalar);
        } else {
            return std::fmod(scalar, val);
        }
    } else if constexpr(op_ == ArithmeticOperation::INTEGER_DIVIDE) {
        if constexpr(right_) {
            return std::floor(val / scalar);
        } else {
            return std::floor(scalar / val);
        }
    }
}

template<ArithmeticOperation op_, bool right_, typename Scalar_, typename Value_>
void delayed_arithmetic_run(Value_& val, Scalar_ scalar) {
    val = delayed_arithmetic<op_, right_>(val, scalar);
}
/**
 * @endcond
 */
}

#endif
