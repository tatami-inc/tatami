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
// We deliberately use an auto type so as to defer a decision on what the output
// type should be; an appropriate coercion is left to the caller classes. 
template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
auto delayed_arithmetic(Value_ val, Scalar_ scalar) { 
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

// Some of the helpers need to divide by zero to compute the fill value. The
// only way to guarantee that division by zero is supported is with IEEE
// floats, otherwise it would be undefined behavior at compile time, and the
// compiler could do anything, including refusing to compile it. So, we hide
// any '0/0' behind a constexpr to ensure that the compiler doesn't see it.
template<ArithmeticOperation op_, typename Value_, typename Scalar_>
constexpr bool has_unsafe_divide_by_zero() {
    if constexpr(std::numeric_limits<Value_>::is_iec559) {
        return false;
    }
    if constexpr(std::numeric_limits<Scalar_>::is_iec559) {
        return false;
    }

    // MODULO (and POWER, for negative powers) also involve division by zero,
    // but they return an implementation-defined value; so they're "safe".
    // - https://en.cppreference.com/w/cpp/numeric/math/pow
    // - https://en.cppreference.com/w/cpp/numeric/math/fmod

    if constexpr(op_ == ArithmeticOperation::DIVIDE) {
        return true;
    }
    if constexpr(op_ == ArithmeticOperation::INTEGER_DIVIDE) {
        return true;
    }

    return false;
}

// COMMENT: if the fill value is visible at a compile time, the coercion to the
// output type is also visible. If the output type cannot hold the fill value,
// we could potentially get more compile-time UB. (This mostly concerns NaNs or
// Infs from divide-by-zero, but could apply to overflows or out-of-range casts
// to integer types.) To get around this, we assume that the compiler supports
// the IEC 559 spec, where - according to Annex F.4 of the C standard - the
// cast from float to integer is merely unspecified, not undefined.

/**
 * @endcond
 */
}

#endif
