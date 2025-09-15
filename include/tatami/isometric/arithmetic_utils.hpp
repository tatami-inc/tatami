#ifndef TATAMI_ARITHMETIC_UTILS_HPP
#define TATAMI_ARITHMETIC_UTILS_HPP

#include <cmath>
#include <limits>

#include "../utils/copy.hpp"

/**
 * @file arithmetic_utils.hpp
 *
 * @brief Utilities for delayed arithmetic operations.
 */

namespace tatami {

/**
 * Type of arithmetic operation.
 *
 * The `INTEGER_DIVIDE` refers to a floored division, which differs from truncation for negative quotients.
 * This choice is based on R's `%/%`, which in turn is based on a recommendation by Donald Knuth. 
 *
 * Similarly, `x MODULO y` is defined as `x - floor(x / y)`, based on the same floored division.
 * Note that this differs from the built-in `%` operator, which performs truncation.
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
auto delayed_arithmetic(const Value_ val, const Scalar_ scalar) { 
    auto left = [&]{
        if constexpr(right_) {
            return val;
        } else {
            return scalar;
        }
    }();

    auto right = [&]{
        if constexpr(right_) {
            return scalar;
        } else {
            return val;
        }
    }();

    if constexpr(op_ == ArithmeticOperation::ADD) {
        return left + right;

    } else if constexpr(op_ == ArithmeticOperation::MULTIPLY) {
        return left * right;

    } else if constexpr(op_ == ArithmeticOperation::SUBTRACT) {
        return left - right;

    } else if constexpr(op_ == ArithmeticOperation::DIVIDE) {
        // Assume that either Value_ is an IEEE-754 float, or that division by
        // zero is impossible in this context. We don't apply manual checks
        // here to avoid performance degradation; we also don't check that the
        // other operations yield a value that doesn't overflow/underflow, so
        // it would be odd to make an exception for div-by-zero errors.
        return left / right;

    } else if constexpr(op_ == ArithmeticOperation::POWER) {
        return std::pow(left, right);

    } else if constexpr(op_ == ArithmeticOperation::MODULO) {
        // Based on a floored divide, so some work is necessary to
        // get the right value when the sign is negative.
        const auto quo = left / right; 
        if constexpr(std::numeric_limits<I<decltype(quo)> >::is_integer) {
            auto rem = left % right;
            return rem + (quo < 0 && rem != 0 ? right : 0);
        } else {
            auto rem = std::fmod(left, right);
            return rem + (quo < 0 && rem != 0 ? right : 0);
        }

    } else if constexpr(op_ == ArithmeticOperation::INTEGER_DIVIDE) {
        const auto out = left / right;
        if constexpr(std::numeric_limits<I<decltype(out)> >::is_integer) {
            // Using a floored divide. This little branch should be optimized
            // away so don't worry too much about it.
            return out - (out < 0 ? (left % right != 0) : 0);
        } else {
            return std::floor(out);
        }
    }
}

// Some of the helpers need to divide by zero to compute the fill value. The
// only way to guarantee that division by zero is supported is with IEEE
// floats, otherwise it would be undefined behavior at compile time, and the
// compiler could do anything, including refusing to compile it. So, we hide
// any '0/0' behind a constexpr to ensure that the compiler doesn't see it.
template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
constexpr bool has_unsafe_divide_by_zero() {
    typedef I<decltype(delayed_arithmetic<op_, right_>(std::declval<Value_>(), std::declval<Scalar_>()))> Product;
    if constexpr(std::numeric_limits<Product>::is_iec559) {
        return false;
    }

    // POWER also involves division by zero for negative powers,
    // but it returns an implementation-defined value; so it's "safe".
    // - https://en.cppreference.com/w/cpp/numeric/math/pow

    if constexpr(op_ == ArithmeticOperation::DIVIDE) {
        return true;
    }
    if constexpr(op_ == ArithmeticOperation::MODULO) {
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
