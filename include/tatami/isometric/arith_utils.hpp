#ifndef TATAMI_ARITH_UTILS_HPP
#define TATAMI_ARITH_UTILS_HPP

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
    DIVIDE
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
    } else {
        // Assume the either Value_ is an IEEE-754 float, or the user has
        // already sanitized the scalar to avoid any division by zero.
        if constexpr(right_) {
            val /= scalar;
        } else {
            val = scalar / val;
        }
    }
}
/**
 * @endcond
 */
}

#endif
