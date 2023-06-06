#ifndef TATAMI_BOOLEAN_UTILS_HPP
#define TATAMI_BOOLEAN_UTILS_HPP

/**
 * @file boolean_utils.hpp
 *
 * @brief Utilities for delayed boolean operations.
 */

namespace tatami {

/**
 * Type of the delayed boolean operation.
 */
enum class DelayedBooleanOp : char {
    AND,
    OR,
    XOR,
    EQUAL
};

/**
 * @cond
 */
template<DelayedBooleanOp op_, typename Value_>
void delayed_boolean_run(Value_& val, bool scalar) {
    if constexpr(op_ == DelayedBooleanOp::AND) {
        val = val && scalar;
    } else if constexpr(op_ == DelayedBooleanOp::OR) {
        val = val || scalar;
    } else if constexpr(op_ == DelayedBooleanOp::XOR) {
        val = static_cast<bool>(val) != scalar;
    } else { // EQUAL.
        val = static_cast<bool>(val) == scalar;
    }
}
/**
 * @endcond
 */

}

#endif
