#ifndef TATAMI_COMPARE_UTILS_HPP
#define TATAMI_COMPARE_UTILS_HPP

/**
 * @file compare_utils.hpp
 *
 * @brief Utilities for delayed comparison operations.
 */

namespace tatami {

/**
 * Type of the delayed comparison operation.
 */
enum class DelayedCompareOp : char { 
    EQUAL, 
    GREATER_THAN, 
    LESS_THAN, 
    GREATER_THAN_OR_EQUAL, 
    LESS_THAN_OR_EQUAL, 
    NOT_EQUAL
};

/**
 * @cond
 */
template<DelayedCompareOp op_, typename Scalar_, typename Value_>
void delayed_compare_run(Value_& val, Scalar_ scalar) {
    if constexpr(op_ == DelayedCompareOp::EQUAL) {
        val = val == scalar;
    } else if constexpr(op_ == DelayedCompareOp::GREATER_THAN) {
        val = val > scalar;
    } else if constexpr(op_ == DelayedCompareOp::LESS_THAN) {
        val = val < scalar;
    } else if constexpr(op_ == DelayedCompareOp::GREATER_THAN_OR_EQUAL) {
        val = val >= scalar;
    } else if constexpr(op_ == DelayedCompareOp::LESS_THAN_OR_EQUAL) {
        val = val <= scalar;
    } else { // NOT EQUAL.
        val = val != scalar;
    }
}
/**
 * @endcond
 */

}

#endif
