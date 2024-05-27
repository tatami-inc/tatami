#ifndef TATAMI_BOOLEAN_UTILS_HPP
#define TATAMI_BOOLEAN_UTILS_HPP

/**
 * @file boolean_utils.hpp
 *
 * @brief Utilities for delayed boolean operations.
 */

namespace tatami {

/**
 * Type of boolean operation.
 */
enum class BooleanOperation : char {
    AND,
    OR,
    XOR,
    EQUAL
};

/**
 * @cond
 */
template<BooleanOperation op_, typename Value_>
bool delayed_boolean(Value_ val, bool scalar) {
    if constexpr(op_ == BooleanOperation::AND) {
        return val && scalar;
    } else if constexpr(op_ == BooleanOperation::OR) {
        return val || scalar;
    } else if constexpr(op_ == BooleanOperation::XOR) {
        return static_cast<bool>(val) != scalar;
    } else { // EQUAL.
        return static_cast<bool>(val) == scalar;
    }
}

template<BooleanOperation op_, typename Value_>
void delayed_boolean_run(Value_& val, bool scalar) {
    val = delayed_boolean<op_>(val, scalar);
}
/**
 * @endcond
 */

}

#endif
