#ifndef TATAMI_COMPARE_UTILS_HPP
#define TATAMI_COMPARE_UTILS_HPP

/**
 * @file compare_utils.hpp
 *
 * @brief Utilities for delayed comparison operations.
 */

namespace tatami {

/**
 * Type of comparison operation.
 */
enum class CompareOperation : char { 
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
template<CompareOperation op_, typename Value_>
bool delayed_compare(Value_ val, Value_ scalar) {
    if constexpr(op_ == CompareOperation::EQUAL) {
        return val == scalar;
    } else if constexpr(op_ == CompareOperation::GREATER_THAN) {
        return val > scalar;
    } else if constexpr(op_ == CompareOperation::LESS_THAN) {
        return val < scalar;
    } else if constexpr(op_ == CompareOperation::GREATER_THAN_OR_EQUAL) {
        return val >= scalar;
    } else if constexpr(op_ == CompareOperation::LESS_THAN_OR_EQUAL) {
        return val <= scalar;
    } else { // NOT EQUAL.
        return val != scalar;
    }
}

template<CompareOperation op_, typename Value_>
void delayed_compare_run(Value_& val, Value_ scalar) {
    val = delayed_compare<op_, Value_>(val, scalar);
}
/**
 * @endcond
 */

/**
 * Type of comparison operation for special IEEE values.
 */
enum class SpecialCompareOperation : char {
    ISNAN,
    ISINF,
    ISFINITE
};

/**
 * @cond
 */
template<SpecialCompareOperation op_, bool pass_, typename Value_>
bool delayed_special_compare(Value_ val) {
    if constexpr(op_ == SpecialCompareOperation::ISNAN) {
        return pass_ == std::isnan(val);
    } else if constexpr(op_ == SpecialCompareOperation::ISINF) {
        return pass_ == std::isinf(val);
    } else {
        return pass_ == std::isfinite(val);
    }
}

template<SpecialCompareOperation op_, bool pass_, typename Value_>
void delayed_special_compare_run(Value_& val) {
    val = delayed_special_compare<op_, pass_, Value_>(val);
}
/**
 * @endcond
 */

}

#endif
