#ifndef TATAMI_BINARY_COMPARE_HELPERS_H
#define TATAMI_BINARY_COMPARE_HELPERS_H

#include "../compare_utils.hpp"
#include "utils.hpp"

/**
 * @file compare_helpers.hpp
 *
 * @brief Helper classes for binary comparison operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedBinaryIsometricOp` class.
 */

namespace tatami {

/**
 * @brief Delayed binary comparison.
 *
 * This should be used as the `OP` in the `DelayedBinaryIsometricOp` class.
 *
 * @tparam op_ The comparison operation.
 */
template<DelayedCompareOp op_>
struct DelayedBinaryCompareHelper {
public:
    /**
     * @cond
     */
    // It's sparse if f(0, 0) == 0.
    static constexpr bool always_sparse = (op_ != DelayedCompareOp::EQUAL && op_ != DelayedCompareOp::GREATER_THAN_OR_EQUAL && op_ != DelayedCompareOp::LESS_THAN_OR_EQUAL);
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            delayed_compare_run<op_>(left_buffer[i], right_buffer[i]);
        }
    }

    template<bool, bool needs_value, bool needs_index, typename Value_, typename Index_>
    Index_ sparse(Index_ idx, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer) const {
        // None of the operations will return zero if one entry is zero and the other entry is non-zero...
        // except for equality, but then, the sparse() method would never even be used.
        return delayed_binary_isometric_sparse_operation<false, needs_value, needs_index>(
            idx, 
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            [](Value_& l, Value_ r) { delayed_compare_run<op_>(l, r); }
        );
    }
    /**
     * @endcond
     */
};

/**
 * @return A helper class for a delayed binary equality comparison.
 */
inline DelayedBinaryCompareHelper<DelayedCompareOp::EQUAL> make_DelayedBinaryEqualHelper() {
    return DelayedBinaryCompareHelper<DelayedCompareOp::EQUAL>();
}

/**
 * @return A helper class for a delayed binary greater-than comparison.
 */
inline DelayedBinaryCompareHelper<DelayedCompareOp::GREATER_THAN> make_DelayedBinaryGreaterThanHelper() {
    return DelayedBinaryCompareHelper<DelayedCompareOp::GREATER_THAN>();
}

/**
 * @return A helper class for a delayed binary less-than comparison.
 */
inline DelayedBinaryCompareHelper<DelayedCompareOp::LESS_THAN> make_DelayedBinaryLessThanHelper() {
    return DelayedBinaryCompareHelper<DelayedCompareOp::LESS_THAN>();
}

/**
 * @return A helper class for a delayed binary greater-than-or-equal comparison.
 */
inline DelayedBinaryCompareHelper<DelayedCompareOp::GREATER_THAN_OR_EQUAL> make_DelayedBinaryGreaterThanOrEqualHelper() {
    return DelayedBinaryCompareHelper<DelayedCompareOp::GREATER_THAN_OR_EQUAL>();
}

/**
 * @return A helper class for a delayed binary less-than-or-equal comparison.
 */
inline DelayedBinaryCompareHelper<DelayedCompareOp::LESS_THAN_OR_EQUAL> make_DelayedBinaryLessThanOrEqualHelper() {
    return DelayedBinaryCompareHelper<DelayedCompareOp::LESS_THAN_OR_EQUAL>();
}

/**
 * @return A helper class for a delayed binary non-equality comparison to a scalar.
 */
inline DelayedBinaryCompareHelper<DelayedCompareOp::NOT_EQUAL> make_DelayedBinaryNotEqualHelper() {
    return DelayedBinaryCompareHelper<DelayedCompareOp::NOT_EQUAL>();
}

}

#endif
