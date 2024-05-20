#ifndef TATAMI_ISOMETRIC_BINARY_COMPARE_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_COMPARE_HELPERS_H

#include "../compare_utils.hpp"
#include "utils.hpp"

/**
 * @file compare_helpers.hpp
 *
 * @brief Helper classes for binary comparison operations.
 */

namespace tatami {

/**
 * @brief Delayed binary isometric comparison.
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 */
template<CompareOperation op_>
struct DelayedBinaryIsometricCompare {
public:
    /**
     * @cond
     */
    // It's sparse if f(0, 0) == 0.
    static constexpr bool known_sparse = (op_ != CompareOperation::EQUAL && 
                                          op_ != CompareOperation::GREATER_THAN_OR_EQUAL && 
                                          op_ != CompareOperation::LESS_THAN_OR_EQUAL);

    static constexpr bool is_basic = false;
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<typename Value_, typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            delayed_compare_run<op_>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0, length = indices.size(); i < length; ++i) {
            delayed_compare_run<op_>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    Index_ sparse(bool, Index_, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
        // None of the operations will return zero if one entry is zero and the other entry is non-zero...
        // except for equality, but then, the sparse() method would never even be used.
        return delayed_binary_isometric_sparse_operation<false>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](Value_& l, Value_ r) { delayed_compare_run<op_>(l, r); }
        );
    }

    template<typename Value_, typename Index_>
    Index_ sparse(bool r, Index_ i, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer) const {
        return sparse<Value_, Index_>(r, i, left, right, value_buffer, index_buffer, true, true);
    }

    template<typename Value_, typename Index_>
    Value_ fill(bool, Index_) const {
        if constexpr(known_sparse) {
            return 0;
        } else {
            return 1;
        }
    }

    bool is_sparse() const {
        return known_sparse;
    }
    /**
     * @endcond
     */
};

/**
 * @return A helper class for a delayed binary equality comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricCompare<CompareOperation::EQUAL> make_DelayedBinaryIsometricEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::EQUAL>();
}

/**
 * @return A helper class for a delayed binary greater-than comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN> make_DelayedBinaryIsometricGreaterThan() {
    return DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN>();
}

/**
 * @return A helper class for a delayed binary less-than comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN> make_DelayedBinaryIsometricLessThan() {
    return DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN>();
}

/**
 * @return A helper class for a delayed binary greater-than-or-equal comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN_OR_EQUAL> make_DelayedBinaryIsometricGreaterThanOrEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN_OR_EQUAL>();
}

/**
 * @return A helper class for a delayed binary less-than-or-equal comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN_OR_EQUAL> make_DelayedBinaryIsometricLessThanOrEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN_OR_EQUAL>();
}

/**
 * @return A helper class for a delayed binary non-equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricCompare<CompareOperation::NOT_EQUAL> make_DelayedBinaryIsometricNotEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::NOT_EQUAL>();
}

}

#endif
