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
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_compare<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_compare<op_>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const {
        Index_ length = indices.size();
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_compare<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_compare<op_>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    template<typename Index_, typename InputValue_, typename OutputValue_>
    Index_ sparse(bool, Index_, const SparseRange<InputValue_, Index_>& left, const SparseRange<InputValue_, Index_>& right, OutputValue_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
        // None of the operations will return zero if one entry is zero and the other entry is non-zero...
        // except for equality, but then, the sparse() method would never even be used.
        return delayed_binary_isometric_sparse_operation<false>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](InputValue_ l, InputValue_ r) -> auto { 
                return delayed_compare<op_>(l, r); 
            }
        );
    }

    template<typename OutputValue_, typename InputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
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
