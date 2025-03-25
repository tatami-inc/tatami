#ifndef TATAMI_ISOMETRIC_BINARY_COMPARE_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_COMPARE_HELPERS_H

#include "../compare_utils.hpp"
#include "utils.hpp"
#include "helper_interface.hpp"

/**
 * @file compare_helpers.hpp
 *
 * @brief Helper class for binary comparison operations.
 */

namespace tatami {

/**
 * @brief Helper for delayed binary isometric comparisons.
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Type of index value.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
struct DelayedBinaryIsometricCompare final : public DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    bool zero_depends_on_row() const { return false; }
    bool zero_depends_on_column() const { return false; }
    bool non_zero_depends_on_row() const { return false; }
    bool non_zero_depends_on_column() const { return false; }

public:
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

    Index_ sparse(
        bool,
        Index_,
        const SparseRange<InputValue_, Index_>& left,
        const SparseRange<InputValue_, Index_>& right,
        OutputValue_* value_buffer,
        Index_* index_buffer,
        bool needs_value,
        bool needs_index)
    const {
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

public:
    /**
     * @cond
     */
    // It's sparse if f(0, 0) == 0.
    static constexpr bool known_sparse = (op_ != CompareOperation::EQUAL && 
                                          op_ != CompareOperation::GREATER_THAN_OR_EQUAL && 
                                          op_ != CompareOperation::LESS_THAN_OR_EQUAL);
    /**
     * @endcond
     */

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
};

/**
 * @cond
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricCompare<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricGreaterThan() {
    return DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricLessThan() {
    return DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricGreaterThanOrEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricLessThanOrEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricCompare<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricNotEqual() {
    return DelayedBinaryIsometricCompare<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_>();
}
/**
 * @endcond
 */

}

#endif
