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
 * @tparam Index_ Integer type for the row/column indices.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
struct DelayedBinaryIsometricCompareHelper final : public DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    bool zero_depends_on_row() const { return false; }
    bool zero_depends_on_column() const { return false; }
    bool non_zero_depends_on_row() const { return false; }
    bool non_zero_depends_on_column() const { return false; }

public:
    void dense(
        const bool,
        const Index_,
        const Index_,
        const Index_ length,
        const InputValue_* const left_buffer,
        const InputValue_* const right_buffer,
        OutputValue_* const output_buffer)
    const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_compare<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_compare<op_>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    void dense(
        const bool,
        const Index_,
        const std::vector<Index_>& indices,
        const InputValue_* const left_buffer,
        const InputValue_* const right_buffer,
        OutputValue_* const output_buffer)
    const {
        const Index_ length = indices.size();
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
        const bool,
        const Index_,
        const SparseRange<InputValue_, Index_>& left,
        const SparseRange<InputValue_, Index_>& right,
        OutputValue_* const value_buffer,
        Index_* const index_buffer,
        const bool needs_value,
        const bool needs_index)
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

    OutputValue_ fill(const bool, const Index_) const {
        if constexpr(known_sparse) {
            return 0;
        } else {
            return 1;
        }
    }

    bool is_sparse() const {
        return known_sparse;
    }

public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }
};

/**
 * Convenient alias for the equality comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricEqualHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "greater than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricGreaterThanHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "greater than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricGreaterThanHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "less than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricLessThanHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "greater than or equal" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricGreaterThanOrEqualHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "less than or equal" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricLessThanOrEqualHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the non-equality comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricNotEqualHelper = DelayedBinaryIsometricCompareHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricEqual() {
    return std::make_shared<DelayedBinaryIsometricEqualHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricGreaterThan() {
    return std::make_shared<DelayedBinaryIsometricGreaterThanHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricLessThan() {
    return std::make_shared<DelayedBinaryIsometricLessThanHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricGreaterThanOrEqual() {
    return std::make_shared<DelayedBinaryIsometricGreaterThanOrEqualHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricLessThanOrEqual() {
    return std::make_shared<DelayedBinaryIsometricLessThanOrEqualHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricNotEqual() {
    return std::make_shared<DelayedBinaryIsometricNotEqualHelper<OutputValue_, InputValue_, Index_> >();
}
/**
 * @endcond
 */

}

#endif
