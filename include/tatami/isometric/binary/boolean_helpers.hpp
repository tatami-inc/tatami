#ifndef TATAMI_ISOMETRIC_BINARY_BOOLEAN_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include "utils.hpp"
#include "helper_interface.hpp"

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper class for binary boolean operations.
 */

namespace tatami {

/**
 * @brief Helper for delayed binary isometric boolean operations.
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOperation` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<BooleanOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
struct DelayedBinaryIsometricBooleanHelper final : public DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
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
                val = delayed_boolean<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_boolean<op_>(left_buffer[i], right_buffer[i]);
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
                val = delayed_boolean<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_boolean<op_>(left_buffer[i], right_buffer[i]);
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
        // Don't bother storing an explicit zero for AND operations when either
        // entry is zero. This should be NaN-safe as NaNs are truthy, so
        // applying AND on that would just be false anyway.
        constexpr bool must_have_both = (op_ == BooleanOperation::AND);
        return delayed_binary_isometric_sparse_operation<must_have_both>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](InputValue_ l, InputValue_ r) -> auto { 
                return delayed_boolean<op_>(l, r); 
            }
        );
    }

public:
    /**
     * @cond
     */
    // It's sparse if f(0, 0) == 0.
    static constexpr bool known_sparse = (op_ != BooleanOperation::EQUAL);
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
 * Convenient alias for the boolean equality helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricBooleanEqualHelper = DelayedBinaryIsometricBooleanHelper<BooleanOperation::EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the boolean AND helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricBooleanAndHelper = DelayedBinaryIsometricBooleanHelper<BooleanOperation::AND, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the boolean OR helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricBooleanOrHelper = DelayedBinaryIsometricBooleanHelper<BooleanOperation::OR, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the boolean XOR helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricBooleanXorHelper = DelayedBinaryIsometricBooleanHelper<BooleanOperation::XOR, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricBooleanEqual() {
    return std::make_shared<DelayedBinaryIsometricBooleanEqualHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricBooleanAnd() {
    return std::make_shared<DelayedBinaryIsometricBooleanAndHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricBooleanOr() {
    return std::make_shared<DelayedBinaryIsometricBooleanOrHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricBooleanXor() {
    return std::make_shared<DelayedBinaryIsometricBooleanXorHelper<OutputValue_, InputValue_, Index_> >();
}
/**
 * @endcond
 */

}

#endif
