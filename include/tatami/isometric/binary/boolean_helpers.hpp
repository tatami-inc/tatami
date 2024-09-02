#ifndef TATAMI_ISOMETRIC_BINARY_BOOLEAN_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include "utils.hpp"

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper classes for binary boolean operations.
 */

namespace tatami {

/**
 * @brief Delayed binary isometric boolean operations.
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOperation` class.
 *
 * @tparam op_ The boolean operation.
 */
template<BooleanOperation op_>
struct DelayedBinaryIsometricBoolean {
public:
    /**
     * @cond
     */
    // It's sparse if f(0, 0) == 0.
    static constexpr bool known_sparse = (op_ != BooleanOperation::EQUAL);

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
                val = delayed_boolean<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_boolean<op_>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const {
        Index_ length = indices.size();
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_boolean<op_>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_boolean<op_>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    template<typename Index_, typename InputValue_, typename OutputValue_>
    Index_ sparse(bool, Index_, const SparseRange<InputValue_, Index_>& left, const SparseRange<InputValue_, Index_>& right, OutputValue_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
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
 * @return A helper class for a delayed binary boolean equivalence operation,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricBoolean<BooleanOperation::EQUAL> make_DelayedBinaryIsometricBooleanEqual() {
    return DelayedBinaryIsometricBoolean<BooleanOperation::EQUAL>();
}

/**
 * @return A helper class for a delayed binary AND comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricBoolean<BooleanOperation::AND> make_DelayedBinaryIsometricBooleanAnd() {
    return DelayedBinaryIsometricBoolean<BooleanOperation::AND>();
}

/**
 * @return A helper class for a delayed binary OR comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricBoolean<BooleanOperation::OR> make_DelayedBinaryIsometricBooleanOr() {
    return DelayedBinaryIsometricBoolean<BooleanOperation::OR>();
}

/**
 * @return A helper class for a delayed binary XOR comparison,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricBoolean<BooleanOperation::XOR> make_DelayedBinaryIsometricBooleanXor() {
    return DelayedBinaryIsometricBoolean<BooleanOperation::XOR>();
}

}

#endif
