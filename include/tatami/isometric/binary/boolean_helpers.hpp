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
    template<typename Value_, typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            delayed_boolean_run<op_>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0, length = indices.size(); i < length; ++i) {
            delayed_boolean_run<op_>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    Index_ sparse(bool, Index_, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
        // Don't bother storing an explicit zero for AND operations when either entry is zero.
        constexpr bool must_have_both = (op_ == BooleanOperation::AND);
        return delayed_binary_isometric_sparse_operation<must_have_both>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](Value_& l, Value_ r) { delayed_boolean_run<op_>(l, r); }
        );
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
