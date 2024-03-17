#ifndef TATAMI_BINARY_BOOLEAN_HELPERS_H
#define TATAMI_BINARY_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include "utils.hpp"

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper classes for binary boolean operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedBinaryIsometricOp` class.
 */

namespace tatami {

/**
 * @brief Delayed binary boolean operations.
 *
 * This should be used as the `OP` in the `DelayedBinaryIsometricOp` class.
 *
 * @tparam op_ The boolean operation.
 */
template<DelayedBooleanOp op_>
struct DelayedBinaryBooleanHelper {
public:
    /**
     * @cond
     */
    // It's sparse if f(0, 0) == 0.
    static constexpr bool is_sparse = (op_ != DelayedBooleanOp::EQUAL);
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
        constexpr bool must_have_both = (op_ == DelayedBooleanOp::AND);
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
    /**
     * @endcond
     */
};

/**
 * @return A helper class for a delayed binary boolean equivalence operation.
 */
inline DelayedBinaryBooleanHelper<DelayedBooleanOp::EQUAL> make_DelayedBinaryBooleanEqualHelper() {
    return DelayedBinaryBooleanHelper<DelayedBooleanOp::EQUAL>();
}

/**
 * @return A helper class for a delayed binary AND comparison.
 */
inline DelayedBinaryBooleanHelper<DelayedBooleanOp::AND> make_DelayedBinaryBooleanAndHelper() {
    return DelayedBinaryBooleanHelper<DelayedBooleanOp::AND>();
}

/**
 * @return A helper class for a delayed binary OR comparison.
 */
inline DelayedBinaryBooleanHelper<DelayedBooleanOp::OR> make_DelayedBinaryBooleanOrHelper() {
    return DelayedBinaryBooleanHelper<DelayedBooleanOp::OR>();
}

/**
 * @return A helper class for a delayed binary XOR comparison.
 */
inline DelayedBinaryBooleanHelper<DelayedBooleanOp::XOR> make_DelayedBinaryBooleanXorHelper() {
    return DelayedBinaryBooleanHelper<DelayedBooleanOp::XOR>();
}

}

#endif
