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
    static constexpr bool always_sparse = (op_ != DelayedBooleanOp::EQUAL);
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
            delayed_boolean_run<op_>(left_buffer[i], right_buffer[i]);
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
