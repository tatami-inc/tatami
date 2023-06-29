#ifndef TATAMI_BINARY_ARITH_HELPERS_H
#define TATAMI_BINARY_ARITH_HELPERS_H

#include "../arith_utils.hpp"
#include "utils.hpp"

/**
 * @file arith_helpers.hpp
 *
 * @brief Helper classes for binary arithmetic operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedBinaryIsometricOp` class.
 */

namespace tatami {

/**
 * @brief Delayed binary arithmetic.
 *
 * This should be used as the `OP` in the `DelayedBinaryIsometricOp` class.
 *
 * @tparam op_ The arithmetic operation.
 */
template<DelayedArithOp op_>
struct DelayedBinaryArithHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_sparse = (op_ == DelayedArithOp::ADD ||
                                           op_ == DelayedArithOp::SUBTRACT ||
                                           op_ == DelayedArithOp::MULTIPLY);
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
            delayed_arith_run<op_, true>(left_buffer[i], right_buffer[i]);
        }
    }

    template<bool, bool needs_value, bool needs_index, typename Value_, typename Index_>
    Index_ sparse(Index_ idx, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer) const {
        // Don't bother storing an explicit zero for MULTIPLY operations when either entry is zero.
        constexpr bool must_have_both = (op_ == DelayedArithOp::MULTIPLY);
        return delayed_binary_isometric_sparse_operation<must_have_both, needs_value, needs_index>(
            idx, 
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            [](Value_& l, Value_ r) { delayed_arith_run<op_, true>(l, r); }
        );
    }
    /**
     * @endcond
     */
};

/**
 * @return A helper class for delayed binary addition.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::ADD> make_DelayedBinaryAddHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::ADD>();
}

/**
 * @return A helper class for delayed binary subtraction.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::SUBTRACT> make_DelayedBinarySubtractHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::SUBTRACT>();
}

/**
 * @return A helper class for delayed binary multiplication.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::MULTIPLY> make_DelayedBinaryMultiplyHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::MULTIPLY>();
}

/**
 * @return A helper class for delayed binary division.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::DIVIDE> make_DelayedBinaryDivideHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::DIVIDE>();
}

/**
 * @return A helper class for delayed binary power.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::POWER> make_DelayedBinaryPowerHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::POWER>();
}

/**
 * @return A helper class for delayed binary modulo.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::MODULO> make_DelayedBinaryModuloHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::MODULO>();
}

/**
 * @return A helper class for delayed binary integer division.
 */
inline DelayedBinaryArithHelper<DelayedArithOp::INTEGER_DIVIDE> make_DelayedBinaryIntegerDivideHelper() {
    return DelayedBinaryArithHelper<DelayedArithOp::INTEGER_DIVIDE>();
}

}

#endif
