#ifndef TATAMI_BINARY_ARITH_HELPERS_H
#define TATAMI_BINARY_ARITH_HELPERS_H

#include "../arith_utils.hpp"
#include "utils.hpp"
#include <limits>
#include <vector>

/**
 * @file arith_helpers.hpp
 *
 * @brief Helper classes for binary arithmetic operations.
 */

namespace tatami {

/**
 * @brief Delayed binary arithmetic.
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOp` class.
 *
 * @tparam op_ The arithmetic operation.
 */
template<DelayedArithOp op_>
struct DelayedBinaryArithHelper {
public:
    /**
     * @cond
     */
    static constexpr bool known_sparse = (op_ == DelayedArithOp::ADD ||
                                          op_ == DelayedArithOp::SUBTRACT ||
                                          op_ == DelayedArithOp::MULTIPLY);

    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;
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
            delayed_arith_run<op_, true>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0, length = indices.size(); i < length; ++i) {
            delayed_arith_run<op_, true>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    Index_ sparse(bool, Index_, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
        // Don't bother storing an explicit zero for MULTIPLY operations when either entry is zero.
        constexpr bool must_have_both = (op_ == DelayedArithOp::MULTIPLY);
        return delayed_binary_isometric_sparse_operation<must_have_both>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](Value_& l, Value_ r) { delayed_arith_run<op_, true>(l, r); }
        );
    }

    template<typename Value_, typename Index_>
    Value_ fill(Index_) const {
        if constexpr(known_sparse) {
            return 0;
        } else if constexpr(op_ == DelayedArithOp::POWER) {
            return 1;
        } else {
            // Zero divided/modulo by zero gives NaN.
            return std::numeric_limits<Value_>::quiet_NaN();
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
