#ifndef TATAMI_BINARY_ARITH_HELPERS_H
#define TATAMI_BINARY_ARITH_HELPERS_H

#include "../arith_utils.hpp"

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
    static constexpr bool always_sparse = (op_ != DelayedArithOp::DIVIDE);
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
        Index_ lcount = 0, rcount = 0, output = 0;
        constexpr bool must_have_both = (op_ == DelayedArithOp::MULTIPLY);

        auto advance_left = [&]() -> void {
            if constexpr(needs_value) {
                value_buffer[output] = left.value[lcount];
                delayed_arith_run<op_, true>(value_buffer[output], 0);
            }
            if constexpr(needs_index) {
                index_buffer[output] = left.index[lcount];
            }
            ++output;
            ++lcount;
        };

        auto advance_right = [&]() -> void {
            if constexpr(needs_value) {
                value_buffer[output] = 0;
                delayed_arith_run<op_, true>(value_buffer[output], right.value[rcount]);
            }
            if constexpr(needs_index) {
                index_buffer[output] = right.index[rcount];
            }
            ++rcount;
            ++output;
        };

        while (lcount < left.number && rcount < right.number) {
            if (left.index[lcount] < right.index[rcount]) {
                if constexpr(!must_have_both) {
                    advance_left();
                } else {
                    ++lcount;
                }

            } else if (left.index[lcount] > right.index[rcount]) {
                if constexpr(!must_have_both) {
                    advance_right();
                } else {
                    ++rcount;
                }

            } else {
                if constexpr(needs_value) {
                    value_buffer[output] = left.value[lcount];
                    delayed_arith_run<op_, true>(value_buffer[output], right.value[rcount]);
                }
                if constexpr(needs_index) {
                    index_buffer[output] = right.index[rcount];
                }
                ++lcount;
                ++rcount;
                ++output;
            }
        }

        if constexpr(!must_have_both) {
            while (lcount < left.number) {
                advance_left();
            }

            while (rcount < right.number) {
                advance_right();
            }
        }

        return output;
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

}

#endif
