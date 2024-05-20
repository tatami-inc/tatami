#ifndef TATAMI_ISOMETRIC_BINARY_ARITHMETIC_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_ARITHMETIC_HELPERS_H

#include "../arithmetic_utils.hpp"
#include "utils.hpp"
#include <limits>
#include <vector>

/**
 * @file arithmetic_helpers.hpp
 *
 * @brief Helper classes for binary arithmetic operations.
 */

namespace tatami {

/**
 * @brief Delayed binary isometric arithmetic. 
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOperation` class.
 *
 * @tparam op_ The arithmetic operation.
 */
template<ArithmeticOperation op_>
class DelayedBinaryIsometricArithmetic {
public:
    /**
     * @cond
     */
    static constexpr bool known_sparse = (op_ == ArithmeticOperation::ADD ||
                                          op_ == ArithmeticOperation::SUBTRACT ||
                                          op_ == ArithmeticOperation::MULTIPLY);

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
            delayed_arithmetic_run<op_, true>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* left_buffer, const Value_* right_buffer) const {
        for (Index_ i = 0, length = indices.size(); i < length; ++i) {
            delayed_arithmetic_run<op_, true>(left_buffer[i], right_buffer[i]);
        }
    }

    template<typename Value_, typename Index_>
    Index_ sparse(bool, Index_, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
        // Don't bother storing an explicit zero for MULTIPLY operations when either entry is zero.
        constexpr bool must_have_both = (op_ == ArithmeticOperation::MULTIPLY);
        return delayed_binary_isometric_sparse_operation<must_have_both>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](Value_& l, Value_ r) { delayed_arithmetic_run<op_, true>(l, r); }
        );
    }

    template<typename Value_, typename Index_>
    Value_ fill(bool, Index_) const {
        if constexpr(known_sparse) {
            return 0;
        } else if constexpr(op_ == ArithmeticOperation::POWER) {
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
 * @return A helper class for delayed binary addition,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::ADD> make_DelayedBinaryIsometricAdd() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::ADD>();
}

/**
 * @return A helper class for delayed binary subtraction,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::SUBTRACT> make_DelayedBinaryIsometricSubtract() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::SUBTRACT>();
}

/**
 * @return A helper class for delayed binary multiplication,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::MULTIPLY> make_DelayedBinaryIsometricMultiply() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::MULTIPLY>();
}

/**
 * @return A helper class for delayed binary division,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::DIVIDE> make_DelayedBinaryIsometricDivide() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::DIVIDE>();
}

/**
 * @return A helper class for delayed binary power,
 * to be used as the `operation` in a `DelayedBinaryIsometricOperation`.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::POWER> make_DelayedBinaryIsometricPower() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::POWER>();
}

/**
 * @return A helper class for delayed binary modulo.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::MODULO> make_DelayedBinaryIsometricModulo() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::MODULO>();
}

/**
 * @return A helper class for delayed binary integer division.
 */
inline DelayedBinaryIsometricArithmetic<ArithmeticOperation::INTEGER_DIVIDE> make_DelayedBinaryIsometricIntegerDivide() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::INTEGER_DIVIDE>();
}

}

#endif
