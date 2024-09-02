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
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_arithmetic<op_, true>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_arithmetic<op_, true>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const {
        Index_ length = indices.size();
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_arithmetic<op_, true>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_arithmetic<op_, true>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    template<typename Index_, typename InputValue_, typename OutputValue_>
    Index_ sparse(bool, Index_, const SparseRange<InputValue_, Index_>& left, const SparseRange<InputValue_, Index_>& right, OutputValue_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const {
        // Technically, MULTIPLY could skip processing if either is a zero.
        // However, this is not correct if the other value is an NaN/Inf, as
        // the product would be a NaN, not a zero; so we have to err on the
        // side of caution of attemping the operation.
        constexpr bool must_have_both = (op_ == ArithmeticOperation::MULTIPLY && 
                                         !std::numeric_limits<InputValue_>::has_quiet_NaN && 
                                         !std::numeric_limits<InputValue_>::has_infinity);

        return delayed_binary_isometric_sparse_operation<must_have_both>(
            left, 
            right, 
            value_buffer, 
            index_buffer, 
            needs_value,
            needs_index,
            [](InputValue_ l, InputValue_ r) -> auto { 
                return delayed_arithmetic<op_, true>(l, r); 
            }
        );
    }

    template<typename OutputValue_, typename InputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        if constexpr(has_unsafe_divide_by_zero<op_, true, InputValue_, InputValue_>()) {
            throw std::runtime_error("division by zero is not supported");
            return 0;
        } else {
            return delayed_arithmetic<op_, true, InputValue_>(0, 0);
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
