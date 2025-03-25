#ifndef TATAMI_ISOMETRIC_BINARY_ARITHMETIC_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_ARITHMETIC_HELPERS_H

#include "../arithmetic_utils.hpp"
#include "utils.hpp"
#include "Helper.hpp"

#include <limits>
#include <vector>

/**
 * @file arithmetic_helpers.hpp
 *
 * @brief Helper class for binary arithmetic operations.
 */

namespace tatami {

/**
 * @brief Helper for delayed binary isometric arithmetic. 
 *
 * This should be used as the `Operation_` in the `DelayedBinaryIsometricOperation` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Type of index value.
 */
template<ArithmeticOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedBinaryIsometricArithmeticHelper final : public DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    bool zero_depends_on_row() const { return false; }
    bool zero_depends_on_column() const { return false; }
    bool non_zero_depends_on_row() const { return false; }
    bool non_zero_depends_on_column() const { return false; }

public:
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

    Index_ sparse(
        bool,
        Index_,
        const SparseRange<InputValue_, Index_>& left,
        const SparseRange<InputValue_, Index_>& right,
        OutputValue_* value_buffer,
        Index_* index_buffer,
        bool needs_value,
        bool needs_index) 
    const {
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

public:
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
        return (
            op_ == ArithmeticOperation::ADD ||
            op_ == ArithmeticOperation::SUBTRACT ||
            op_ == ArithmeticOperation::MULTIPLY
        );
    }
};

/**
 * @cond
 */
// Provided for back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::ADD, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricAdd() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::ADD, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::SUBTRACT, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricSubtract() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::SUBTRACT, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::MULTIPLY, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricMultiply() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::MULTIPLY, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::DIVIDE, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricDivide() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::DIVIDE, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::POWER, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricPower() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::POWER, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::MODULO, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricModulo() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::MODULO, OutputValue_, InputValue_, Index_>();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
DelayedBinaryIsometricArithmetic<ArithmeticOperation::INTEGER_DIVIDE, OutputValue_, InputValue_, Index_> make_DelayedBinaryIsometricIntegerDivide() {
    return DelayedBinaryIsometricArithmetic<ArithmeticOperation::INTEGER_DIVIDE, OutputValue_, InputValue_, Index_>();
}
/**
 * @endcond
 */

}

#endif
