#ifndef TATAMI_ISOMETRIC_BINARY_ARITHMETIC_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_ARITHMETIC_HELPERS_H

#include "../arithmetic_utils.hpp"
#include "utils.hpp"
#include "helper_interface.hpp"

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
 * @tparam Index_ Integer type for the row/column indices.
 */
template<ArithmeticOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedBinaryIsometricArithmeticHelper final : public DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    bool zero_depends_on_row() const { return false; }
    bool zero_depends_on_column() const { return false; }
    bool non_zero_depends_on_row() const { return false; }
    bool non_zero_depends_on_column() const { return false; }

public:
    void dense(
        const bool,
        const Index_,
        const Index_,
        const Index_ length,
        const InputValue_* const left_buffer,
        const InputValue_* const right_buffer,
        OutputValue_* const output_buffer)
    const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output_buffer[i];
                val = delayed_arithmetic<op_, true>(val, right_buffer[i]);
            } else {
                output_buffer[i] = delayed_arithmetic<op_, true>(left_buffer[i], right_buffer[i]);
            }
        }
    }

    void dense(
        const bool,
        const Index_,
        const std::vector<Index_>& indices,
        const InputValue_* const left_buffer,
        const InputValue_* const right_buffer,
        OutputValue_* const output_buffer)
    const {
        const Index_ length = indices.size();
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
        const bool,
        const Index_,
        const SparseRange<InputValue_, Index_>& left,
        const SparseRange<InputValue_, Index_>& right,
        OutputValue_* const value_buffer,
        Index_* const index_buffer,
        const bool needs_value,
        const bool needs_index) 
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
    OutputValue_ fill(const bool, const Index_) const {
        if constexpr(has_unsafe_divide_by_zero<op_, true, InputValue_, InputValue_>()) {
            throw std::runtime_error("division by zero is not supported");
            return 0;
        } else {
            return delayed_arithmetic<op_, true, InputValue_, InputValue_>(0, 0);
        }
    }

    bool is_sparse() const {
        return (
            op_ == ArithmeticOperation::ADD ||
            op_ == ArithmeticOperation::SUBTRACT ||
            op_ == ArithmeticOperation::MULTIPLY
        );
    }

public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }
};

/**
 * Convenient alias for the addition helper.
 *
 * @tparam OutputValue_ Type of the result of the addition.
 * @tparam InputValue_ Type of the matrix value used in the addition.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricAddHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::ADD, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the subtraction helper.
 *
 * @tparam OutputValue_ Type of the result of the subtraction.
 * @tparam InputValue_ Type of the matrix value used in the subtraction.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricSubtractHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::SUBTRACT, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the multiplication helper.
 *
 * @tparam OutputValue_ Type of the result of the multiplication.
 * @tparam InputValue_ Type of the matrix value used in the multiplication.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricMultiplyHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::MULTIPLY, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the division helper.
 *
 * @tparam OutputValue_ Type of the result of the division.
 * @tparam InputValue_ Type of the matrix value used in the division.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricDivideHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::DIVIDE, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the power operation helper.
 *
 * @tparam OutputValue_ Type of the result of the power operation.
 * @tparam InputValue_ Type of the matrix value used in the power operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricPowerHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::POWER, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the modulo operation helper.
 *
 * @tparam OutputValue_ Type of the result of the modulo operation.
 * @tparam InputValue_ Type of the matrix value used in the modulo operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricModuloHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::MODULO, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the integer division helper.
 *
 * @tparam OutputValue_ Type of the result of the integer division.
 * @tparam InputValue_ Type of the matrix value used in the integer division.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
using DelayedBinaryIsometricIntegerDivideHelper = DelayedBinaryIsometricArithmeticHelper<ArithmeticOperation::INTEGER_DIVIDE, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
// Provided for back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricAdd() {
    return std::make_shared<DelayedBinaryIsometricAddHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubtract() {
    return std::make_shared<DelayedBinaryIsometricSubtractHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricMultiply() {
    return std::make_shared<DelayedBinaryIsometricMultiplyHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricDivide() {
    return std::make_shared<DelayedBinaryIsometricDivideHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricPower() {
    return std::make_shared<DelayedBinaryIsometricPowerHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricModulo() {
    return std::make_shared<DelayedBinaryIsometricModuloHelper<OutputValue_, InputValue_, Index_> >();
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricIntegerDivide() {
    return std::make_shared<DelayedBinaryIsometricIntegerDivideHelper<OutputValue_, InputValue_, Index_> >();
}
/**
 * @endcond
 */

}

#endif
