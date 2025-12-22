#ifndef TATAMI_MATH_HELPERS_H
#define TATAMI_MATH_HELPERS_H

/**
 * @file math_helpers.hpp
 *
 * @brief Helpers for unary math operations.
 */

#include "helper_interface.hpp"
#include <cmath>

namespace tatami {

/**
 * @brief Helper for delayed calculation of the sign of each matrix entry.
 *
 * This class takes the absolute value of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAbsHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::abs(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the sign of each matrix entry.
 *
 * This class takes the sign of each element of a `Matrix`, returning -1, 0 or 1 for negative, zero or positive values, respectively.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * This operation will report NaNs in the input as NaNs in the output if supported by the `OutputValue_` type, otherwise they are set to 0.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSignHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            auto val = input[i];
            if ([&]{
                if constexpr(std::numeric_limits<InputValue_>::has_quiet_NaN) {
                    return !std::isnan(val);
                } else {
                    return true;
                }
            }()) {
                output[i] = (static_cast<InputValue_>(0) < val) - (val < static_cast<InputValue_>(0));
            } else if constexpr(std::numeric_limits<OutputValue_>::has_quiet_NaN) {
                output[i] = std::numeric_limits<OutputValue_>::quiet_NaN();
            } else {
                output[i] = 0;
            }
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of a logarithm of each matrix entry.
 *
 * This class takes the logarithm of each element of a `Matrix` with a custom base.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam base_ Base of the logarithm.
 * This can either be 2, 10 or -1 (for the natural log).
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<int base_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricFixedLogHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    static auto logify(InputValue_ x) {
        if constexpr(base_ == 2) {
            return std::log2(x);
        } else if constexpr(base_ == 10) {
            return std::log10(x);
        } else {
            return std::log(x);
        }
    }

    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = logify(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // Use the implementation-defined value.
        return logify(0);
    }
};
/**
 * @endcond
 */

/**
 * Convenient alias for a delayed natural logarithm.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricLogHelper = DelayedUnaryIsometricFixedLogHelper<-1, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for a delayed logarithm with base 2.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricLog2Helper = DelayedUnaryIsometricFixedLogHelper<2, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for a delayed logarithm with base 10.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricLog10Helper = DelayedUnaryIsometricFixedLogHelper<10, OutputValue_, InputValue_, Index_>;

/**
 * @brief Helper for delayed calculation of a custom logarithm of each matrix entry.
 *
 * This class takes the logarithm of each element of a `Matrix` with a custom base.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Base_ Type of the base of the logarithm.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Base_>
class DelayedUnaryIsometricCustomLogHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param base Base of the logarithm.
     */
    DelayedUnaryIsometricCustomLogHelper(Base_ base) : my_base(std::log(base)) {}

public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    Base_ my_base;

    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::log(input[i]) / my_base;
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // Use the implementation-defined value.
        return std::log(static_cast<InputValue_>(0)) / my_base;
    }
};

/**
 * @brief Helper for delayed calculation of the square root of each matrix entry.
 *
 * This class takes the square root of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSqrtHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::sqrt(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the ceiling of each matrix entry.
 *
 * This class takes the ceiling of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricCeilingHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::ceil(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the floor of each matrix entry.
 *
 * This class takes the floor of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricFloorHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::floor(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed truncation of each matrix entry to an integer.
 *
 * This class performs an integer truncation of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricTruncHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::trunc(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for the delayed calculation of the logarithm of each matrix entry plus 1.
 *
 * This class computes the `log1p` of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricLog1pHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::log1p(input[i]); 
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed rounding of each matrix entry to the nearest integer.
 *
 * This class rounds each element of a `Matrix` to the nearest integer.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricRoundHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::round(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the exponent function for each matrix entry.
 *
 * This class computes `exp()` on each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricExpHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::exp(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 1;
    }
};

/**
 * @brief Helper for delayed calculation of the exponential function of each matrix entry minus 1.
 *
 * This class computes `expm1()` on each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricExpm1Helper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::expm1(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the inverse cosine of each matrix entry.
 *
 * This class computes the inverse cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAcosHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::acos(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // Use the implementation-defined special value.
        return std::acos(0);
    }
};

/**
 * @brief Helper for delayed calculation of the inverse hyperbolic cosine of each matrix entry.
 *
 * This class computes the inverse hyperbolic cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAcoshHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::acosh(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // Use the implementation-defined special value.
        return std::acosh(static_cast<InputValue_>(0));
    }
};

/**
 * @brief Helper for delayed calculation of the inverse sine of each matrix entry.
 *
 * This class computes the inverse sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAsinHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_>  {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::asin(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the inverse hyperbolic sine of each matrix entry.
 *
 * This class computes the inverse hyperbolic sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAsinhHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::asinh(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the inverse tangent of each matrix entry.
 *
 * This class computes the inverse tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAtanHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::atan(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the inverse hyperbolic tangent of each matrix entry.
 *
 * This class computes the inverse hyperbolic tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricAtanhHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::atanh(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the cosine of a matrix entry.
 *
 * This class computes the cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricCosHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::cos(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 1;
    }
};

/**
 * @brief Helper for delayed calculation of the hyperbolic cosine of each matrix entry.
 *
 * This class computes the hyperbolic cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricCoshHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::cosh(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 1;
    }
};

/**
 * @brief Helper for delayed calculation of the sine of each matrix entry.
 *
 * This class computes the sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSinHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::sin(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the hyperbolic sine of each matrix entry.
 *
 * This class computes the hyperbolic sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSinhHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::sinh(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the tangent of each matrix entry.
 *
 * This class computes the tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricTanHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::tan(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Helper for delayed calculation of the hyperbolic tangent of each matrix entry.
 *
 * This class computes the hyperbolic tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricTanhHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::tanh(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        return 0;
    }
};

/**
 * @brief Apply the gamma function to a matrix entry.
 *
 * This class applies the gamma function to each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricGammaHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::tgamma(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // Use the implementation-defined special value.
        return std::tgamma(static_cast<InputValue_>(0));
    }
};

/**
 * @brief Apply the log-gamma function to a matrix entry.
 *
 * This class applies the log-gamma function to each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricLgammaHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

private:
    void core(const InputValue_* input, const Index_ length, OutputValue_* const output) const {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            input = output; // basically an assertion to the compiler to allow it to skip aliasing protection.
        }
        for (Index_ i = 0; i < length; ++i) {
            output[i] = std::lgamma(input[i]);
        }
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        core(input, length, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        core(input, indices.size(), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input, const Index_* const, OutputValue_* const output) const {
        core(input, number, output);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // Use the implementation-defined special value.
        return std::lgamma(static_cast<InputValue_>(0));
    }
};

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAbs = DelayedUnaryIsometricAbsHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricSign = DelayedUnaryIsometricSignHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Base_ = InputValue_> 
using DelayedUnaryIsometricLog = DelayedUnaryIsometricCustomLogHelper<OutputValue_, InputValue_, Index_, Base_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricSqrt = DelayedUnaryIsometricSqrtHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricCeiling = DelayedUnaryIsometricCeilingHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricFloor = DelayedUnaryIsometricFloorHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricTrunc = DelayedUnaryIsometricTruncHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Base_ = InputValue_> 
using DelayedUnaryIsometricLog1p = DelayedUnaryIsometricLog1pHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricRound = DelayedUnaryIsometricRoundHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricExp = DelayedUnaryIsometricExpHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricExpm1 = DelayedUnaryIsometricExpm1Helper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAcos = DelayedUnaryIsometricAcosHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAcosh = DelayedUnaryIsometricAcoshHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAsin = DelayedUnaryIsometricAsinHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAsinh = DelayedUnaryIsometricAsinhHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAtan = DelayedUnaryIsometricAtanHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricAtanh = DelayedUnaryIsometricAtanhHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricCos = DelayedUnaryIsometricCosHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricCosh = DelayedUnaryIsometricCoshHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricSin = DelayedUnaryIsometricSinHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricSinh = DelayedUnaryIsometricSinhHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricTan = DelayedUnaryIsometricTanHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricTanh = DelayedUnaryIsometricTanhHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricGamma = DelayedUnaryIsometricGammaHelper<OutputValue_, InputValue_, Index_>;

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int> 
using DelayedUnaryIsometricLgamma = DelayedUnaryIsometricLgammaHelper<OutputValue_, InputValue_, Index_>;
/**
 * @endcond
 */

}

#endif
