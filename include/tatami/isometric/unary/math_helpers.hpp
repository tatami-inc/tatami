#ifndef TATAMI_MATH_HELPERS_H
#define TATAMI_MATH_HELPERS_H

/**
 * @file math_helpers.hpp
 *
 * @brief Helpers for unary math operations.
 */

#include <cmath>

namespace tatami {

/**
 * @brief Take the absolute value of a matrix entry.
 *
 * This class takes the absolute value of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAbs {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }

    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::abs(val);
            } else {
                output[i] = std::abs(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core<Index_>(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the sign of a matrix entry.
 *
 * This class takes the sign of each element of a `Matrix`, returning -1, 0 or 1 for negative, zero or positive values, respectively.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * This operation will report NaNs in the input as NaNs in the output if supported by the `OutputValue_` type, otherwise they are set to 0.
 *
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricSign {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                if (!std::isnan(val)) {
                    val = (static_cast<InputValue_>(0) < val) - (val < static_cast<InputValue_>(0));
                }
            } else {
                auto val = input[i];
                if (!std::isnan(val)) {
                    output[i] = (static_cast<InputValue_>(0) < val) - (val < static_cast<InputValue_>(0));
                } else if constexpr(std::numeric_limits<OutputValue_>::has_quiet_NaN) {
                    output[i] = std::numeric_limits<OutputValue_>::quiet_NaN();
                } else {
                    output[i] = 0;
                }
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core<Index_>(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry.
 *
 * This class takes the logarithm of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Base_ Numeric type for the log base.
 */
template<typename InputValue_ = double, typename Base_ = InputValue_>
class DelayedUnaryIsometricLog {
public:
    /**
     * Defaults to the natural log.
     */
    DelayedUnaryIsometricLog() : my_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedUnaryIsometricLog(Base_ base) : my_base(std::log(base)) {}

public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    Base_ my_base;

    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::log(val) / my_base;
            } else {
                output[i] = std::log(input[i]) / my_base;
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core<Index_>(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // Use the implementation-defined value.
        return std::log(static_cast<OutputValue_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the square root of a matrix entry.
 *
 * This class takes the square root of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricSqrt {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::sqrt(val);
            } else {
                output[i] = std::sqrt(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the ceiling of a matrix entry.
 *
 * This class takes the ceiling of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricCeiling {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::ceil(val);
            } else {
                output[i] = std::ceil(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core<Index_>(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the floor of a matrix entry.
 *
 * This class takes the floor of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricFloor {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::floor(val);
            } else {
                output[i] = std::floor(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Integer truncation of a matrix entry.
 *
 * This class performs an integer truncation of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricTrunc {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::trunc(val);
            } else {
                output[i] = std::trunc(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry plus 1.
 *
 * This class computes the `log1p` of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Base_ Numeric type for the log base.
 */
template<typename InputValue_ = double, typename Base_ = InputValue_>
class DelayedUnaryIsometricLog1p {
public:
    /**
     * Defaults to the natural log.
     */
    DelayedUnaryIsometricLog1p() : my_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedUnaryIsometricLog1p(Base_ base) : my_base(std::log(base)) {}

public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::log1p(val) / my_base;
            } else {
                output[i] = std::log1p(input[i]) / my_base;
            }
        }
    }

    Base_ my_base;

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Round a matrix entry to the nearest integer.
 *
 * This class rounds each element of a `Matrix` to the nearest integer.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricRound {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::round(val);
            } else {
                output[i] = std::round(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent.
 *
 * This class computes `exp()` on each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricExp {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::exp(val);
            } else {
                output[i] = std::exp(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 1;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent minus 1.
 *
 * This class computes `expm1()` on each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricExpm1 {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::expm1(val);
            } else {
                output[i] = std::expm1(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse cosine of a matrix entry.
 *
 * This class computes the inverse cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAcos {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::acos(val);
            } else {
                output[i] = std::acos(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::acos(0);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic cosine of a matrix entry.
 *
 * This class computes the inverse hyperbolic cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAcosh {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::acosh(val);
            } else {
                output[i] = std::acosh(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::acosh(static_cast<InputValue_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse sine of a matrix entry.
 *
 * This class computes the inverse sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAsin {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::asin(val);
            } else {
                output[i] = std::asin(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic sine of a matrix entry.
 *
 * This class computes the inverse hyperbolic sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAsinh {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::asinh(val);
            } else {
                output[i] = std::asinh(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse tangent of a matrix entry.
 *
 * This class computes the inverse tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAtan {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::atan(val);
            } else {
                output[i] = std::atan(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic tangent of a matrix entry.
 *
 * This class computes the inverse hyperbolic tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricAtanh {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::atanh(val);
            } else {
                output[i] = std::atanh(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the cosine of a matrix entry.
 *
 * This class computes the cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricCos {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::cos(val);
            } else {
                output[i] = std::cos(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic cosine of a matrix entry.
 *
 * This class computes the hyperbolic cosine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricCosh {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::cosh(val);
            } else {
                output[i] = std::cosh(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the sine of a matrix entry.
 *
 * This class computes the sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricSin {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::sin(val);
            } else {
                output[i] = std::sin(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic sine of a matrix entry.
 *
 * This class computes the hyperbolic sine of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricSinh {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::sinh(val);
            } else {
                output[i] = std::sinh(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the tangent of a matrix entry.
 *
 * This class computes the tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricTan {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::tan(val);
            } else {
                output[i] = std::tan(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic tangent of a matrix entry.
 *
 * This class computes the hyperbolic tangent of each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricTanh {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return true;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::tanh(val);
            } else {
                output[i] = std::tanh(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Apply the gamma function to a matrix entry.
 *
 * This class applies the gamma function to each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricGamma {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::tgamma(val);
            } else {
                output[i] = std::tgamma(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::tgamma(static_cast<InputValue_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Apply the log-gamma function to a matrix entry.
 *
 * This class applies the log-gamma function to each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 * 
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
class DelayedUnaryIsometricLgamma {
public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool is_sparse() const {
        return false;
    }
    /**
     * @endcond
     */

private:
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = std::lgamma(val);
            } else {
                output[i] = std::lgamma(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, indices.size(), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        core(input, number, output);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::lgamma(static_cast<InputValue_>(0));
    }
    /**
     * @endcond
     */
};

}

#endif
