#ifndef TATAMI_MATH_HELPERS_H
#define TATAMI_MATH_HELPERS_H

/**
 * @file math_helpers.hpp
 *
 * @brief Helpers for unary math operations.
 *
 * These classes should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 */

#include <cmath>

namespace tatami {

/**
 * @brief Take the absolute value of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAbsHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::abs(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the sign of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedSignHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            if (!std::isnan(buffer[i])) {
                buffer[i] = (static_cast<Value_>(0) < buffer[i]) - (buffer[i] < static_cast<Value_>(0));
            }
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry.
 * 
 * @tparam Value_ Type of the data value.
 * @tparam Base_ Numeric type for the log base.
 */
template<typename Value_ = double, typename Base_ = Value_>
struct DelayedLogHelper {
    /**
     * Defaults to the natural log.
     */
    DelayedLogHelper() : log_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedLogHelper(Base_ base) : log_base(std::log(base)) {}

public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static const bool needs_row = false;

    static const bool needs_column = false;
    /**
     * @endcond
     */

private:
    const Base_ log_base;

    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::log(buffer[i]) / log_base;
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        // Use the implementation-defined value.
        return std::log(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the square root of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedSqrtHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::sqrt(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the ceiling of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedCeilingHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::ceil(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the floor of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedFloorHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::floor(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the trunc of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedTruncHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::trunc(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry plus 1.
 *
 * @tparam Value_ Type of the data value.
 * @tparam Base_ Numeric type for the log base.
 */
template<typename Value_ = double, typename Base_ = Value_>
struct DelayedLog1pHelper {
    /**
     * Defaults to the natural log.
     */
    DelayedLog1pHelper() : log_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedLog1pHelper(Base_ base) : log_base(std::log(base)) {}

public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::log1p(buffer[i]) / log_base;
        }
    }

    const Base_ log_base;

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Round a matrix entry to the nearest integer.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedRoundHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::round(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedExpHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::exp(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent minus 1.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedExpm1Helper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::expm1(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the arc cosine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAcosHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::acos(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return std::acos(0);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic cosine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAcoshHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::acosh(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return std::acosh(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the arc sine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAsinHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::asin(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic sine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAsinhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::asinh(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the arc tangent of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAtanHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::atan(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic tangent of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedAtanhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::atanh(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the cosine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedCosHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::cos(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic cosine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedCoshHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::cosh(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the sine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedSinHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::sin(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic sine of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedSinhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::sinh(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the tangent of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedTanHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::tan(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic tangent of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedTanhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = true;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::tanh(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the gamma of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedGammaHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::tgamma(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return std::tgamma(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of the gamma of a matrix entry.
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedLgammaHelper {
public:
    /**
     * @cond
     */
    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;

    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::lgamma(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        return std::lgamma(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

}

#endif
