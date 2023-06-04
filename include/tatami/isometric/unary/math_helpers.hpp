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
 */
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
    template<typename Value_, typename Index_>
    void core (Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::abs(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry.
 */
template<typename T = double>
struct DelayedLogHelper {
    /**
     * Defaults to the natural log.
     */
    DelayedLogHelper() : log_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedLogHelper(double base) : log_base(std::log(base)) {}

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
    template<typename Value_, typename Index_>
    void core (Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::log(buffer[i]) / log_base;
        }
    }

    const double log_base;

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the square root of a matrix entry.
 */
template<typename T = double>
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
    template<typename Value_, typename Index_>
    void core (Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::sqrt(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry plus 1.
 */
template<typename T = double>
struct DelayedLog1pHelper {
    /**
     * Defaults to the natural log.
     */
    DelayedLog1pHelper() : log_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedLog1pHelper(double base) : log_base(std::log(base)) {}

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
    template<typename Value_, typename Index_>
    void core (Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::log1p(buffer[i]) / log_base;
        }
    }

    const double log_base;

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Round a matrix entry to the nearest integer.
 */
template<typename T = double>
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
    template<typename Value_, typename Index_>
    void core (Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::round(buffer[i]);
        }
    }

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent.
 */
template<typename T = double>
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
    template<typename Value_, typename Index_>
    void core (Index_ length, Value_* buffer) const {
    }

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::exp(buffer[i]);
        }
    }

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            if (buffer[i]) {
                buffer[i] = std::exp(buffer[i]);
            } else {
                buffer[i] = 1;
            }
        }
    }
    /**
     * @endcond
     */
};

}

#endif
