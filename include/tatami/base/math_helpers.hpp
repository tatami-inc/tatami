#ifndef TATAMI_MATH_HELPERS_H
#define TATAMI_MATH_HELPERS_H

/**
 * @file math_helpers.hpp
 *
 * Helper functions focusing on unary math operations,
 * to be used as the `OP` in the `DelayedIsometricOp` class.
 */

#include <cmath>

namespace tatami {

/**
 * @brief Take the absolute value of a matrix entry.
 */
template<typename T = double>
struct DelayedAbsHelper {
    /**
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value.
     *
     * @return Absolute value of `val`.
     */
    T operator()(size_t r, size_t c, T val) const {
        return std::abs(val);
    }

    /**
     * Sparsity is always preserved.
     */
    static const bool sparse = true;
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

    /**
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value.
     *
     * @return Logarithm of `val` with the specified base.
     */
    T operator()(size_t r, size_t c, T val) const {
        return std::log(val)/log_base;
    }

    /**
     * Sparsity is always discarded
     */
    static const bool sparse = false;
private:
    const double log_base;
};

/**
 * @brief Take the square root of a matrix entry.
 */
template<typename T = double>
struct DelayedSqrtHelper {
    /**
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value.
     *
     * @return Square root of `val`.
     */
    T operator()(size_t r, size_t c, T val) const {
        return std::sqrt(val);
    }

    /**
     * Sparsity is always preserved.
     */
    static const bool sparse = true;
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

    /**
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value.
     *
     * @return Logarithm of `val + 1` with the specified base.
     */
    T operator()(size_t r, size_t c, T val) const {
        return std::log1p(val)/log_base;
    }

    /**
     * Sparsity is always preserved.
     */
    static const bool sparse = true;
private:
    const double log_base;
};

/**
 * @brief Round a matrix entry to the nearest integer.
 */
template<typename T = double>
struct DelayedRoundHelper {
    /**
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value.
     *
     * @return `val` rounded to the closest integer value.
     */
    T operator()(size_t r, size_t c, T val) const {
        return std::round(val);
    }

    /**
     * Sparsity is always preserved.
     */
    static const bool sparse = true;
};

/**
 * @brief Use a matrix entry as an exponent.
 */
template<typename T = double>
struct DelayedExpHelper {
    /**
     * @param r Row index, ignored.
     * @param c Column index, ignored.
     * @param val Matrix value.
     *
     * @return `e` to the power of `val`.
     */
    T operator()(size_t r, size_t c, T val) const {
        return std::exp(val);
    }

    /**
     * Sparsity is always discarded.
     */
    static const bool sparse = false;
};

}

#endif
