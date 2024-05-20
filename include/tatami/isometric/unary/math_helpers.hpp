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
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core<Index_>(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the sign of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core<Index_>(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 * @tparam Base_ Numeric type for the log base.
 */
template<typename Value_ = double, typename Base_ = Value_>
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
    const Base_ my_base;

    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::log(buffer[i]) / my_base;
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core<Index_>(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        // Use the implementation-defined value.
        return std::log(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the square root of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the ceiling of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core<Index_>(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the floor of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the trunc of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of a matrix entry plus 1.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 * @tparam Base_ Numeric type for the log base.
 */
template<typename Value_ = double, typename Base_ = Value_>
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
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = std::log1p(buffer[i]) / my_base;
        }
    }

    const Base_ my_base;

public:
    /**
     * @cond
     */
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Round a matrix entry to the nearest integer.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Use a matrix entry as an exponent minus 1.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the arc cosine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
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
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::acosh(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the arc sine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic sine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the arc tangent of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the inverse hyperbolic tangent of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the cosine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic cosine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 1.0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the sine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic sine of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the tangent of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the hyperbolic tangent of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return 0;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the gamma of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::tgamma(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

/**
 * @brief Take the logarithm of the gamma of a matrix entry.
 *
 * This can be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * 
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
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
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        core(length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        core(indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        core(number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        // Use the implementation-defined special value.
        return std::lgamma(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

}

#endif
