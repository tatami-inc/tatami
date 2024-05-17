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
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
class DelayedAbsHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedSignHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedLogHelper {
public:
    /**
     * Defaults to the natural log.
     */
    DelayedLogHelper() : my_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedLogHelper(Base_ base) : my_base(std::log(base)) {}

public:
    /**
     * @cond
     */
    static const bool zero_depends_on_row = false;

    static const bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
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
class DelayedSqrtHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedCeilingHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedFloorHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedTruncHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedLog1pHelper {
public:
    /**
     * Defaults to the natural log.
     */
    DelayedLog1pHelper() : my_base(1) {}

    /**
     * @param base Base of the logarithm.
     */
    DelayedLog1pHelper(Base_ base) : my_base(std::log(base)) {}

public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedRoundHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedExpHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
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
class DelayedExpm1Helper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedAcosHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        // Use the implementation-defined special value.
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
class DelayedAcoshHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        // Use the implementation-defined special value.
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
class DelayedAsinHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedAsinhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedAtanHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedAtanhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedCosHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
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
class DelayedCoshHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
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
class DelayedSinHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedSinhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedTanHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedTanhHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        return 0;
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
class DelayedGammaHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        // Use the implementation-defined special value.
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
class DelayedLgammaHelper {
public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

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
    Value_ fill(Index_) const {
        // Use the implementation-defined special value.
        return std::lgamma(static_cast<Value_>(0));
    }
    /**
     * @endcond
     */
};

}

#endif
