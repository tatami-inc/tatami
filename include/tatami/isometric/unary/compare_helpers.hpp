#ifndef TATAMI_COMPARE_HELPERS_H
#define TATAMI_COMPARE_HELPERS_H

#include "../compare_utils.hpp"

/**
 * @file compare_helpers.hpp
 *
 * @brief Helper classes for delayed unary comparison operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 */

namespace tatami {

/**
 * @cond
 */
template<DelayedCompareOp op_, typename Scalar_, typename Value_, typename Index_>
void delayed_compare_run_simple(Scalar_ scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_compare_run<op_>(buffer[i], scalar);
    }
}

template<DelayedCompareOp op_, typename Scalar_>
bool delayed_compare_actual_sparse(Scalar_ scalar) {
    if constexpr(op_ == DelayedCompareOp::EQUAL) {
        return scalar != 0;
    } else if constexpr(op_ == DelayedCompareOp::GREATER_THAN) {
        return scalar >= 0;
    } else if constexpr(op_ == DelayedCompareOp::LESS_THAN) {
        return scalar <= 0;
    } else if constexpr(op_ == DelayedCompareOp::GREATER_THAN_OR_EQUAL) {
        return scalar > 0;
    } else if constexpr(op_ == DelayedCompareOp::LESS_THAN_OR_EQUAL) {
        return scalar < 0;
    } else { // NOT_EQUAL.
        return scalar == 0;
    }
}
/**
 * @endcond
 */

/**
 * @brief Delayed scalar comparison.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam Scalar_ Type of the scalar value.
 */
template<DelayedCompareOp op_, typename Scalar_>
struct DelayedCompareScalarHelper {
    /**
     * @param s Scalar value to be added.
     */
    DelayedCompareScalarHelper(Scalar_ s) : scalar(s) {}

private:
    const Scalar_ scalar;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = false;

    bool actual_sparse() const {
        return delayed_compare_actual_sparse<op_>(scalar);
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_compare_run_simple<op_>(scalar, length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_compare_run_simple<op_>(scalar, number, buffer);
    }

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_compare_run_simple<op_>(scalar, length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector comparisons.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam margin_ Matrix dimension along which the operation is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and that value is subtracted from all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam Vector_ Type of the vector.
 */
template<DelayedCompareOp op_, int margin_, typename Vector_>
struct DelayedCompareVectorHelper {
    /**
     * @param v Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedCompareVectorHelper(Vector_ v) : vec(std::move(v)) {
        for (auto x : vec) {
             if (!delayed_compare_actual_sparse<op_>(x)) {
                 still_sparse = false;
                 break;
             }
        }
    }

private:
    const Vector_ vec;
    bool still_sparse = true;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = (margin_ == 0);

    static constexpr bool needs_column = (margin_ == 1);

    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = false;

    bool actual_sparse() const {
        return still_sparse;
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool accrow_, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_ idx, ExtractType_ start, Index_ length, Value_* buffer) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_compare_run_simple<op_>(vec[idx], length, buffer);

        } else if constexpr(std::is_same<ExtractType_, Index_>::value) {
            for (Index_ i = 0; i < length; ++i) {
                delayed_compare_run<op_>(buffer[i], vec[i + start]);
            }

        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_compare_run<op_>(buffer[i], vec[start[i]]);
            }
        }
    }

    template<bool accrow_, typename Value_, typename Index_>
    void sparse(Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_compare_run_simple<op_>(vec[idx], number, buffer);

        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_compare_run<op_>(buffer[i], vec[indices[i]]);
            }
        }
    }

    template<bool accrow_, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_ idx, ExtractType_&& start, Index_ length, Value_* buffer) const {
        dense<accrow_>(idx, std::forward<ExtractType_>(start), length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be compared.
 * @return A helper class for a delayed equality comparison to a scalar.
 */
template<typename Scalar_>
DelayedCompareScalarHelper<DelayedCompareOp::EQUAL, Scalar_> make_DelayedEqualScalarHelper(Scalar_ s) {
    return DelayedCompareScalarHelper<DelayedCompareOp::EQUAL, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be compared.
 * @return A helper class for a delayed greater-than comparison to a scalar.
 */
template<typename Scalar_>
DelayedCompareScalarHelper<DelayedCompareOp::GREATER_THAN, Scalar_> make_DelayedGreaterThanScalarHelper(Scalar_ s) {
    return DelayedCompareScalarHelper<DelayedCompareOp::GREATER_THAN, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be compared.
 * @return A helper class for a delayed less-than comparison to a scalar.
 */
template<typename Scalar_>
DelayedCompareScalarHelper<DelayedCompareOp::LESS_THAN, Scalar_> make_DelayedLessThanScalarHelper(Scalar_ s) {
    return DelayedCompareScalarHelper<DelayedCompareOp::LESS_THAN, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be compared.
 * @return A helper class for a delayed greater-than-or-equal comparison to a scalar.
 */
template<typename Scalar_>
DelayedCompareScalarHelper<DelayedCompareOp::GREATER_THAN_OR_EQUAL, Scalar_> make_DelayedGreaterThanOrEqualScalarHelper(Scalar_ s) {
    return DelayedCompareScalarHelper<DelayedCompareOp::GREATER_THAN_OR_EQUAL, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be compared.
 * @return A helper class for a delayed less-than-or-equal comparison to a scalar.
 */
template<typename Scalar_>
DelayedCompareScalarHelper<DelayedCompareOp::LESS_THAN_OR_EQUAL, Scalar_> make_DelayedLessThanOrEqualScalarHelper(Scalar_ s) {
    return DelayedCompareScalarHelper<DelayedCompareOp::LESS_THAN_OR_EQUAL, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be compared.
 * @return A helper class for a delayed non-equality comparison to a scalar.
 */
template<typename Scalar_>
DelayedCompareScalarHelper<DelayedCompareOp::NOT_EQUAL, Scalar_> make_DelayedNotEqualScalarHelper(Scalar_ s) {
    return DelayedCompareScalarHelper<DelayedCompareOp::NOT_EQUAL, Scalar_>(std::move(s));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedCompareVectorHelper`.
 * @param v Vector of values to be compared.
 * @return A helper class for a delayed equality comparison to a vector.
 */
template<int margin_, typename Vector_>
DelayedCompareVectorHelper<DelayedCompareOp::EQUAL, margin_, Vector_> make_DelayedEqualVectorHelper(Vector_ v) {
    return DelayedCompareVectorHelper<DelayedCompareOp::EQUAL, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedCompareVectorHelper`.
 * @param v Vector of values to be compared.
 * @return A helper class for a delayed greater-than comparison to a vector.
 */
template<int margin_, typename Vector_>
DelayedCompareVectorHelper<DelayedCompareOp::GREATER_THAN, margin_, Vector_> make_DelayedGreaterThanVectorHelper(Vector_ v) {
    return DelayedCompareVectorHelper<DelayedCompareOp::GREATER_THAN, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedCompareVectorHelper`.
 * @param v Vector of values to be compared.
 * @return A helper class for a delayed less-than comparison to a vector.
 */
template<int margin_, typename Vector_>
DelayedCompareVectorHelper<DelayedCompareOp::LESS_THAN, margin_, Vector_> make_DelayedLessThanVectorHelper(Vector_ v) {
    return DelayedCompareVectorHelper<DelayedCompareOp::LESS_THAN, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedCompareVectorHelper`.
 * @param v Vector of values to be compared.
 * @return A helper class for a delayed greater-than-or-equal comparison to a vector.
 */
template<int margin_, typename Vector_>
DelayedCompareVectorHelper<DelayedCompareOp::GREATER_THAN_OR_EQUAL, margin_, Vector_> make_DelayedGreaterThanOrEqualVectorHelper(Vector_ v) {
    return DelayedCompareVectorHelper<DelayedCompareOp::GREATER_THAN_OR_EQUAL, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedCompareVectorHelper`.
 * @param v Vector of values to be compared.
 * @return A helper class for a delayed less-than-or-equal comparison to a vector.
 */
template<int margin_, typename Vector_>
DelayedCompareVectorHelper<DelayedCompareOp::LESS_THAN_OR_EQUAL, margin_, Vector_> make_DelayedLessThanOrEqualVectorHelper(Vector_ v) {
    return DelayedCompareVectorHelper<DelayedCompareOp::LESS_THAN_OR_EQUAL, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedCompareVectorHelper`.
 * @param v Vector of values to be compared.
 * @return A helper class for a delayed non-equality comparison to a vector.
 */
template<int margin_, typename Vector_>
DelayedCompareVectorHelper<DelayedCompareOp::NOT_EQUAL, margin_, Vector_> make_DelayedNotEqualVectorHelper(Vector_ v) {
    return DelayedCompareVectorHelper<DelayedCompareOp::NOT_EQUAL, margin_, Vector_>(std::move(v));
}

}

#endif
