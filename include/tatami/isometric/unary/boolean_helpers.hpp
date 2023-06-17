#ifndef TATAMI_BOOLEAN_HELPERS_H
#define TATAMI_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include <vector>

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper classes for delayed unary boolean operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 */

namespace tatami {

/**
 * @cond
 */
template<DelayedBooleanOp op_, typename Value_, typename Index_>
void delayed_boolean_run_simple(bool scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_boolean_run<op_>(buffer[i], scalar);
    }
}

template<DelayedBooleanOp op_, typename Value_>
bool delayed_boolean_actual_sparse(bool scalar) {
    Value_ output = 0;
    delayed_boolean_run<op_>(output, scalar);
    return output == 0;
}
/**
 * @endcond
 */

/**
 * @brief Delayed scalar boolean operation.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam Value_ Type of the data value.
 */
template<DelayedBooleanOp op_, typename Value_ = double>
struct DelayedBooleanScalarHelper {
    /**
     * @param s Scalar value.
     */
    DelayedBooleanScalarHelper(bool s) : scalar(s) {
        still_sparse = delayed_boolean_actual_sparse<op_, Value_>(scalar);
    }

private:
    const bool scalar;
    bool still_sparse;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = (op_ == DelayedBooleanOp::AND);

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
    template<bool, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_boolean_run_simple<op_>(scalar, length, buffer);
    }

    template<bool, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_boolean_run_simple<op_>(scalar, number, buffer);
    }

    template<bool, typename Index_>
    Value_ zero(Index_) const {
        Value_ output = 0;
        delayed_boolean_run<op_>(output, scalar);
        return output;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed boolean NOT operation.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
struct DelayedBooleanNotHelper {
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool always_dense = true;

    static constexpr bool always_sparse = false;
    /**
     * @endcond
     */

private:
    template<typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = !static_cast<bool>(buffer[i]);
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
        return 1;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector boolean operations.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam margin_ Matrix dimension along which the operation is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and that value is subtracted from all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 */
template<DelayedBooleanOp op_, int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
struct DelayedBooleanVectorHelper {
    /**
     * @param v Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `margin_ = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedBooleanVectorHelper(Vector_ v) : vec(std::move(v)) {
        for (auto x : vec) {
             if (!delayed_boolean_actual_sparse<op_, Value_>(x)) {
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

    static constexpr bool always_sparse = (op_ == DelayedBooleanOp::AND);

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
    template<bool accrow_, typename Index_, typename ExtractType_>
    void dense(Index_ idx, ExtractType_ start, Index_ length, Value_* buffer) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_boolean_run_simple<op_>(vec[idx], length, buffer);

        } else if constexpr(std::is_same<ExtractType_, Index_>::value) {
            for (Index_ i = 0; i < length; ++i) {
                delayed_boolean_run<op_>(buffer[i], vec[i + start]);
            }

        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_boolean_run<op_>(buffer[i], vec[start[i]]);
            }
        }
    }

    template<bool accrow_, typename Index_>
    void sparse(Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_boolean_run_simple<op_>(vec[idx], number, buffer);

        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_boolean_run<op_>(buffer[i], vec[indices[i]]);
            }
        }
    }

    template<bool, typename Index_>
    Value_ zero(Index_ idx) const {
        Value_ output = 0;
        delayed_boolean_run<op_>(output, vec[idx]);
        return output;
    }
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Type of the data value.
 * @return A helper class for a delayed NOT operation.
 */
template<typename Value_ = double>
DelayedBooleanNotHelper<Value_> make_DelayedBooleanNotHelper() {
    return DelayedBooleanNotHelper<Value_>();
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to use in the operation.
 * @return A helper class for a delayed AND operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::AND, Value_> make_DelayedBooleanAndScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::AND, Value_>(s);
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to use in the operation.
 * @return A helper class for a delayed OR operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::OR> make_DelayedBooleanOrScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::OR, Value_>(s);
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to be used in the operation.
 * @return A helper class for a delayed XOR operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::XOR> make_DelayedBooleanXorScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::XOR, Value_>(s);
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::EQUAL> make_DelayedBooleanEqualScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::EQUAL, Value_>(s);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed AND operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::AND, margin_, Value_, Vector_> make_DelayedBooleanAndVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::AND, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed OR operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::OR, margin_, Value_, Vector_> make_DelayedBooleanOrVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::OR, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed XOR operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::XOR, margin_, Value_, Vector_> make_DelayedBooleanXorVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::XOR, margin_, Value_, Vector_>(std::move(v));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::EQUAL, margin_, Value_, Vector_> make_DelayedBooleanEqualVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::EQUAL, margin_, Value_, Vector_>(std::move(v));
}

}

#endif
