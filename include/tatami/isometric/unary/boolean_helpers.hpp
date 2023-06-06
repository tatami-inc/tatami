#ifndef TATAMI_BOOLEAN_HELPERS_H
#define TATAMI_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"

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

template<DelayedBooleanOp op_>
bool delayed_boolean_actual_sparse(bool scalar) {
    if constexpr(op_ ==DelayedBooleanOp::OR || op_ == DelayedBooleanOp::XOR) {
        // for both or and xor, if the scalar is true, then
        // the result will not be sparse, as zeros + scalar will be true.
        return !scalar;
    } else { // i.e., EQUAL.
        return scalar;
    }
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
 */
template<DelayedBooleanOp op_>
struct DelayedBooleanScalarHelper {
    /**
     * @param s Scalar value.
     */
    DelayedBooleanScalarHelper(bool s) : scalar(s) {}

private:
    const bool scalar;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool always_dense = false;

    static constexpr bool always_sparse = (op_ == DelayedBooleanOp::AND);

    bool actual_sparse() const {
        return delayed_boolean_actual_sparse<op_>(scalar);
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
        delayed_boolean_run_simple<op_>(scalar, length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_boolean_run_simple<op_>(scalar, number, buffer);
    }

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_boolean_run_simple<op_>(scalar, length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed boolean NOT operation.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 */
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
    template<typename Value_, typename Index_>
    void core(Index_ length, Value_* buffer) const {
        for (Index_ i = 0; i < length; ++i) {
            buffer[i] = !static_cast<bool>(buffer[i]);
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

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        core(length, buffer);
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
 * @tparam Vector_ Type of the vector.
 */
template<DelayedBooleanOp op_, int margin_, typename Vector_>
struct DelayedBooleanVectorHelper {
    /**
     * @param v Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `margin_ = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedBooleanVectorHelper(Vector_ v) : vec(std::move(v)) {
        for (auto x : vec) {
             if (!delayed_boolean_actual_sparse<op_>(x)) {
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
    template<bool accrow_, typename Value_, typename Index_, typename ExtractType_>
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

    template<bool accrow_, typename Value_, typename Index_>
    void sparse(Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_boolean_run_simple<op_>(vec[idx], number, buffer);

        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_boolean_run<op_>(buffer[i], vec[indices[i]]);
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
 * @return A helper class for a delayed NOT operation.
 */
inline DelayedBooleanNotHelper make_DelayedBooleanNotHelper() {
    return DelayedBooleanNotHelper();
}

/**
 * @param s Scalar value to use in the operation.
 * @return A helper class for a delayed AND operation with a scalar.
 */
inline DelayedBooleanScalarHelper<DelayedBooleanOp::AND> make_DelayedBooleanAndScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::AND>(s);
}

/**
 * @param s Scalar value to use in the operation.
 * @return A helper class for a delayed OR operation with a scalar.
 */
inline DelayedBooleanScalarHelper<DelayedBooleanOp::OR> make_DelayedBooleanOrScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::OR>(s);
}

/**
 * @param s Scalar value to be used in the operation.
 * @return A helper class for a delayed XOR operation with a scalar.
 */
inline DelayedBooleanScalarHelper<DelayedBooleanOp::XOR> make_DelayedBooleanXorScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::XOR>(s);
}

/**
 * @param s Scalar value to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a scalar.
 */
inline DelayedBooleanScalarHelper<DelayedBooleanOp::EQUAL> make_DelayedBooleanEqualScalarHelper(bool s) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::EQUAL>(s);
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed AND operation with a vector.
 */
template<int margin_, typename Vector_>
DelayedBooleanVectorHelper<DelayedBooleanOp::AND, margin_, Vector_> make_DelayedBooleanAndVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::AND, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed OR operation with a vector.
 */
template<int margin_, typename Vector_>
DelayedBooleanVectorHelper<DelayedBooleanOp::OR, margin_, Vector_> make_DelayedBooleanOrVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::OR, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed XOR operation with a vector.
 */
template<int margin_, typename Vector_>
DelayedBooleanVectorHelper<DelayedBooleanOp::XOR, margin_, Vector_> make_DelayedBooleanXorVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::XOR, margin_, Vector_>(std::move(v));
}

/**
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param v Vector of values to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a vector.
 */
template<int margin_, typename Vector_>
DelayedBooleanVectorHelper<DelayedBooleanOp::EQUAL, margin_, Vector_> make_DelayedBooleanEqualVectorHelper(Vector_ v) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::EQUAL, margin_, Vector_>(std::move(v));
}

}

#endif
