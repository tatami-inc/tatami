#ifndef TATAMI_BOOLEAN_HELPERS_H
#define TATAMI_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include <vector>

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper classes for delayed unary boolean operations.
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
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam Value_ Type of the data value.
 */
template<DelayedBooleanOp op_, typename Value_ = double>
class DelayedBooleanScalarHelper {
public:
    /**
     * @param s Scalar value.
     */
    DelayedBooleanScalarHelper(bool scalar) : my_scalar(scalar) {
        my_sparse = delayed_boolean_actual_sparse<op_, Value_>(my_scalar);
    }

private:
    const bool my_scalar;
    bool my_sparse;

public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = false;

    static constexpr bool zero_depends_on_column = false;

    static constexpr bool non_zero_depends_on_row = false;

    static constexpr bool non_zero_depends_on_column = false;

    bool is_sparse() const {
        return my_sparse;
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<typename Index_> 
    void dense(bool, Index_, Index_, Index_ length, Value_* buffer) const {
        delayed_boolean_run_simple<op_>(my_scalar, length, buffer);
    }

    template<typename Index_> 
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        delayed_boolean_run_simple<op_>(my_scalar, indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_boolean_run_simple<op_>(my_scalar, number, buffer);
    }

    template<typename Index_>
    Value_ fill(Index_) const {
        Value_ output = 0;
        delayed_boolean_run<op_>(output, my_scalar);
        return output;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed boolean NOT operation.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
class DelayedBooleanNotHelper {
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
            buffer[i] = !static_cast<bool>(buffer[i]);
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
        return 1;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector boolean operations.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam margin_ Matrix dimension along which the operation is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and that value is subtracted from all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 */
template<DelayedBooleanOp op_, int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
class DelayedBooleanVectorHelper {
public:
    /**
     * @param vector Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `margin_ = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedBooleanVectorHelper(Vector_ vector) : my_vector(std::move(vector)) {
        for (auto x : my_vector) {
             if (!delayed_boolean_actual_sparse<op_, Value_>(x)) {
                 my_sparse = false;
                 break;
             }
        }
    }

private:
    const Vector_ my_vector;
    bool my_sparse = true;

public:
    /**
     * @cond
     */
    static constexpr bool zero_depends_on_row = (margin_ == 0);

    static constexpr bool zero_depends_on_column = (margin_ == 1);

    static constexpr bool non_zero_depends_on_row = (margin_ == 0);

    static constexpr bool non_zero_depends_on_column = (margin_ == 1);

    bool is_sparse() const {
        return my_sparse;
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<typename Index_>
    void dense(bool row, Index_ idx, Index_ start, Index_ length, Value_* buffer) const {
        if (row == (margin_ == 0)) {
            delayed_boolean_run_simple<op_>(my_vector[idx], length, buffer);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_boolean_run<op_>(buffer[i], my_vector[i + start]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, Value_* buffer) const {
        if (row == (margin_ == 0)) {
            delayed_boolean_run_simple<op_>(my_vector[idx], indices.size(), buffer);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                delayed_boolean_run<op_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if (row == (margin_ == 0)) {
            delayed_boolean_run_simple<op_>(my_vector[idx], number, buffer);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_boolean_run<op_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    Value_ fill(Index_ idx) const {
        Value_ output = 0;
        delayed_boolean_run<op_>(output, my_vector[idx]);
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
DelayedBooleanScalarHelper<DelayedBooleanOp::AND, Value_> make_DelayedBooleanAndScalarHelper(bool scalar) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::AND, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to use in the operation.
 * @return A helper class for a delayed OR operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::OR> make_DelayedBooleanOrScalarHelper(bool scalar) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::OR, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to be used in the operation.
 * @return A helper class for a delayed XOR operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::XOR> make_DelayedBooleanXorScalarHelper(bool scalar) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::XOR, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @param s Scalar value to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a scalar.
 */
template<typename Value_ = double>
DelayedBooleanScalarHelper<DelayedBooleanOp::EQUAL> make_DelayedBooleanEqualScalarHelper(bool scalar) {
    return DelayedBooleanScalarHelper<DelayedBooleanOp::EQUAL, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param vector Vector of values to be used in the operation.
 * @return A helper class for a delayed AND operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::AND, margin_, Value_, Vector_> make_DelayedBooleanAndVectorHelper(Vector_ vector) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::AND, margin_, Value_, Vector_>(std::move(vector));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param vector Vector of values to be used in the operation.
 * @return A helper class for a delayed OR operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::OR, margin_, Value_, Vector_> make_DelayedBooleanOrVectorHelper(Vector_ vector) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::OR, margin_, Value_, Vector_>(std::move(vector));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param vector Vector of values to be used in the operation.
 * @return A helper class for a delayed XOR operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::XOR, margin_, Value_, Vector_> make_DelayedBooleanXorVectorHelper(Vector_ vector) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::XOR, margin_, Value_, Vector_>(std::move(vector));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @tparam margin_ Matrix dimension along which the comparison is to occur, see `DelayedBooleanVectorHelper`.
 * @param vector Vector of values to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a vector.
 */
template<int margin_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedBooleanVectorHelper<DelayedBooleanOp::EQUAL, margin_, Value_, Vector_> make_DelayedBooleanEqualVectorHelper(Vector_ vector) {
    return DelayedBooleanVectorHelper<DelayedBooleanOp::EQUAL, margin_, Value_, Vector_>(std::move(vector));
}

}

#endif
