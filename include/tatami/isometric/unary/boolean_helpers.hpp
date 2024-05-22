#ifndef TATAMI_ISOMETRIC_UNARY_BOOLEAN_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include <vector>

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric boolean operations.
 */

namespace tatami {

/**
 * @cond
 */
template<BooleanOperation op_, typename Value_, typename Index_>
void delayed_boolean_run_simple(bool scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_boolean_run<op_>(buffer[i], scalar);
    }
}

template<BooleanOperation op_, typename Value_>
bool delayed_boolean_actual_sparse(bool scalar) {
    Value_ output = 0;
    delayed_boolean_run<op_>(output, scalar);
    return output == 0;
}
/**
 * @endcond
 */

/**
 * @brief Delayed unary isometric scalar boolean operation.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam Value_ Type of the data value.
 */
template<BooleanOperation op_, typename Value_ = double>
class DelayedUnaryIsometricBooleanScalar {
public:
    /**
     * @param scalar Scalar value.
     */
    DelayedUnaryIsometricBooleanScalar(bool scalar) : my_scalar(scalar) {
        my_sparse = delayed_boolean_actual_sparse<op_, Value_>(my_scalar);
    }

private:
    const bool my_scalar;
    bool my_sparse;

public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

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
    Value_ fill(bool, Index_) const {
        Value_ output = 0;
        delayed_boolean_run<op_>(output, my_scalar);
        return output;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed unary isometric boolean NOT operation.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam Value_ Type of the data value.
 */
template<typename Value_ = double>
class DelayedUnaryIsometricBooleanNot {
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
    Value_ fill(bool, Index_) const {
        return 1;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed unary isometric vector boolean operations.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 */
template<BooleanOperation op_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
class DelayedUnaryIsometricBooleanVector {
public:
    /**
     * @param vector Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * @param by_row Whether `vector` corresponds to the rows.
     * If true, each element of the vector is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vector is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricBooleanVector(Vector_ vector, bool by_row) : my_vector(std::move(vector)), my_by_row(by_row) {
        for (auto x : my_vector) {
             if (!delayed_boolean_actual_sparse<op_, Value_>(x)) {
                 my_sparse = false;
                 break;
             }
        }
    }

private:
    const Vector_ my_vector;
    bool my_by_row;
    bool my_sparse = true;

public:
    /**
     * @cond
     */
    static constexpr bool is_basic = false;

    bool zero_depends_on_row() const {
        return my_by_row;
    }

    bool zero_depends_on_column() const {
        return !my_by_row;
    }

    bool non_zero_depends_on_row() const {
        return my_by_row;
    }

    bool non_zero_depends_on_column() const {
        return !my_by_row;
    }

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
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(my_vector[idx], length, buffer);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_boolean_run<op_>(buffer[i], my_vector[i + start]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, Value_* buffer) const {
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(my_vector[idx], indices.size(), buffer);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                delayed_boolean_run<op_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(my_vector[idx], number, buffer);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_boolean_run<op_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    Value_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            Value_ output = 0;
            delayed_boolean_run<op_>(output, my_vector[idx]);
            return output;
        } else {
            // We should only get to this point if it's sparse, otherwise no
            // single fill value would work across the length of my_vector.
            return 0;
        }
    }
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Type of the data value.
 * @return A helper class for a delayed NOT operation,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricBooleanNot<Value_> make_DelayedUnaryIsometricBooleanNot() {
    return DelayedUnaryIsometricBooleanNot<Value_>();
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to use in the operation.
 * @return A helper class for a delayed AND operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::AND, Value_> make_DelayedUnaryIsometricBooleanAndScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::AND, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to use in the operation.
 * @return A helper class for a delayed OR operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::OR> make_DelayedUnaryIsometricBooleanOrScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::OR, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be used in the operation.
 * @return A helper class for a delayed XOR operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::XOR> make_DelayedUnaryIsometricBooleanXorScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::XOR, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::EQUAL> make_DelayedUnaryIsometricBooleanEqualScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::EQUAL, Value_>(scalar);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed AND operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricBooleanVector<BooleanOperation::AND, Value_, Vector_> make_DelayedUnaryIsometricBooleanAndVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::AND, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed OR operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricBooleanVector<BooleanOperation::OR, Value_, Vector_> make_DelayedUnaryIsometricBooleanOrVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::OR, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed XOR operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricBooleanVector<BooleanOperation::XOR, Value_, Vector_> make_DelayedUnaryIsometricBooleanXorVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::XOR, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed boolean equality operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricBooleanVector<BooleanOperation::EQUAL, Value_, Vector_> make_DelayedUnaryIsometricBooleanEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::EQUAL, Value_, Vector_>(std::move(vector), by_row);
}

}

#endif
