#ifndef TATAMI_ISOMETRIC_UNARY_COMPARE_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_COMPARE_HELPERS_H

#include "../compare_utils.hpp"
#include <vector>

/**
 * @file compare_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric comparison operations.
 */

namespace tatami {

/**
 * @cond
 */
template<CompareOperation op_, typename Scalar_, typename Value_, typename Index_>
void delayed_compare_run_simple(Scalar_ scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_compare_run<op_>(buffer[i], scalar);
    }
}

template<CompareOperation op_, typename Value_, typename Scalar_>
bool delayed_compare_actual_sparse(Scalar_ scalar) {
    Value_ output = 0;
    delayed_compare_run<op_>(output, scalar);
    return output == 0;
}
/**
 * @endcond
 */

/**
 * @brief Delayed scalar comparison.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar value.
 */
template<CompareOperation op_, typename Value_ = double, typename Scalar_ = Value_>
class DelayedUnaryIsometricCompareScalar {
public:
    /**
     * @param scalar Scalar value to be added.
     */
    DelayedUnaryIsometricCompareScalar(Scalar_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_compare_actual_sparse<op_, Value_>(my_scalar);
    }

private:
    const Scalar_ my_scalar;
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
        delayed_compare_run_simple<op_>(my_scalar, length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        delayed_compare_run_simple<op_>(my_scalar, indices.size(), buffer);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_compare_run_simple<op_>(my_scalar, number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        Value_ output = 0;
        delayed_compare_run<op_>(output, my_scalar);
        return output;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector comparisons.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 */
template<CompareOperation op_, typename Value_ = double, typename Vector_ = std::vector<Value_> >
class DelayedUnaryIsometricCompareVector {
public:
    /**
     * @param vector Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * @param by_row Whether `vector` corresponds to the rows.
     * If true, each element of the vector is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vector is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricCompareVector(Vector_ vector, bool by_row) : my_vector(std::move(vector)), my_by_row(by_row) {
        for (auto x : my_vector) {
             if (!delayed_compare_actual_sparse<op_, Value_>(x)) {
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
            delayed_compare_run_simple<op_>(my_vector[idx], length, buffer);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_compare_run<op_>(buffer[i], my_vector[i + start]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, Value_* buffer) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_>(my_vector[idx], indices.size(), buffer);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                delayed_compare_run<op_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_>(my_vector[idx], number, buffer);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_compare_run<op_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    Value_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            Value_ output = 0;
            delayed_compare_run<op_>(output, my_vector[idx]);
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
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricCompareScalar<CompareOperation::EQUAL, Value_, Scalar_> make_DelayedUnaryIsometricEqualScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::EQUAL, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed greater-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN, Value_, Scalar_> make_DelayedUnaryIsometricGreaterThanScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed less-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN, Value_, Scalar_> make_DelayedUnaryIsometricLessThanScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed greater-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Scalar_> make_DelayedUnaryIsometricGreaterThanOrEqualScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed less-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Scalar_> make_DelayedUnaryIsometricLessThanOrEqualScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed non-equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricCompareScalar<CompareOperation::NOT_EQUAL, Value_, Scalar_> make_DelayedUnaryIsometricNotEqualScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::NOT_EQUAL, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed equality comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricCompareVector<CompareOperation::EQUAL, Value_, Vector_> make_DelayedUnaryIsometricEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::EQUAL, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed greater-than comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN, Value_, Vector_> make_DelayedUnaryIsometricGreaterThanVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed less-than comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN, Value_, Vector_> make_DelayedUnaryIsometricLessThanVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed greater-than-or-equal comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Vector_> make_DelayedUnaryIsometricGreaterThanOrEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed less-than-or-equal comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Vector_> make_DelayedUnaryIsometricLessThanOrEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed non-equality comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricCompareVector<CompareOperation::NOT_EQUAL, Value_, Vector_> make_DelayedUnaryIsometricNotEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::NOT_EQUAL, Value_, Vector_>(std::move(vector), by_row);
}

}

#endif
