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
template<CompareOperation op_, typename Value_, typename Index_>
void delayed_compare_run_simple(Value_* buffer, Index_ length, Value_ scalar) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_compare_run<op_, Value_>(buffer[i], scalar);
    }
}

template<CompareOperation op_, typename Value_>
bool delayed_compare_actual_sparse(Value_ scalar) {
    return !delayed_compare<op_, Value_>(0, scalar);
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
 */
template<CompareOperation op_, typename Value_ = double>
class DelayedUnaryIsometricCompareScalar {
public:
    /**
     * @param scalar Scalar to be compared to the matrix values.
     * The matrix value is assumed to be on the left hand side of the comparison, while `scalar` is on the right.
     */
    DelayedUnaryIsometricCompareScalar(Value_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_compare_actual_sparse<op_, Value_>(my_scalar);
    }

private:
    const Value_ my_scalar;
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
        delayed_compare_run_simple<op_, Value_>(buffer, length, my_scalar);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        delayed_compare_run_simple<op_, Value_>(buffer, indices.size(), my_scalar);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_compare_run_simple<op_, Value_>(buffer, number, my_scalar);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        Value_ output = 0;
        delayed_compare_run<op_, Value_>(output, my_scalar);
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
     * @param vector Vector to use in the comparison with the matrix values.
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * The matrix value from each row/column is assumed to be on the left hand side of the comparison, while the corresponding value of `vector` is on the right.
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
            delayed_compare_run_simple<op_, Value_>(buffer, length, my_vector[idx]);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_compare_run<op_, Value_>(buffer[i], my_vector[i + start]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, Value_* buffer) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, Value_>(buffer, indices.size(), my_vector[idx]);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                delayed_compare_run<op_, Value_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, Value_>(buffer, number, my_vector[idx]);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_compare_run<op_, Value_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    Value_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            Value_ output = 0;
            delayed_compare_run<op_, Value_>(output, my_vector[idx]);
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
 * @param scalar Value to be compared.
 * @return A helper class for a delayed equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::EQUAL, Value_> make_DelayedUnaryIsometricEqualScalar(Value_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::EQUAL, Value_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed greater-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN, Value_> make_DelayedUnaryIsometricGreaterThanScalar(Value_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN, Value_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed less-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN, Value_> make_DelayedUnaryIsometricLessThanScalar(Value_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN, Value_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed greater-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN_OR_EQUAL, Value_> make_DelayedUnaryIsometricGreaterThanOrEqualScalar(Value_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN_OR_EQUAL, Value_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed less-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN_OR_EQUAL, Value_> make_DelayedUnaryIsometricLessThanOrEqualScalar(Value_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN_OR_EQUAL, Value_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed non-equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::NOT_EQUAL, Value_> make_DelayedUnaryIsometricNotEqualScalar(Value_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::NOT_EQUAL, Value_>(std::move(scalar));
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

/**
 * @cond
 */
template<SpecialCompareOperation op_, bool pass_, typename Value_, typename Index_>
void delayed_special_compare_run_simple(Value_* buffer, Index_ length) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_special_compare_run<op_, pass_, Value_>(buffer[i]);
    }
}

template<SpecialCompareOperation op_, bool pass_, typename Value_>
bool delayed_special_compare_actual_sparse() {
    return !delayed_special_compare<op_, pass_, Value_>(0);
}
/**
 * @endcond
 */

/**
 * @brief Delayed special value comparison.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The special comparison operation.
 * @tparam pass_ Whether to return true if the special comparison is true.
 * @tparam Value_ Floating-point type of the data value.
 */
template<SpecialCompareOperation op_, bool pass_, typename Value_ = double>
class DelayedUnaryIsometricSpecialCompare {
public:
    DelayedUnaryIsometricSpecialCompare() {
        my_sparse = !delayed_special_compare<op_, pass_, Value_>(0);
    }

private:
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
        delayed_special_compare_run_simple<op_, pass_>(buffer, length);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        delayed_special_compare_run_simple<op_, pass_>(buffer, static_cast<Index_>(indices.size()));
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_special_compare_run_simple<op_, pass_>(buffer, number);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return static_cast<Value_>(!my_sparse);
    }
    /**
     * @endcond
     */
};

/**
 * @tparam pass_ Whether to return truthy if the matrix value is NaN.
 * @tparam Value_ Type of the data value.
 * @return A helper class for a delayed NaN check,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool pass_ = true, typename Value_ = double>
DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISNAN, pass_, Value_> make_DelayedUnaryIsometricIsnan() {
    return DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISNAN, pass_, Value_>();
}

/**
 * @tparam pass_ Whether to return truthy if the matrix value is infinite.
 * @tparam Value_ Type of the data value.
 * @return A helper class for a delayed check for infinity (positive or negative),
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool pass_ = true, typename Value_ = double>
DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISINF, pass_, Value_> make_DelayedUnaryIsometricIsinf() {
    return DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISINF, pass_, Value_>();
}

/**
 * @tparam pass_ Whether to return truthy if the matrix value is finite.
 * @tparam Value_ Type of the data value.
 * @return A helper class for a delayed check for finite values,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool pass_ = true, typename Value_ = double>
DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISFINITE, pass_, Value_> make_DelayedUnaryIsometricIsfinite() {
    return DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISFINITE, pass_, Value_>();
}

}

#endif
