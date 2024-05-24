#ifndef TATAMI_ISOMETRIC_UNARY_SUBSTITUTE_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_SUBSTITUTE_HELPERS_H

#include "../arithmetic_utils.hpp"
#include <vector>
#include <limits>

/**
 * @file substitute_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric substitution.
 */

namespace tatami {

/**
 * @cond
 */
template<CompareOperation op_, typename Value_, typename Scalar_>
bool delayed_substitute_is_sparse(Value_ compared, Scalar_ substitute) {
    return !delayed_compare<op_, Value_>(0, compared) || substitute == 0;
}

template<CompareOperation op_, typename Value_, typename Scalar_>
void delayed_substitute_run(Value_& val, Value_ compared, Scalar_ substitute) {
    if (delayed_compare<op_, Value_>(val, compared)) {
        val = substitute;
    }
}

template<CompareOperation op_, typename Value_, typename Index_, typename Scalar_>
void delayed_substitute_run_simple(Value_* buffer, Index_ length, Value_ compared, Scalar_ substitute) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_substitute_run<op_, Value_>(buffer[i], compared, substitute);
    }
}
/**
 * @endcond
 */


/**
 * @brief Delayed scalar substitution.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar value.
 */
template<CompareOperation op_, typename Value_ = double, typename Scalar_ = Value_>
class DelayedUnaryIsometricSubstituteScalar {
public:
    /**
     * @param compared Scalar value to be compared to the matrix values.
     * The matrix value is assumed to be on the left hand side of the comparison, while `compared` is on the right.
     * @param substitue Scalar value to substitute into the matrix for every element where the comparison to `compared` is true.
     */
    DelayedUnaryIsometricSubstituteScalar(Scalar_ compared, Scalar_ substitute) : my_compared(compared), my_substitute(substitute) {
        my_sparse = delayed_substitute_is_sparse<op_, Value_>(my_compared, my_substitute);
    }

private:
    Scalar_ my_compared, my_substitute;
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
        delayed_substitute_run_simple<op_, Value_>(buffer, length, my_compared, my_substitute);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        delayed_substitute_run_simple<op_, Value_>(buffer, static_cast<Index_>(indices.size()), my_compared, my_substitute);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_substitute_run_simple<op_, Value_>(buffer, number, my_compared, my_substitute);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        if (my_sparse) {
            return 0;
        } else {
            return my_substitute;
        }
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
class DelayedUnaryIsometricSubstituteVector {
public:
    /**
     * @param compared Vector to use in the comparison with the matrix values.
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * The matrix value from each row/column is assumed to be on the left hand side of the comparison, while the corresponding value of `compared` is on the right.
     * @param substitute Vector containing values to be substituted into the matrix at every element where the comparison to the corresponding element of `compared` is true.
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * @param by_row Whether `compared` and `substitute` corresponds to the rows.
     * If true, each element of the vectors is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vectors is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricSubstituteVector(Vector_ compared, Vector_ substitute, bool by_row) : 
        my_compared(std::move(compared)), my_substitute(std::move(substitute)), my_by_row(by_row) 
    {
        for (size_t i = 0, end = my_compared.size(); i < end; ++i) {
            if (!delayed_substitute_is_sparse<op_, Value_>(my_compared[i], my_substitute[i])) {
                 my_sparse = false;
                 break;
             }
        }
    }

private:
    Vector_ my_compared, my_substitute;
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
            delayed_substitute_run_simple<op_, Value_>(buffer, length, my_compared[idx], my_substitute[idx]);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                Index_ is = i + start;
                delayed_substitute_run<op_, Value_>(buffer[i], my_compared[is], my_substitute[is]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, Value_* buffer) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_, Value_>(buffer, static_cast<Index_>(indices.size()), my_compared[idx], my_substitute[idx]);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                auto ii = indices[i];
                delayed_substitute_run<op_, Value_>(buffer[i], my_compared[ii], my_substitute[ii]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_, Value_>(buffer, number, my_compared[idx], my_substitute[idx]);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                auto ii = indices[i];
                delayed_substitute_run<op_, Value_>(buffer[i], my_compared[ii], my_substitute[ii]);
            }
        }
    }

    template<typename Index_>
    Value_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            auto sub = my_substitute[idx];
            if (!delayed_compare<op_, Value_>(0, my_compared[idx])) {
                return 0;
            } else {
                return sub;
            }
        } else {
            // We should only get to this point if it's sparse, otherwise no
            // single fill value would work across the length of my_compared.
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
 * @param compared Scalar to be compared to the matrix values.
 * @param substitute Scalar value to substitute into the matrix when the comparison is true.
 * @return A helper class for a delayed equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricSubstituteScalar<CompareOperation::EQUAL, Value_, Scalar_> 
    make_DelayedUnaryIsometricSubstituteEqualScalar(Scalar_ compared, Scalar_ substitute) 
{
    return DelayedUnaryIsometricSubstituteScalar<CompareOperation::EQUAL, Value_, Scalar_>(std::move(compared), std::move(substitute));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param compared Scalar to be compared to the matrix values.
 * @param substitute Scalar value to substitute into the matrix when the comparison is true.
 * @return A helper class for a delayed greater-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricSubstituteScalar<CompareOperation::GREATER_THAN, Value_, Scalar_> 
    make_DelayedUnaryIsometricSubstituteGreaterThanScalar(Scalar_ compared, Scalar_ substitute)
{
    return DelayedUnaryIsometricSubstituteScalar<CompareOperation::GREATER_THAN, Value_, Scalar_>(std::move(compared), std::move(substitute));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param compared Scalar to be compared to the matrix values.
 * @param substitute Scalar value to substitute into the matrix when the comparison is true.
 * @return A helper class for a delayed less-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricSubstituteScalar<CompareOperation::LESS_THAN, Value_, Scalar_>
    make_DelayedUnaryIsometricSubstituteLessThanScalar(Scalar_ compared, Scalar_ substitute)
{
    return DelayedUnaryIsometricSubstituteScalar<CompareOperation::LESS_THAN, Value_, Scalar_>(std::move(compared), std::move(substitute));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param compared Scalar to be compared to the matrix values.
 * @param substitute Scalar value to substitute into the matrix when the comparison is true.
 * @return A helper class for a delayed greater-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricSubstituteScalar<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Scalar_> 
    make_DelayedUnaryIsometricSubstituteGreaterThanOrEqualScalar(Scalar_ compared, Scalar_ substitute)
{
    return DelayedUnaryIsometricSubstituteScalar<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Scalar_>(std::move(compared), std::move(substitute));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param compared Scalar to be compared to the matrix values.
 * @param substitute Scalar to substitute into the matrix when the comparison is true.
 * @return A helper class for a delayed less-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricSubstituteScalar<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Scalar_>
    make_DelayedUnaryIsometricSubstituteLessThanOrEqualScalar(Scalar_ compared, Scalar_ substitute)
{
    return DelayedUnaryIsometricSubstituteScalar<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Scalar_>(std::move(compared), std::move(substitute));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param compared Scalar to be compared to the matrix values.
 * @param substitute Scalar value to substitute into the matrix when the comparison is true.
 * @return A helper class for a delayed non-equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricSubstituteScalar<CompareOperation::NOT_EQUAL, Value_, Scalar_>
    make_DelayedUnaryIsometricSubstituteNotEqualScalar(Scalar_ compared, Scalar_ substitute) 
{
    return DelayedUnaryIsometricSubstituteScalar<CompareOperation::NOT_EQUAL, Value_, Scalar_>(std::move(compared), std::move(substitute));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param compared Vector to be compared to the matrix values.
 * @param substitute Vector containing values to substitute into the matrix when the comparison is true.
 * @param by_row Whether each element of `compared` and `substitute` corresponds to a row, see `DelayedUnaryIsometricSubstituteVector`.
 * @return A helper class for a delayed equality comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricSubstituteVector<CompareOperation::EQUAL, Value_, Vector_> 
    make_DelayedUnaryIsometricSubstituteEqualVector(Vector_ compared, Vector_ substitute, bool by_row) 
{
    return DelayedUnaryIsometricSubstituteVector<CompareOperation::EQUAL, Value_, Vector_>(std::move(compared), std::move(substitute), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param compared Vector to be compared to the matrix values.
 * @param substitute Vector containing values to substitute into the matrix when the comparison is true.
 * @param by_row Whether each element of `compared` and `substitute` corresponds to a row, see `DelayedUnaryIsometricSubstituteVector`.
 * @return A helper class for a delayed greater-than comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricSubstituteVector<CompareOperation::GREATER_THAN, Value_, Vector_> 
    make_DelayedUnaryIsometricSubstituteGreaterThanVector(Vector_ compared, Vector_ substitute, bool by_row)
{
    return DelayedUnaryIsometricSubstituteVector<CompareOperation::GREATER_THAN, Value_, Vector_>(std::move(compared), std::move(substitute), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param compared Vector to be compared to the matrix values.
 * @param substitute Vector containing values to substitute into the matrix when the comparison is true.
 * @param by_row Whether each element of `compared` and `substitute` corresponds to a row, see `DelayedUnaryIsometricSubstituteVector`.
 * @return A helper class for a delayed less-than comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricSubstituteVector<CompareOperation::LESS_THAN, Value_, Vector_> 
    make_DelayedUnaryIsometricSubstituteLessThanVector(Vector_ compared, Vector_ substitute, bool by_row) 
{
    return DelayedUnaryIsometricSubstituteVector<CompareOperation::LESS_THAN, Value_, Vector_>(std::move(compared), std::move(substitute), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param compared Vector to be compared to the matrix values.
 * @param substitute Vector containing values to substitute into the matrix when the comparison is true.
 * @param by_row Whether each element of `compared` and `substitute` corresponds to a row, see `DelayedUnaryIsometricSubstituteVector`.
 * @return A helper class for a delayed greater-than-or-equal comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricSubstituteVector<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Vector_> 
    make_DelayedUnaryIsometricSubstituteGreaterThanOrEqualVector(Vector_ compared, Vector_ substitute, bool by_row) 
{
    return DelayedUnaryIsometricSubstituteVector<CompareOperation::GREATER_THAN_OR_EQUAL, Value_, Vector_>(std::move(compared), std::move(substitute), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param compared Vector to be compared to the matrix values.
 * @param substitute Vector containing values to substitute into the matrix when the comparison is true.
 * @param by_row Whether each element of `compared` and `substitute` corresponds to a row, see `DelayedUnaryIsometricSubstituteVector`.
 * @return A helper class for a delayed less-than-or-equal comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricSubstituteVector<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Vector_> 
    make_DelayedUnaryIsometricSubstituteLessThanOrEqualVector(Vector_ compared, Vector_ substitute, bool by_row) 
{
    return DelayedUnaryIsometricSubstituteVector<CompareOperation::LESS_THAN_OR_EQUAL, Value_, Vector_>(std::move(compared), std::move(substitute), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 * @param compared Vector to be compared to the matrix values.
 * @param substitute Vector containing values to substitute into the matrix when the comparison is true.
 * @param by_row Whether each element of `compared` and `substitute` corresponds to a row, see `DelayedUnaryIsometricSubstituteVector`.
 * @return A helper class for a delayed non-equality comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<Value_> >
DelayedUnaryIsometricSubstituteVector<CompareOperation::NOT_EQUAL, Value_, Vector_> 
    make_DelayedUnaryIsometricSubstituteNotEqualVector(Vector_ compared, Vector_ substitute, bool by_row) 
{
    return DelayedUnaryIsometricSubstituteVector<CompareOperation::NOT_EQUAL, Value_, Vector_>(std::move(compared), std::move(substitute), by_row);
}

}

#endif
