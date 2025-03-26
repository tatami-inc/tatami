#ifndef TATAMI_ISOMETRIC_UNARY_SUBSTITUTE_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_SUBSTITUTE_HELPERS_H

#include "../compare_utils.hpp"
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
template<CompareOperation op_, typename InputValue_, typename OutputValue_>
bool delayed_substitute_is_sparse(InputValue_ compared, OutputValue_ substitute) {
    return !delayed_compare<op_, InputValue_>(0, compared) || substitute == 0;
}

template<CompareOperation op_, typename InputValue_, typename OutputValue_>
OutputValue_ delayed_substitute_run(InputValue_ val, InputValue_ compared, OutputValue_ substitute) {
    if (delayed_compare<op_, InputValue_>(val, compared)) {
        return substitute;
    } else {
        return val;
    }
}

template<CompareOperation op_, typename InputValue_, typename OutputValue_, typename Index_>
void delayed_substitute_run_simple(const InputValue_* buffer, Index_ length, InputValue_ compared, OutputValue_* output, OutputValue_ substitute) {
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_substitute_run<op_, Value_>(buffer[i], compared, substitute);
    }
}
/**
 * @endcond
 */

/**
 * @brief Helper for delayed scalar substitution.
 *
 * This class compares each element of a `Matrix` to a scalar;
 * when this comparison is true, it replaces the matrix element with another scalar value.
 *
 * @tparam op_ The comparison operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSubstituteScalarHelper final : public DelayeUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param compared Scalar to be compared to the matrix values.
     * The matrix value is assumed to be on the left hand side of the comparison, while `compared` is on the right.
     * @param substitute Scalar to substitute into the matrix for every element where the comparison to `compared` is true.
     */
    DelayedUnaryIsometricSubstituteScalar(InputValue_ compared, OutputValue_ substitute) :
        my_compared(compared),
        my_substitute(substitute)
    {
        my_sparse = delayed_substitute_is_sparse<op_>(my_compared, my_substitute);
    }

private:
    InputValue_ my_compared;
    outputValue_ my_substitute;
    bool my_sparse;

public:
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_substitute_run_simple<op_>(input, length, my_compared, output, my_substitute);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_substitute_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_compared, output, my_substitute);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_substitute_run_simple<op_, Value_>(input, number, my_compared, output_value, my_substitute);
    }

    OutputValue_ fill(bool, Index_) const {
        if (my_sparse) {
            return 0;
        } else {
            return my_substitute;
        }
    }
};

/**
 * Convenient alias for the scalar equality substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteEqualScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the scalar "greater than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteGreaterThanScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the scalar "greater than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteGreaterThanScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the scalar "less than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteLessThanScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the scalar "greater than or equal" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteGreaterThanOrEqualScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the scalar "less than or equal" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteLessThanOrEqualScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the scalar non-equality substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedBinaryIsometricSubstituteNotEqualScalarHelper = DelayedBinaryIsometricSubstituteScalarHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteEqualScalar(OutputValue_ compared, InputValue_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteEqualScalarHelper<OutputValue_, InputValue_, Index_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteGreaterThanScalar(OutputValue_ compared, InputValue_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteGreaterThanScalarHelper<OutputValue_, InputValue_, Index_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteLessThanScalar(OutputValue_ compared, InputValue_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteLessThanScalarHelper<OutputValue_, InputValue_, Index_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteGreaterThanOrEqualScalar(OutputValue_ compared, InputValue_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteGreaterThanOrEqualScalarHelper<OutputValue_, InputValue_, Index_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricLessSubstituteThanOrEqualScalar(OutputValue_ compared, InputValue_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteLessThanOrEqualScalarHelper<OutputValue_, InputValue_, Index_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteNotEqualScalar(OutputValue_ compared, InputValue_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteNotEqualScalarHelper<OutputValue_, InputValue_, Index_> >(compared, substitute);
}
/**
 * @endcond
 */

/**
 * @brief Delayed vector comparisons.
 *
 * This class compares each element of a `Matrix` to a row/column-specific value;
 * when this comparison is true, it replaces the matrix element with another row/column-specific value.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class, and only when `InputValue_ == OutputValue_`.
 *
 * @tparam op_ The comparison operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
class DelayedUnaryIsometricSubstituteVectorHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
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
    DelayedUnaryIsometricSubstituteVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) : 
        my_compared(std::move(compared)),
        my_substitute(std::move(substitute)),
        my_by_row(by_row) 
    {
        for (size_t i = 0, end = my_compared.size(); i < end; ++i) {
            if (!delayed_substitute_is_sparse<op_>(my_compared[i], my_substitute[i])) {
                 my_sparse = false;
                 break;
             }
        }
    }

private:
    ComparedVector_ my_compared;
    SubstituteVector_ my_substitute;
    bool my_by_row;
    bool my_sparse = true;

public:
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

public:
    void dense(bool row, Index_ idx, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_>(input, length, my_compared[idx], output, my_substitute[idx]);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                Index_ is = i + start;
                output[i] = delayed_substitute_run<op_>(input[i], my_compared[is], my_substitute[is]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const Value_*, Value_* output) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_compared[idx], output, my_substitute[idx]);
        } else {
            Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                auto ii = indices[i];
                output[i] = delayed_substitute_run<op_>(input[i], my_compared[ii], my_substitute[ii]);
            }
        }
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input_value, const Index_* indices, OutputValue_* output_value) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_>(input_value, number, my_compared[idx], output_value, my_substitute[idx]);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                auto ii = indices[i];
                output_value[i] = delayed_substitute_run<op_>(input_value[i], my_compared[ii], my_substitute[ii]);
            }
        }
    }

    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            auto sub = my_substitute[idx];
            if (!delayed_compare<op_, InputValue_>(0, my_compared[idx])) {
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
};

/**
 * Convenient alias for the scalar equality substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteEqualVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * Convenient alias for the scalar "greater than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteGreaterThanVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * Convenient alias for the scalar "greater than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteGreaterThanVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * Convenient alias for the scalar "less than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteLessThanVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * Convenient alias for the scalar "greater than or equal" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteGreaterThanOrEqualVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * Convenient alias for the scalar "less than or equal" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteLessThanOrEqualVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * Convenient alias for the scalar non-equality substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam ComparedVector_ Type of the vector containing values to compare to the input matrix.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename ComparedVector_, typename SubstituteVector_>
using DelayedBinaryIsometricSubstituteNotEqualVectorHelper = DelayedBinaryIsometricSubstituteVectorHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteEqualVector(ComparedVector_ compared, SubstituteVector_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteGreaterThanVector(ComparedVector_ compared, SubstituteVector_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteGreaterThanVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteLessThanVector(ComparedVector_ compared, SubstituteVector_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteLessThanVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteGreaterThanOrEqualVector(ComparedVector_ compared, SubstituteVector_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteGreaterThanOrEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricLessSubstituteThanOrEqualVector(ComparedVector_ compared, SubstituteVector_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteLessThanOrEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedBinaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedBinaryIsometricSubstituteNotEqualVector(ComparedVector_ compared, SubstituteVector_ substitute) {
    return std::make_shared<DelayedBinaryIsometricSubstituteNotEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute));
}
/**
 * @endcond
 */

/**
 * @cond
 */
template<SpecialCompareOperation op_, bool pass_, typename OutputValue_>
bool delayed_special_substitute_is_sparse(OutputValue_ substitute) {
    return !delayed_special_compare<op_, pass_, OutputValue_>(0) || substitute == 0;
}

template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename OutputValue_>
OutputValue_ delayed_special_substitute_run(InputValue_ val, Value_ substitute) {
    if (delayed_special_compare<op_, pass_, InputValue_>(val)) {
        return substitute;
    } else {
        return val;
    }
}

template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename Index_, typename OutputValue_>
void delayed_special_substitute_run_simple(const InputValue_* buffer, Index_ length, OutputValue_ substitute, OutputValue_* output) {
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_special_substitute_run<op_, pass_, Value_>(buffer[i], substitute);
    }
}
/**
 * @endcond
 */

/**
 * @brief Delayed special value substitution.
 *
 * This class checks whether each element of a `Matrix` is an IEEE special value, and if so, replaces it with another scalar value.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class, and only when `InputValue_ == OutputValue_`.
 *
 * @tparam op_ The special comparison operation.
 * @tparam pass_ Whether to perform the substitution if the special comparison is true.
 * Otherwise the substitution is only performed if the comparison is false.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<SpecialCompareOperation op_, bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSpecialSubstituteHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param substitute Scalar to substitute into the matrix for every element where the special comparison is true (if `pass_ = true`) or false (otherwise).
     */
    DelayedUnaryIsometricSpecialSubstitute(OutputValue_ substitute) : my_substitute(substitute) {
        my_sparse = delayed_special_substitute_is_sparse<op_, pass_, Value_>(my_substitute);
    }

private:
    Value_ my_substitute;
    bool my_sparse;

public:
    template<typename Index_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_special_substitute_run_simple<op_, pass_, Value_>(input, length, my_substitute, output);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_special_substitute_run_simple<op_, pass_, Value_>(input, static_cast<Index_>(indices.size()), my_substitute, output);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_special_substitute_run_simple<op_, pass_, Value_>(input_value, number, my_substitute, output_value);
    }

    Value_ fill(bool, Index_) const {
        if (my_sparse) {
            return 0;
        } else {
            return my_substitute;
        }
    }
};

/**
 * Convenient alias for the NaN substitution helper.
 *
 * @tparam pass_ Whether to substitute if the input value is NaN.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricSubstituteIsnan = DelayedUnaryIsometricSpecialSubstituteHelper<SpecialCompareOperation::ISNAN, pass_, OutputValue_, InputValue_, Index_>

/**
 * Convenient alias for the infinity substitution helper.
 *
 * @tparam pass_ Whether to substitute if the input value is infinite.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricSubstituteIsinf = DelayedUnaryIsometricSpecialSubstituteHelper<SpecialCompareOperation::ISINF, pass_, OutputValue_, InputValue_, Index_>

/**
 * Convenient alias for the finite value substitution helper.
 *
 * @tparam pass_ Whether to substitute if the input value is finite.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricSubstituteIsfinite = DelayedUnaryIsometricSpecialSubstituteHelper<SpecialCompareOperation::ISFINITE, pass_, OutputValue_, InputValue_, Index_>

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteIsnan(OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteIsnanHelper<OutputValue_, InputValue_, Index_> >(substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteIsinf(OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteIsinfHelper<OutputValue_, InputValue_, Index_> >(substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteIsfinite(OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteIsfiniteVectorHelper<OutputValue_, InputValue_, Index_> >(substitute);
}
/**
 * @endcond
 */

}

#endif
