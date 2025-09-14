#ifndef TATAMI_ISOMETRIC_UNARY_SUBSTITUTE_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_SUBSTITUTE_HELPERS_H

#include "../compare_utils.hpp"
#include "../../utils/Index_to_container.hpp"
#include "helper_interface.hpp"

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
template<CompareOperation op_, typename InputValue_, typename Scalar_, typename OutputValue_>
bool delayed_substitute_is_sparse(const Scalar_ compared, const OutputValue_ substitute) {
    return !delayed_compare<op_, InputValue_>(0, compared) || substitute == 0;
}

template<CompareOperation op_, typename InputValue_, typename Scalar_, typename OutputValue_>
OutputValue_ delayed_substitute_run(const InputValue_ val, const Scalar_ compared, const OutputValue_ substitute) {
    if (delayed_compare<op_>(val, compared)) {
        return substitute;
    } else {
        return val;
    }
}

template<CompareOperation op_, typename InputValue_, typename Index_, typename Scalar_, typename OutputValue_>
void delayed_substitute_run_simple(const InputValue_* input, const Index_ length, const Scalar_ compared, OutputValue_* const output, const OutputValue_ substitute) {
    if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
        input = output; // basically an assertion to the compiler to skip aliasing protection.
    }
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_substitute_run<op_>(input[i], compared, substitute);
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
 * @tparam Scalar_ Type of the scalar value.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
class DelayedUnaryIsometricSubstituteScalarHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param compared Scalar to be compared to the matrix values.
     * The matrix value is assumed to be on the left hand side of the comparison, while `compared` is on the right.
     * @param substitute Scalar to substitute into the matrix for every element where the comparison to `compared` is true.
     */
    DelayedUnaryIsometricSubstituteScalarHelper(const Scalar_ compared, const OutputValue_ substitute) :
        my_compared(compared),
        my_substitute(substitute)
    {
        my_sparse = delayed_substitute_is_sparse<op_, InputValue_>(my_compared, my_substitute);
    }

private:
    Scalar_ my_compared;
    OutputValue_ my_substitute;
    bool my_sparse;

public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        delayed_substitute_run_simple<op_>(input, length, my_compared, output, my_substitute);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        delayed_substitute_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_compared, output, my_substitute);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input_value, const Index_* const, OutputValue_* const output_value) const {
        delayed_substitute_run_simple<op_>(input_value, number, my_compared, output_value, my_substitute);
    }

    OutputValue_ fill(const bool, const Index_) const {
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
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteEqualScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "greater than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteGreaterThanScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "greater than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteGreaterThanScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "less than" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteLessThanScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "greater than or equal" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteGreaterThanOrEqualScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "less than or equal" substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteLessThanOrEqualScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar non-equality substitution helper.
 *
 * @tparam OutputValue_ Type of the result of the substitution.
 * @tparam InputValue_ Type of the matrix value used in the substitution.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubstituteNotEqualScalarHelper = DelayedUnaryIsometricSubstituteScalarHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteEqualScalar(Scalar_ compared, OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteGreaterThanScalar(Scalar_ compared, OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteGreaterThanScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteLessThanScalar(Scalar_ compared, OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteLessThanScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteGreaterThanOrEqualScalar(Scalar_ compared, OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteGreaterThanOrEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteLessThanOrEqualScalar(Scalar_ compared, OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteLessThanOrEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(compared, substitute);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteNotEqualScalar(Scalar_ compared, OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteNotEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(compared, substitute);
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
 * This should have a `[]` accessor method and a `size()` method.
 * @tparam SubstituteVector_ Type of the vector containing values to substitute in the output matrix.
 * This should have a `[]` accessor method and a `size()` method.
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
    DelayedUnaryIsometricSubstituteVectorHelper(ComparedVector_ compared, SubstituteVector_ substitute, const bool by_row) : 
        my_compared(std::move(compared)),
        my_substitute(std::move(substitute)),
        my_by_row(by_row) 
    {
        const auto ncompared = my_compared.size();
        if (!safe_non_negative_equal(ncompared, my_substitute.size())) {
            throw std::runtime_error("'compared' and 'substitute' should have the same length");
        }
        for (I<decltype(ncompared)> i = 0; i < ncompared; ++i) {
            if (!delayed_substitute_is_sparse<op_, InputValue_>(my_compared[i], my_substitute[i])) {
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
    std::optional<Index_> nrow() const {
        if (my_by_row) {
            return my_compared.size();
        } else {
            return std::nullopt;
        }
    }

    std::optional<Index_> ncol() const {
        if (my_by_row) {
            return std::nullopt;
        } else {
            return my_compared.size();
        }
    }

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
    void dense(const bool row, const Index_ idx, const Index_ start, const Index_ length, const InputValue_* input, OutputValue_* const output) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_>(input, length, my_compared[idx], output, my_substitute[idx]);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to skip aliasing protection.
            }
            for (Index_ i = 0; i < length; ++i) {
                const Index_ is = i + start;
                output[i] = delayed_substitute_run<op_>(input[i], my_compared[is], my_substitute[is]);
            }
        }
    }

    void dense(const bool row, const Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* const output) const {
        if (row == my_by_row) {
            delayed_substitute_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_compared[idx], output, my_substitute[idx]);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to skip aliasing protection.
            }
            const Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                const auto ii = indices[i];
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
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input_value = output_value; // basically an assertion to the compiler to skip aliasing protection.
            }
            for (Index_ i = 0; i < number; ++i) {
                const auto ii = indices[i];
                output_value[i] = delayed_substitute_run<op_>(input_value[i], my_compared[ii], my_substitute[ii]);
            }
        }
    }

    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            const auto sub = my_substitute[idx];
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
using DelayedUnaryIsometricSubstituteEqualVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

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
using DelayedUnaryIsometricSubstituteGreaterThanVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

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
using DelayedUnaryIsometricSubstituteGreaterThanVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

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
using DelayedUnaryIsometricSubstituteLessThanVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

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
using DelayedUnaryIsometricSubstituteGreaterThanOrEqualVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

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
using DelayedUnaryIsometricSubstituteLessThanOrEqualVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

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
using DelayedUnaryIsometricSubstituteNotEqualVectorHelper = DelayedUnaryIsometricSubstituteVectorHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteEqualVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubstituteEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteGreaterThanVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubstituteGreaterThanVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteLessThanVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubstituteLessThanVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteGreaterThanOrEqualVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubstituteGreaterThanOrEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricLessSubstituteThanOrEqualVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubstituteLessThanOrEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename ComparedVector_, typename SubstituteVector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteNotEqualVector(ComparedVector_ compared, SubstituteVector_ substitute, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubstituteNotEqualVectorHelper<OutputValue_, InputValue_, Index_, ComparedVector_, SubstituteVector_> >(std::move(compared), std::move(substitute), by_row);
}
/**
 * @endcond
 */

/**
 * @cond
 */
template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename OutputValue_>
bool delayed_special_substitute_is_sparse(const OutputValue_ substitute) {
    return !delayed_special_compare<op_, pass_, InputValue_>(0) || substitute == 0;
}

template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename OutputValue_>
OutputValue_ delayed_special_substitute_run(const InputValue_ val, const OutputValue_ substitute) {
    if (delayed_special_compare<op_, pass_>(val)) {
        return substitute;
    } else {
        return val;
    }
}

template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename Index_, typename OutputValue_>
void delayed_special_substitute_run_simple(const InputValue_* input, const Index_ length, const OutputValue_ substitute, OutputValue_* const output) {
    if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
        input = output; // basically an assertion to the compiler to skip aliasing protection.
    }
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_special_substitute_run<op_, pass_>(input[i], substitute);
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
    DelayedUnaryIsometricSpecialSubstituteHelper(const OutputValue_ substitute) : my_substitute(substitute) {
        my_sparse = delayed_special_substitute_is_sparse<op_, pass_, InputValue_>(my_substitute);
    }

private:
    OutputValue_ my_substitute;
    bool my_sparse;

public:
    std::optional<Index_> nrow() const {
        return std::nullopt;
    }

    std::optional<Index_> ncol() const {
        return std::nullopt;
    }

public:
    bool zero_depends_on_row() const {
        return false;
    }

    bool zero_depends_on_column() const {
        return false;
    }

    bool non_zero_depends_on_row() const {
        return false;
    }

    bool non_zero_depends_on_column() const {
        return false;
    }

public:
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        delayed_special_substitute_run_simple<op_, pass_>(input, length, my_substitute, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        delayed_special_substitute_run_simple<op_, pass_>(input, static_cast<Index_>(indices.size()), my_substitute, output);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input_value, const Index_* const, OutputValue_* const output_value) const {
        delayed_special_substitute_run_simple<op_, pass_>(input_value, number, my_substitute, output_value);
    }

    OutputValue_ fill(const bool, const Index_) const {
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
using DelayedUnaryIsometricSubstituteIsnanHelper = DelayedUnaryIsometricSpecialSubstituteHelper<SpecialCompareOperation::ISNAN, pass_, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the infinity substitution helper.
 *
 * @tparam pass_ Whether to substitute if the input value is infinite.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricSubstituteIsinfHelper = DelayedUnaryIsometricSpecialSubstituteHelper<SpecialCompareOperation::ISINF, pass_, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the finite value substitution helper.
 *
 * @tparam pass_ Whether to substitute if the input value is finite.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricSubstituteIsfiniteHelper = DelayedUnaryIsometricSpecialSubstituteHelper<SpecialCompareOperation::ISFINITE, pass_, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
// Back-compatibility only.
template<bool pass_ = true, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteIsnan(OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteIsnanHelper<pass_, OutputValue_, InputValue_, Index_> >(substitute);
}

template<bool pass_ = true, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteIsinf(OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteIsinfHelper<pass_, OutputValue_, InputValue_, Index_> >(substitute);
}

template<bool pass_ = true, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubstituteIsfinite(OutputValue_ substitute) {
    return std::make_shared<DelayedUnaryIsometricSubstituteIsfiniteHelper<pass_, OutputValue_, InputValue_, Index_> >(substitute);
}
/**
 * @endcond
 */

}

#endif
