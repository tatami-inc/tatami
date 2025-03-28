#ifndef TATAMI_ISOMETRIC_UNARY_COMPARE_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_COMPARE_HELPERS_H

#include "../compare_utils.hpp"
#include "helper_interface.hpp"
#include <vector>
#include <type_traits>

/**
 * @file compare_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric comparison operations.
 */

namespace tatami {

/**
 * @cond
 */
template<CompareOperation op_, typename InputValue_, typename Index_, typename Scalar_, typename OutputValue_>
void delayed_compare_run_simple(const InputValue_* input, Index_ length, Scalar_ scalar, OutputValue_* output) {
    if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
        input = output; // basically an assertion to the compiler to enable optimizations.
    }
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_compare<op_>(input[i], scalar);
    }
}

template<CompareOperation op_, typename InputValue_, typename Scalar_>
bool delayed_compare_actual_sparse(Scalar_ scalar) {
    return !delayed_compare<op_, InputValue_>(0, scalar);
}
/**
 * @endcond
 */

/**
 * @brief Helper for delayed scalar comparisons.
 *
 * This class compares each element of a `Matrix` to a scalar.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
class DelayedUnaryIsometricCompareScalarHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param scalar Scalar to be compared to the matrix values.
     * The matrix value is assumed to be on the left hand side of the comparison, while `scalar` is on the right.
     */
    DelayedUnaryIsometricCompareScalarHelper(Scalar_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_compare_actual_sparse<op_, InputValue_>(my_scalar);
    }

private:
    Scalar_ my_scalar;
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
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_compare_run_simple<op_>(input, length, my_scalar, output);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_compare_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_compare_run_simple<op_>(input_value, number, my_scalar, output_value);
    }

    OutputValue_ fill(bool, Index_) const {
        return delayed_compare<op_, InputValue_>(0, my_scalar);
    }
};

/**
 * Convenient alias for the scalar equality comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricEqualScalarHelper = DelayedUnaryIsometricCompareScalarHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "greater than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricGreaterThanScalarHelper = DelayedUnaryIsometricCompareScalarHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "less than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricLessThanScalarHelper = DelayedUnaryIsometricCompareScalarHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "greater than or equal" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricGreaterThanOrEqualScalarHelper = DelayedUnaryIsometricCompareScalarHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar "less than or equal" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricLessThanOrEqualScalarHelper = DelayedUnaryIsometricCompareScalarHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar non-equality comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricNotEqualScalarHelper = DelayedUnaryIsometricCompareScalarHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricEqualScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricGreaterThanScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricGreaterThanScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricLessThanScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricLessThanScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricGreaterThanOrEqualScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricGreaterThanOrEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricLessThanOrEqualScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricLessThanOrEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricNotEqualScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricNotEqualScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}
/**
 * @endcond
 */

/**
 * @brief Helper for delayed vector comparisons.
 *
 * This class compares each element of a `Matrix` to a row/column-specific value.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The comparison operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<CompareOperation op_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
class DelayedUnaryIsometricCompareVectorHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param vector Vector to use in the comparison with the matrix values.
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * The matrix value from each row/column is assumed to be on the left hand side of the comparison, while the corresponding value of `vector` is on the right.
     * @param by_row Whether `vector` corresponds to the rows.
     * If true, each element of the vector is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vector is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricCompareVectorHelper(Vector_ vector, bool by_row) : 
        my_vector(std::move(vector)),
        my_by_row(by_row)
    {
        for (auto x : my_vector) {
             if (!delayed_compare_actual_sparse<op_, InputValue_>(x)) {
                 my_sparse = false;
                 break;
             }
        }
    }

private:
    Vector_ my_vector;
    bool my_by_row;
    bool my_sparse = true;

public:
    std::optional<Index_> nrow() const {
        if (my_by_row) {
            return my_vector.size();
        } else {
            return std::nullopt;
        }
    }

    std::optional<Index_> ncol() const {
        if (my_by_row) {
            return std::nullopt;
        } else {
            return my_vector.size();
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
    void dense(bool row, Index_ idx, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, InputValue_>(input, length, my_vector[idx], output);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to enable optimizations.
            }
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_compare<op_, InputValue_>(input[i], my_vector[i + start]);
            }
        }
    }

    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, InputValue_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to enable optimizations.
            }
            Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_compare<op_, InputValue_>(input[i], my_vector[indices[i]]);
            }
        }
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input_value, const Index_* indices, OutputValue_* output_value) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, InputValue_>(input_value, number, my_vector[idx], output_value);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input_value = output_value; // basically an assertion to the compiler to enable optimizations.
            }
            for (Index_ i = 0; i < number; ++i) {
                output_value[i] = delayed_compare<op_, InputValue_>(input_value[i], my_vector[indices[i]]);
            }
        }
    }

    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            return delayed_compare<op_, InputValue_>(0, my_vector[idx]);
        } else {
            // We should only get to this point if it's sparse, otherwise no
            // single fill value would work across the length of my_vector.
            return 0;
        }
    }
};

/**
 * Convenient alias for the vector equality comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricEqualVectorHelper = DelayedUnaryIsometricCompareVectorHelper<CompareOperation::EQUAL, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector "greater than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricGreaterThanVectorHelper = DelayedUnaryIsometricCompareVectorHelper<CompareOperation::GREATER_THAN, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector "less than" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricLessThanVectorHelper = DelayedUnaryIsometricCompareVectorHelper<CompareOperation::LESS_THAN, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector "greater than or equal" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricGreaterThanOrEqualVectorHelper = DelayedUnaryIsometricCompareVectorHelper<CompareOperation::GREATER_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector "less than or equal" comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricLessThanOrEqualVectorHelper = DelayedUnaryIsometricCompareVectorHelper<CompareOperation::LESS_THAN_OR_EQUAL, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector non-equality comparison helper.
 *
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricNotEqualVectorHelper = DelayedUnaryIsometricCompareVectorHelper<CompareOperation::NOT_EQUAL, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricEqualVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricEqualVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricGreaterThanVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricGreaterThanVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricLessThanVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricLessThanVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricGreaterThanOrEqualVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricGreaterThanOrEqualVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricLessThanOrEqualVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricLessThanOrEqualVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricNotEqualVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricNotEqualVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}
/**
 * @endcond
 */

/**
 * @cond
 */
template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename Index_, typename OutputValue_>
void delayed_special_compare_run_simple(const InputValue_* input, Index_ length, OutputValue_* output) {
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_special_compare<op_, pass_, InputValue_>(input[i]);
    }
}

template<SpecialCompareOperation op_, bool pass_, typename InputValue_>
bool delayed_special_compare_actual_sparse() {
    return !delayed_special_compare<op_, pass_, InputValue_>(0);
}
/**
 * @endcond
 */

/**
 * @brief Delayed special value comparison.
 *
 * This class checks whether each element of a `Matrix` is one of the IEEE special values, e.g., NaN, Inf.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The special comparison operation.
 * @tparam pass_ Whether to return true if the special comparison is true.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<SpecialCompareOperation op_, bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricSpecialCompareHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * Default constructor.
     */
    DelayedUnaryIsometricSpecialCompareHelper() {
        my_sparse = !delayed_special_compare<op_, pass_, InputValue_>(0);
    }

private:
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
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_special_compare_run_simple<op_, pass_>(input, length, output);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_special_compare_run_simple<op_, pass_>(input, static_cast<Index_>(indices.size()), output);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_special_compare_run_simple<op_, pass_>(input_value, number, output_value);
    }

    OutputValue_ fill(bool, Index_) const {
        return !my_sparse;
    }
};

/**
 * Convenient alias for the "comparison to NaN" helper.
 *
 * @tparam pass_ Whether to return true if the input value is NaN.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricIsnanHelper = DelayedUnaryIsometricSpecialCompareHelper<SpecialCompareOperation::ISNAN, pass_, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "comparison to infinity" helper.
 *
 * @tparam pass_ Whether to return true if the input value is infinite.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricIsinfHelper = DelayedUnaryIsometricSpecialCompareHelper<SpecialCompareOperation::ISINF, pass_, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the "is finite" helper.
 *
 * @tparam pass_ Whether to return true if the input value is finite.
 * @tparam OutputValue_ Type of the result of the comparison.
 * @tparam InputValue_ Type of the matrix value used in the comparison.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<bool pass_, typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricIsfiniteHelper = DelayedUnaryIsometricSpecialCompareHelper<SpecialCompareOperation::ISFINITE, pass_, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
// Back-compatibility only.
template<bool pass_ = true, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricIsnan() {
    return std::make_shared<DelayedUnaryIsometricIsnanHelper<pass_, OutputValue_, InputValue_, Index_> >();
}

template<bool pass_ = true, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricIsinf() {
    return std::make_shared<DelayedUnaryIsometricIsinfHelper<pass_, OutputValue_, InputValue_, Index_> >();
}

template<bool pass_ = true, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricIsfinite() {
    return std::make_shared<DelayedUnaryIsometricIsfiniteHelper<pass_, OutputValue_, InputValue_, Index_> >();
}
/**
 * @endcond
 */

}

#endif
