#ifndef TATAMI_ISOMETRIC_UNARY_COMPARE_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_COMPARE_HELPERS_H

#include "../compare_utils.hpp"
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
template<CompareOperation op_, typename InputValue_, typename Index_, typename OutputValue_>
void delayed_compare_run_simple(const InputValue_* input, Index_ length, InputValue_ scalar, OutputValue_* output) {
    for (Index_ i = 0; i < length; ++i) {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            auto& val = output[i];
            val = delayed_compare<op_>(val, scalar);
        } else {
            output[i] = delayed_compare<op_>(input[i], scalar);
        }
    }
}

template<CompareOperation op_, typename InputValue_>
bool delayed_compare_actual_sparse(InputValue_ scalar) {
    return !delayed_compare<op_, InputValue_>(0, scalar);
}
/**
 * @endcond
 */

/**
 * @brief Delayed scalar comparison.
 *
 * This class compares each element of a `Matrix` to a scalar.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The comparison operation.
 * @tparam InputValue_ Type of the matrix value to use in the comparison.
 */
template<CompareOperation op_, typename InputValue_>
class DelayedUnaryIsometricCompareScalar {
public:
    /**
     * @param scalar Scalar to be compared to the matrix values.
     * The matrix value is assumed to be on the left hand side of the comparison, while `scalar` is on the right.
     */
    DelayedUnaryIsometricCompareScalar(InputValue_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_compare_actual_sparse<op_>(my_scalar);
    }

private:
    InputValue_ my_scalar;
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
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_compare_run_simple<op_>(input, length, my_scalar, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_compare_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_compare_run_simple<op_>(input_value, number, my_scalar, output_value);
    }

    template<typename OutputValue_, typename, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return delayed_compare<op_, InputValue_>(0, my_scalar);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector comparisons.
 *
 * This class compares each element of a `Matrix` to a row/column-specific value.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The comparison operation.
 * @tparam InputValue_ Type of the matrix value to use in the comparison.
 * @tparam Vector_ Type of the vector.
 */
template<CompareOperation op_, typename InputValue_, typename Vector_>
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
    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, InputValue_>(input, length, my_vector[idx], output);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_compare<op_, InputValue_>(val, my_vector[i + start]);
                } else {
                    output[i] = delayed_compare<op_, InputValue_>(input[i], my_vector[i + start]);
                }
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, InputValue_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_compare<op_, InputValue_>(val, my_vector[indices[i]]);
                } else {
                    output[i] = delayed_compare<op_, InputValue_>(input[i], my_vector[indices[i]]);
                }
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input, const Index_* indices, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_compare_run_simple<op_, InputValue_>(input, number, my_vector[idx], output);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_compare<op_, InputValue_>(val, my_vector[indices[i]]);
                } else {
                    output[i] = delayed_compare<op_, InputValue_>(input[i], my_vector[indices[i]]);
                }
            }
        }
    }

    template<typename OutputValue_, typename, typename Index_>
    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            return delayed_compare<op_, InputValue_>(0, my_vector[idx]);
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Value to be compared.
 * @return A helper class for a delayed equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::EQUAL, InputValue_> make_DelayedUnaryIsometricEqualScalar(InputValue_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::EQUAL, InputValue_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed greater-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN, InputValue_> make_DelayedUnaryIsometricGreaterThanScalar(InputValue_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN, InputValue_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed less-than comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN, InputValue_> make_DelayedUnaryIsometricLessThanScalar(InputValue_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN, InputValue_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed greater-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN_OR_EQUAL, InputValue_> make_DelayedUnaryIsometricGreaterThanOrEqualScalar(InputValue_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::GREATER_THAN_OR_EQUAL, InputValue_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed less-than-or-equal comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN_OR_EQUAL, InputValue_> make_DelayedUnaryIsometricLessThanOrEqualScalar(InputValue_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::LESS_THAN_OR_EQUAL, InputValue_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be compared.
 * @return A helper class for a delayed non-equality comparison to a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricCompareScalar<CompareOperation::NOT_EQUAL, InputValue_> make_DelayedUnaryIsometricNotEqualScalar(InputValue_ scalar) {
    return DelayedUnaryIsometricCompareScalar<CompareOperation::NOT_EQUAL, InputValue_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed equality comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricCompareVector<CompareOperation::EQUAL, InputValue_, Vector_> make_DelayedUnaryIsometricEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::EQUAL, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed greater-than comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN, InputValue_, Vector_> make_DelayedUnaryIsometricGreaterThanVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed less-than comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN, InputValue_, Vector_> make_DelayedUnaryIsometricLessThanVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed greater-than-or-equal comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN_OR_EQUAL, InputValue_, Vector_> make_DelayedUnaryIsometricGreaterThanOrEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::GREATER_THAN_OR_EQUAL, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed less-than-or-equal comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN_OR_EQUAL, InputValue_, Vector_> make_DelayedUnaryIsometricLessThanOrEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::LESS_THAN_OR_EQUAL, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be compared.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricCompareVector`.
 * @return A helper class for a delayed non-equality comparison to a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricCompareVector<CompareOperation::NOT_EQUAL, InputValue_, Vector_> make_DelayedUnaryIsometricNotEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricCompareVector<CompareOperation::NOT_EQUAL, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @cond
 */
template<SpecialCompareOperation op_, bool pass_, typename InputValue_, typename Index_>
void delayed_special_compare_run_simple(InputValue_* buffer, Index_ length) {
    for (Index_ i = 0; i < length; ++i) {
        auto& val = buffer[i];
        val = delayed_special_compare<op_, pass_, InputValue_>(val);
    }
}

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
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The special comparison operation.
 * @tparam pass_ Whether to return true if the special comparison is true.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<SpecialCompareOperation op_, bool pass_, typename InputValue_>
class DelayedUnaryIsometricSpecialCompare {
public:
    DelayedUnaryIsometricSpecialCompare() {
        my_sparse = !delayed_special_compare<op_, pass_, InputValue_>(0);
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
    void dense(bool, Index_, Index_, Index_ length, InputValue_* buffer) const {
        delayed_special_compare_run_simple<op_, pass_>(buffer, length);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_special_compare_run_simple<op_, pass_>(input, length, output);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, InputValue_* buffer) const {
        delayed_special_compare_run_simple<op_, pass_>(buffer, static_cast<Index_>(indices.size()));
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_special_compare_run_simple<op_, pass_>(input, static_cast<Index_>(indices.size()), output);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, InputValue_* buffer, const Index_*) const {
        delayed_special_compare_run_simple<op_, pass_>(buffer, number);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input, const Index_*, OutputValue_* output) const {
        delayed_special_compare_run_simple<op_, pass_>(input, number, output);
    }

    template<typename OutputValue_, typename, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return !my_sparse;
    }
    /**
     * @endcond
     */
};

/**
 * @tparam pass_ Whether to return truthy if the matrix value is NaN.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @return A helper class for a delayed NaN check,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool pass_ = true, typename InputValue_ = double>
DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISNAN, pass_, InputValue_> make_DelayedUnaryIsometricIsnan() {
    return DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISNAN, pass_, InputValue_>();
}

/**
 * @tparam pass_ Whether to return truthy if the matrix value is infinite.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @return A helper class for a delayed check for infinity (positive or negative),
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool pass_ = true, typename InputValue_ = double>
DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISINF, pass_, InputValue_> make_DelayedUnaryIsometricIsinf() {
    return DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISINF, pass_, InputValue_>();
}

/**
 * @tparam pass_ Whether to return truthy if the matrix value is finite.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @return A helper class for a delayed check for finite values,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool pass_ = true, typename InputValue_ = double>
DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISFINITE, pass_, InputValue_> make_DelayedUnaryIsometricIsfinite() {
    return DelayedUnaryIsometricSpecialCompare<SpecialCompareOperation::ISFINITE, pass_, InputValue_>();
}

}

#endif
