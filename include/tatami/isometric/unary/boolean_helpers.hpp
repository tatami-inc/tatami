#ifndef TATAMI_ISOMETRIC_UNARY_BOOLEAN_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
#include "helper_interface.hpp"
#include <vector>
#include <type_traits>

/**
 * @file boolean_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric boolean operations.
 */

namespace tatami {

/**
 * @cond
 */
template<typename InputValue_, typename Index_, typename OutputValue_>
void delayed_boolean_cast(const InputValue_* input, Index_ length, OutputValue_* output) {
    if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
        input = output; // basically an assertion to the compiler to skip aliasing protection.
    }
    for (Index_ i = 0; i < length; ++i) {
        output[i] = static_cast<bool>(input[i]);
    }
}

template<typename InputValue_, typename Index_, typename OutputValue_>
void delayed_boolean_not(const InputValue_* input, Index_ length, OutputValue_* output) {
    if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
        input = output; // basically an assertion to the compiler to skip aliasing protection.
    }
    for (Index_ i = 0; i < length; ++i) {
        output[i] = !static_cast<bool>(input[i]);
    }
}

template<BooleanOperation op_, typename InputValue_, typename Index_, typename OutputValue_>
void delayed_boolean_run_simple(const InputValue_* input, Index_ length, bool scalar, OutputValue_* output) {
    if constexpr(op_ == BooleanOperation::AND) {
        if (scalar) {
            delayed_boolean_cast(input, length, output); 
        } else {
            std::fill_n(output, length, 0);
        }
    } else if constexpr(op_ == BooleanOperation::OR) {
        if (scalar) {
            std::fill_n(output, length, 1);
        } else {
            delayed_boolean_cast(input, length, output); 
        }
    } else if constexpr(op_ == BooleanOperation::XOR) {
        if (scalar) {
            delayed_boolean_not(input, length, output);
        } else {
            delayed_boolean_cast(input, length, output);
        }
    } else { // EQUAL
        if (scalar) {
            delayed_boolean_cast(input, length, output);
        } else {
            delayed_boolean_not(input, length, output);
        }
    }
}

template<BooleanOperation op_>
bool delayed_boolean_actual_sparse(bool scalar) {
    return delayed_boolean<op_>(0, scalar) == static_cast<bool>(0);
}
/**
 * @endcond
 */

/**
 * @brief Helper for delayed unary isometric scalar boolean operations.
 *
 * This class applies a boolean operation to each element of a `Matrix` where the other operand is a scalar.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<BooleanOperation op_, typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricBooleanScalarHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param scalar Scalar value.
     */
    DelayedUnaryIsometricBooleanScalarHelper(bool scalar) : my_scalar(scalar) {
        my_sparse = delayed_boolean_actual_sparse<op_>(my_scalar);
    }

private:
    bool my_scalar;
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
        delayed_boolean_run_simple<op_>(input, length, my_scalar, output);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_boolean_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_boolean_run_simple<op_>(input_value, number, my_scalar, output_value);
    }

    OutputValue_ fill(bool, Index_) const {
        // Remember, the operation needs to be performed on the InputValue_
        // to use in casting it to the OutputValue_.
        return delayed_boolean<op_>(0, my_scalar);
    }
};

/**
 * Convenient alias for the boolean equality scalar helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricBooleanEqualScalarHelper = DelayedUnaryIsometricBooleanScalarHelper<BooleanOperation::EQUAL, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the boolean AND scalar helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricBooleanAndScalarHelper = DelayedUnaryIsometricBooleanScalarHelper<BooleanOperation::AND, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the boolean OR scalar helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricBooleanOrScalarHelper = DelayedUnaryIsometricBooleanScalarHelper<BooleanOperation::OR, OutputValue_, InputValue_, Index_>;

/**
 * Convenient alias for the boolean XOR scalar helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
using DelayedUnaryIsometricBooleanXorScalarHelper = DelayedUnaryIsometricBooleanScalarHelper<BooleanOperation::XOR, OutputValue_, InputValue_, Index_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanEqualScalar(bool scalar) {
    return std::make_shared<DelayedUnaryIsometricBooleanEqualScalarHelper<OutputValue_, InputValue_, Index_> >(scalar);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanAndScalar(bool scalar) {
    return std::make_shared<DelayedUnaryIsometricBooleanAndScalarHelper<OutputValue_, InputValue_, Index_> >(scalar);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanOrScalar(bool scalar) {
    return std::make_shared<DelayedUnaryIsometricBooleanOrScalarHelper<OutputValue_, InputValue_, Index_> >(scalar);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanXorScalar(bool scalar) {
    return std::make_shared<DelayedUnaryIsometricBooleanXorScalarHelper<OutputValue_, InputValue_, Index_> >(scalar);
}
/**
 * @endcond
 */

/**
 * @brief Helper for a delayed unary isometric boolean NOT operation.
 *
 * This class applies a boolean NOT operation to each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricBooleanNotHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
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
        delayed_boolean_not(input, length, output);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_boolean_not(input, static_cast<Index_>(indices.size()), output);
    }

public:
    bool is_sparse() const {
        return false;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_boolean_not(input_value, number, output_value);
    }

    OutputValue_ fill(bool, Index_) const {
        return 1;
    }
};

/**
 * @brief Delayed unary isometric boolean cast.
 *
 * This class casts each element of a `Matrix` to a boolean 1 or 0, equivalent to the old `!!` trick.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricBooleanCastHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
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
        delayed_boolean_cast(input, length, output);
    }

    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_boolean_cast(input, static_cast<Index_>(indices.size()), output);
    }

public:
    bool is_sparse() const {
        return true;
    }

    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_boolean_cast(input_value, number, output_value);
    }

    OutputValue_ fill(bool, Index_) const {
        return 0;
    }
};

/**
 * @cond
 */
// Back-compatibility only.
typedef DelayedUnaryIsometricBooleanCastHelper<double, double, int> DelayedUnaryIsometricBooleanCast;
typedef DelayedUnaryIsometricBooleanNotHelper<double, double, int> DelayedUnaryIsometricBooleanNot;
/**
 * @endcond
 */

/**
 * @brief Helper for delayed unary isometric vector boolean operations.
 *
 * This class applies the specified boolean operation to each element of a `Matrix` where the other operand is a row/column-specific value.
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The boolean operation.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<BooleanOperation op_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
class DelayedUnaryIsometricBooleanVectorHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param vector Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * @param by_row Whether `vector` corresponds to the rows.
     * If true, each element of the vector is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vector is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricBooleanVectorHelper(Vector_ vector, bool by_row) : my_vector(std::move(vector)), my_by_row(by_row) {
        for (auto x : my_vector) {
             if (!delayed_boolean_actual_sparse<op_>(x)) {
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
            delayed_boolean_run_simple<op_>(input, length, my_vector[idx], output);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to skip aliasing protection.
            }
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_boolean<op_>(input[i], my_vector[i + start]);
            }
        }
    }

    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to skip aliasing protection.
            }
            Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_boolean<op_>(input[i], my_vector[indices[i]]);
            }
        }
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input_value, const Index_* index, OutputValue_* output_value) const {
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(input_value, number, my_vector[idx], output_value);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input_value = output_value; // basically an assertion to the compiler to skip aliasing protection.
            }
            for (Index_ i = 0; i < number; ++i) {
                output_value[i] = delayed_boolean<op_>(input_value[i], my_vector[index[i]]);
            }
        }
    }

    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            return delayed_boolean<op_>(0, my_vector[idx]);
        } else {
            // We should only get to this point if it's sparse, otherwise no
            // single fill value would work across the length of my_vector.
            return 0;
        }
    }
};

/**
 * Convenient alias for the boolean equality vector helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricBooleanEqualVectorHelper = DelayedUnaryIsometricBooleanVectorHelper<BooleanOperation::EQUAL, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the boolean AND vector helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricBooleanAndVectorHelper = DelayedUnaryIsometricBooleanVectorHelper<BooleanOperation::AND, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the boolean OR vector helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricBooleanOrVectorHelper = DelayedUnaryIsometricBooleanVectorHelper<BooleanOperation::OR, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the boolean XOR vector helper.
 *
 * @tparam OutputValue_ Type of the result of the boolean operation.
 * @tparam InputValue_ Type of the matrix value used in the boolean operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricBooleanXorVectorHelper = DelayedUnaryIsometricBooleanVectorHelper<BooleanOperation::XOR, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanEqualVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricBooleanEqualVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanAndVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricBooleanAndVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanOrVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricBooleanOrVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricBooleanXorVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricBooleanXorVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}
/**
 * @endcond
 */


}

#endif
