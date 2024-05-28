#ifndef TATAMI_ISOMETRIC_UNARY_BOOLEAN_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_BOOLEAN_HELPERS_H

#include "../boolean_utils.hpp"
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
template<BooleanOperation op_, typename InputValue_, typename Index_, typename OutputValue_>
void delayed_boolean_run_simple(const InputValue_* input, Index_ length, bool scalar, OutputValue_* output) {
    for (Index_ i = 0; i < length; ++i) {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            auto& val = output[i];
            val = delayed_boolean<op_>(val, scalar);
        } else {
            output[i] = delayed_boolean<op_>(input[i], scalar);
        }
    }
}

template<BooleanOperation op_, typename InputValue_>
bool delayed_boolean_actual_sparse(bool scalar) {
    return delayed_boolean<op_>(0, scalar) == static_cast<bool>(0);
}
/**
 * @endcond
 */

/**
 * @brief Delayed unary isometric scalar boolean operation.
 *
 * This class applies a boolean operation to each element of a `Matrix` where the other operand is a scalar.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The boolean operation.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<BooleanOperation op_, typename InputValue_>
class DelayedUnaryIsometricBooleanScalar {
public:
    /**
     * @param scalar Scalar value.
     */
    DelayedUnaryIsometricBooleanScalar(bool scalar) : my_scalar(scalar) {
        my_sparse = delayed_boolean_actual_sparse<op_, InputValue_>(my_scalar);
    }

private:
    bool my_scalar;
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
        delayed_boolean_run_simple<op_>(input, length, my_scalar, output);
    }

    template<typename Index_, typename OutputValue_> 
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_boolean_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_boolean_run_simple<op_>(input_value, number, my_scalar, output_value);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // Remember, the operation needs to be performed on the InputValue_
        // to use in casting it to the OutputValue_.
        return delayed_boolean<op_>(0, my_scalar);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed unary isometric boolean NOT operation.
 *
 * This class applies a boolean NOT operation to each element of a `Matrix`.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 */
template<typename InputValue_ = double>
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
    template<typename Index_, typename OutputValue_>
    void core(const InputValue_* input, Index_ length, OutputValue_* output) const {
        for (Index_ i = 0; i < length; ++i) {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                auto& val = output[i];
                val = !static_cast<bool>(val);
            } else {
                output[i] = !static_cast<bool>(input[i]);
            }
        }
    }

public:
    /**
     * @cond
     */
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        core(input, length, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        core(input, static_cast<Index_>(indices.size()), output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        core(input_value, number, output_value);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        return 1;
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed unary isometric vector boolean operations.
 *
 * This class applies the specified boolean operation to each element of a `Matrix` where the other operand is a row/column-specific value.
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The boolean operation.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 */
template<BooleanOperation op_, typename InputValue_, typename Vector_>
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
             if (!delayed_boolean_actual_sparse<op_, InputValue_>(x)) {
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
            delayed_boolean_run_simple<op_>(input, length, my_vector[idx], output);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_boolean<op_>(val, my_vector[i + start]);
                } else {
                    output[i] = delayed_boolean<op_>(input[i], my_vector[i + start]);
                }
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_boolean<op_>(val, my_vector[indices[i]]);
                } else {
                    output[i] = delayed_boolean<op_>(input[i], my_vector[indices[i]]);
                }
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input_value, const Index_* index, OutputValue_* output_value) const {
        if (row == my_by_row) {
            delayed_boolean_run_simple<op_>(input_value, number, my_vector[idx], output_value);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output_value[i];
                    val = delayed_boolean<op_>(val, my_vector[index[i]]);
                } else {
                    output_value[i] = delayed_boolean<op_>(input_value[i], my_vector[index[i]]);
                }
            }
        }
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            return delayed_boolean<op_>(0, my_vector[idx]);
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
 * @return A helper class for a delayed NOT operation,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricBooleanNot<InputValue_> make_DelayedUnaryIsometricBooleanNot() {
    return DelayedUnaryIsometricBooleanNot<InputValue_>();
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to use in the operation.
 * @return A helper class for a delayed AND operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::AND, InputValue_> make_DelayedUnaryIsometricBooleanAndScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::AND, InputValue_>(scalar);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to use in the operation.
 * @return A helper class for a delayed OR operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::OR, InputValue_> make_DelayedUnaryIsometricBooleanOrScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::OR, InputValue_>(scalar);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be used in the operation.
 * @return A helper class for a delayed XOR operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::XOR, InputValue_> make_DelayedUnaryIsometricBooleanXorScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::XOR, InputValue_>(scalar);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @param scalar Scalar value to be used in the operation.
 * @return A helper class for a delayed boolean equality operation with a scalar,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double>
DelayedUnaryIsometricBooleanScalar<BooleanOperation::EQUAL, InputValue_> make_DelayedUnaryIsometricBooleanEqualScalar(bool scalar) {
    return DelayedUnaryIsometricBooleanScalar<BooleanOperation::EQUAL, InputValue_>(scalar);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed AND operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricBooleanVector<BooleanOperation::AND, InputValue_, Vector_> make_DelayedUnaryIsometricBooleanAndVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::AND, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed OR operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricBooleanVector<BooleanOperation::OR, InputValue_, Vector_> make_DelayedUnaryIsometricBooleanOrVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::OR, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed XOR operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricBooleanVector<BooleanOperation::XOR, InputValue_, Vector_> make_DelayedUnaryIsometricBooleanXorVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::XOR, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector of values to be used in the operation.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricBooleanVector`.
 * @return A helper class for a delayed boolean equality operation with a vector,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricBooleanVector<BooleanOperation::EQUAL, InputValue_, Vector_> make_DelayedUnaryIsometricBooleanEqualVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricBooleanVector<BooleanOperation::EQUAL, InputValue_, Vector_>(std::move(vector), by_row);
}

}

#endif
