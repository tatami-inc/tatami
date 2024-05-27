#ifndef TATAMI_ISOMETRIC_UNARY_ARITHMETIC_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_ARITHMETIC_HELPERS_H

#include "../arithmetic_utils.hpp"
#include <vector>
#include <limits>

/**
 * @file arithmetic_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric arithmetic.
 */

namespace tatami {

/**
 * @cond
 */
template<ArithmeticOperation op_, bool right_, typename Scalar_, typename Value_, typename Index_>
void delayed_arithmetic_run_simple(Value_* buffer, Index_ length, Scalar_ scalar) {
    for (Index_ i = 0; i < length; ++i) {
        auto& val = buffer[i];
        val = delayed_arithmetic<op_, right_>(val, scalar);
    }
}

template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Index_, typename Scalar_, typename OutputValue_>
void delayed_arithmetic_run_simple(const InputValue_* input, Index_ length, Scalar_ scalar, OutputValue_* output) {
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_arithmetic<op_, right_>(input[i], scalar);
    }
}

template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
constexpr bool delayed_arithmetic_unsupported_division_by_zero() {
    return !std::numeric_limits<Value_>::is_iec559 && op_ == ArithmeticOperation::DIVIDE && !right_;
}

template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
bool delayed_arithmetic_actual_sparse(Scalar_ scalar) {
    if constexpr(delayed_arithmetic_unsupported_division_by_zero<op_, right_, Value_, Scalar_>()) {
        // If we didn't catch this case, the else() condition would be dividing
        // by zero in a Value_ that doesn't support it, and that would be visible
        // at compile time - possibly resulting in compiler warnings. So we
        // declare that this is always non-sparse, and hope that the equivalent
        // zero() method doesn't get called.
        return false;
    } else {
        // Empirically testing this, to accommodate special values (e.g., NaN, Inf) for scalars.
        return delayed_arithmetic<op_, right_, Value_>(0, scalar) == 0;
    }
}

template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
Value_ delayed_arithmetic_zero(Scalar_ scalar) {
    if constexpr(delayed_arithmetic_unsupported_division_by_zero<op_, right_, Value_, Scalar_>()) {
        // Avoid potential problems with division by zero that can be detected
        // at compile time (e.g., resulting in unnecessary compiler warnings).
        throw std::runtime_error("division by zero is not supported with IEEE-754 floats");
    } else {
        return delayed_arithmetic<op_, right_, Value_>(0, scalar);
    }
}
/**
 * @endcond
 */

/**
 * @brief Delayed unary isometric scalar arithmetic.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the scalar should be on the right hand side of the arithmetic operation.
 * Ignored for commutative operations, e.g., `ADD` and `MULTIPLY`.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar value.
 */
template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Scalar_>
class DelayedUnaryIsometricArithmeticScalar {
public:
    /**
     * @param scalar Scalar value to be used in the operation.
     */
    DelayedUnaryIsometricArithmeticScalar(Scalar_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_arithmetic_actual_sparse<op_, right_, InputValue_>(my_scalar);
    }

private:
    Scalar_ my_scalar;
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
        delayed_arithmetic_run_simple<op_, right_>(buffer, length, my_scalar);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_arithmetic_run_simple<op_, right_>(input, length, my_scalar, output);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, InputValue_* buffer) const {
        delayed_arithmetic_run_simple<op_, right_>(buffer, static_cast<Index_>(indices.size()), my_scalar);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_arithmetic_run_simple<op_, right_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

    template<typename Index_>
    void sparse(bool, Index_, Index_ number, InputValue_* value, const Index_*) const {
        delayed_arithmetic_run_simple<op_, right_>(value, number, my_scalar);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_arithmetic_run_simple<op_, right_>(input_value, number, my_scalar, output_value);
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // We perform the operation with the InputValue_ before casting it to
        // the OutputValue_, which is consistent with the behavior of all other
        // methods. This has some interesting implications if only one or the
        // other supports, e.g., division by zero, but that's not my problem.
        return delayed_arithmetic_zero<op_, right_, InputValue_>(my_scalar);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed unary isometric vector arithmetic.
 *
 * This should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the vector's values should be on the right hand side of the arithmetic operation.
 * Ignored for some `op_`.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 */
template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Vector_>
class DelayedUnaryIsometricArithmeticVector {
public:
    /**
     * @param vector Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `by_row = true`, otherwise it should be of length equal to the number of columns.
     * @param by_row Whether `vector` corresponds to the rows.
     * If true, each element of the vector is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vector is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricArithmeticVector(Vector_ vector, bool by_row) : my_vector(std::move(vector)), my_by_row(by_row) {
        for (auto x : my_vector) {
            if (!delayed_arithmetic_actual_sparse<op_, right_, InputValue_>(x)) {
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
    template<typename Index_>
    void dense(bool row, Index_ idx, Index_ start, Index_ length, InputValue_* buffer) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(buffer, length, my_vector[idx]);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                auto& val = buffer[i];
                val = delayed_arithmetic<op_, right_>(val, my_vector[i + start]);
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input, length, my_vector[idx], output);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_arithmetic<op_, right_>(input[i], my_vector[i + start]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, InputValue_* buffer) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(buffer, static_cast<Index_>(indices.size()), my_vector[idx]);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                auto& val = buffer[i];
                val = delayed_arithmetic<op_, right_>(val, my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                output[i] = delayed_arithmetic<op_, right_>(input[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, InputValue_* value, const Index_* indices) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(value, number, my_vector[idx]);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                auto& val = value[i];
                val = delayed_arithmetic<op_, right_>(val, my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input_value, const Index_* indices, OutputValue_* output_value) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input_value, number, my_vector[idx], output_value);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                output_value[i] = delayed_arithmetic<op_, right_>(input_value[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename OutputValue_, typename Index_>
    OutputValue_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            return delayed_arithmetic_zero<op_, right_, InputValue_>(my_vector[idx]);
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
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be added.
 * @return A helper class for delayed scalar addition,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::ADD, true, InputValue_, Scalar_> make_DelayedUnaryIsometricAddScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::ADD, true, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be subtracted.
 * @return A helper class for delayed scalar subtraction,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::SUBTRACT, right_, InputValue_, Scalar_> make_DelayedUnaryIsometricSubtractScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::SUBTRACT, right_, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be multiplied.
 * @return A helper class for delayed scalar multiplication,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MULTIPLY, true, InputValue_, Scalar_> make_DelayedUnaryIsometricMultiplyScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MULTIPLY, true, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be divided.
 * @return A helper class for delayed scalar division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::DIVIDE, right_, InputValue_, Scalar_> make_DelayedUnaryIsometricDivideScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::DIVIDE, right_, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the power transformation.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be power transformed.
 * @return A helper class for delayed scalar power transformation,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::POWER, right_, InputValue_, Scalar_> make_DelayedUnaryIsometricPowerScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::POWER, right_, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the modulus.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be modulo transformed.
 * @return A helper class for delayed scalar modulus,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MODULO, right_, InputValue_, Scalar_> make_DelayedUnaryIsometricModuloScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MODULO, right_, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the integer division.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be integer divided.
 * @return A helper class for delayed scalar integer division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Scalar_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::INTEGER_DIVIDE, right_, InputValue_, Scalar_> make_DelayedUnaryIsometricIntegerDivideScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::INTEGER_DIVIDE, right_, InputValue_, Scalar_>(std::move(scalar));
}

/**
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to be added to the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector addition,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::ADD, true, InputValue_, Vector_> make_DelayedUnaryIsometricAddVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::ADD, true, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to subtract from (or be subtracted by) the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector subtraction,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::SUBTRACT, right_, InputValue_, Vector_> make_DelayedUnaryIsometricSubtractVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::SUBTRACT, right_, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to multiply the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector multiplication,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MULTIPLY, true, InputValue_, Vector_> make_DelayedUnaryIsometricMultiplyVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MULTIPLY, true, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to divide (or be divided by) the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::DIVIDE, right_, InputValue_, Vector_> make_DelayedUnaryIsometricDivideVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::DIVIDE, right_, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the power transformation.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to use in the power transformation of the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector power transformation,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::POWER, right_, InputValue_, Vector_> make_DelayedUnaryIsometricPowerVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::POWER, right_, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the modulus.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to use in the modulus of the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector modulus,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MODULO, right_, InputValue_, Vector_> make_DelayedUnaryIsometricModuloVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MODULO, right_, InputValue_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the integer division.
 * @tparam InputValue_ Type of the matrix value before the operation.
 * @tparam Vector_ Type of the vector.
 * @param vector Vector to integer divide (or be integer divided by) the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename InputValue_ = double, typename Vector_>
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::INTEGER_DIVIDE, right_, InputValue_, Vector_> make_DelayedUnaryIsometricIntegerDivideVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::INTEGER_DIVIDE, right_, InputValue_, Vector_>(std::move(vector), by_row);
}

}

#endif
