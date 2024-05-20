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
void delayed_arithmetic_run_simple(Scalar_ scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        delayed_arithmetic_run<op_, right_>(buffer[i], scalar);
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
        Value_ output = 0;
        delayed_arithmetic_run<op_, right_>(output, scalar);
        return output == 0;
    }
}

template<ArithmeticOperation op_, bool right_, typename Value_, typename Scalar_>
Value_ delayed_arithmetic_zero(Scalar_ scalar) {
    if constexpr(delayed_arithmetic_unsupported_division_by_zero<op_, right_, Value_, Scalar_>()) {
        // Avoid potential problems with division by zero that can be detected
        // at compile time (e.g., resulting in unnecessary compiler warnings).
        throw std::runtime_error("division by zero is not supported with IEEE-754 floats");
    } else {
        Value_ output = 0;
        delayed_arithmetic_run<op_, right_>(output, scalar);
        return output;
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
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar value.
 */
template<ArithmeticOperation op_, bool right_, typename Value_ = double, typename Scalar_ = Value_>
class DelayedUnaryIsometricArithmeticScalar {
public:
    /**
     * @param scalar Scalar value to be used in the operation.
     */
    DelayedUnaryIsometricArithmeticScalar(Scalar_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_arithmetic_actual_sparse<op_, right_, Value_>(my_scalar);
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
        delayed_arithmetic_run_simple<op_, right_>(my_scalar, length, buffer);
    }

    template<typename Index_>
    void dense(bool, Index_, const std::vector<Index_>& indices, Value_* buffer) const {
        delayed_arithmetic_run_simple<op_, right_>(my_scalar, indices.size(), buffer);
    }


    template<typename Index_>
    void sparse(bool, Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_arithmetic_run_simple<op_, right_>(my_scalar, number, buffer);
    }

    template<typename Index_>
    Value_ fill(bool, Index_) const {
        return delayed_arithmetic_zero<op_, right_, Value_>(my_scalar);
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
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 */
template<ArithmeticOperation op_, bool right_, typename Value_ = double, typename Vector_ = std::vector<double> >
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
            if (!delayed_arithmetic_actual_sparse<op_, right_, Value_>(x)) {
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
            delayed_arithmetic_run_simple<op_, right_>(my_vector[idx], length, buffer);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                delayed_arithmetic_run<op_, right_>(buffer[i], my_vector[i + start]);
            }
        }
    }

    template<typename Index_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, Value_* buffer) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(my_vector[idx], indices.size(), buffer);
        } else {
            for (Index_ i = 0, length = indices.size(); i < length; ++i) {
                delayed_arithmetic_run<op_, right_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    void sparse(bool row, Index_ idx, Index_ number, Value_* buffer, const Index_* indices) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(my_vector[idx], number, buffer);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                delayed_arithmetic_run<op_, right_>(buffer[i], my_vector[indices[i]]);
            }
        }
    }

    template<typename Index_>
    Value_ fill(bool row, Index_ idx) const {
        if (row == my_by_row) {
            return delayed_arithmetic_zero<op_, right_, Value_>(my_vector[idx]);
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
 * @param scalar Scalar value to be added.
 * @return A helper class for delayed scalar addition,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::ADD, true, Value_, Scalar_> make_DelayedUnaryIsometricAddScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::ADD, true, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be subtracted.
 * @return A helper class for delayed scalar subtraction,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::SUBTRACT, right_, Value_, Scalar_> make_DelayedUnaryIsometricSubtractScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::SUBTRACT, right_, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be multiplied.
 * @return A helper class for delayed scalar multiplication,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MULTIPLY, true, Value_, Scalar_> make_DelayedUnaryIsometricMultiplyScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MULTIPLY, true, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be divided.
 * @return A helper class for delayed scalar division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::DIVIDE, right_, Value_, Scalar_> make_DelayedUnaryIsometricDivideScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::DIVIDE, right_, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the power transformation.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be power transformed.
 * @return A helper class for delayed scalar power transformation,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::POWER, right_, Value_, Scalar_> make_DelayedUnaryIsometricPowerScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::POWER, right_, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the modulus.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be modulo transformed.
 * @return A helper class for delayed scalar modulus,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MODULO, right_, Value_, Scalar_> make_DelayedUnaryIsometricModuloScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::MODULO, right_, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the integer division.
 * @tparam Value_ Type of the data value.
 * @tparam Scalar_ Type of the scalar.
 * @param scalar Scalar value to be integer divided.
 * @return A helper class for delayed scalar integer division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Scalar_ = Value_>
DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::INTEGER_DIVIDE, right_, Value_, Scalar_> make_DelayedUnaryIsometricIntegerDivideScalar(Scalar_ scalar) {
    return DelayedUnaryIsometricArithmeticScalar<ArithmeticOperation::INTEGER_DIVIDE, right_, Value_, Scalar_>(std::move(scalar));
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to be added to the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector addition,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::ADD, true, Value_, Vector_> make_DelayedUnaryIsometricAddVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::ADD, true, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to subtract from (or be subtracted by) the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector subtraction,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::SUBTRACT, right_, Value_, Vector_> make_DelayedUnaryIsometricSubtractVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::SUBTRACT, right_, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to multiply the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector multiplication,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MULTIPLY, true, Value_, Vector_> make_DelayedUnaryIsometricMultiplyVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MULTIPLY, true, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to divide (or be divided by) the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::DIVIDE, right_, Value_, Vector_> make_DelayedUnaryIsometricDivideVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::DIVIDE, right_, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the power transformation.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to use in the power transformation of the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector power transformation,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::POWER, right_, Value_, Vector_> make_DelayedUnaryIsometricPowerVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::POWER, right_, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the modulus.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to use in the modulus of the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector modulus,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MODULO, right_, Value_, Vector_> make_DelayedUnaryIsometricModuloVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::MODULO, right_, Value_, Vector_>(std::move(vector), by_row);
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the integer division.
 * @tparam Value_ Type of the data value.
 * @tparam Vector_ Type of the vector.
 *
 * @param vector Vector to integer divide (or be integer divided by) the rows/columns.
 * @param by_row Whether each element of `vector` corresponds to a row, see `DelayedUnaryIsometricArithmeticVector`.
 * @return A helper class for delayed vector division,
 * to be used as the `operation` in a `DelayedUnaryIsometricOperation`.
 */
template<bool right_, typename Value_ = double, typename Vector_ = std::vector<double> >
DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::INTEGER_DIVIDE, right_, Value_, Vector_> make_DelayedUnaryIsometricIntegerDivideVector(Vector_ vector, bool by_row) {
    return DelayedUnaryIsometricArithmeticVector<ArithmeticOperation::INTEGER_DIVIDE, right_, Value_, Vector_>(std::move(vector), by_row);
}

}

#endif
