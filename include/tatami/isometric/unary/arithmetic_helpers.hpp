#ifndef TATAMI_ISOMETRIC_UNARY_ARITHMETIC_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_ARITHMETIC_HELPERS_H

#include "../arithmetic_utils.hpp"
#include <vector>
#include <limits>
#include <type_traits>

/**
 * @file arithmetic_helpers.hpp
 *
 * @brief Helper classes for delayed unary isometric arithmetic.
 */

namespace tatami {

/**
 * @cond
 */
template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Index_, typename Scalar_, typename OutputValue_>
void delayed_arithmetic_run_simple(const InputValue_* input, Index_ length, Scalar_ scalar, OutputValue_* output) {
    for (Index_ i = 0; i < length; ++i) {
        if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
            auto& val = output[i];
            val = delayed_arithmetic<op_, right_>(val, scalar);
        } else {
            output[i] = delayed_arithmetic<op_, right_>(input[i], scalar);
        }
    }
}

// The '*_actual_sparse' and '*_zero' functions should be mirrors of each other;
// we enforce this by putting their logic all in the same place.
template<bool check_only_, ArithmeticOperation op_, bool right_, typename InputValue_, typename Scalar_>
auto delayed_arithmetic_zeroish(Scalar_ scalar) {
    if constexpr(has_unsafe_divide_by_zero<op_, right_, InputValue_, Scalar_>()) {
        if constexpr(right_) {
            if (scalar) {
                auto val = delayed_arithmetic<op_, right_, InputValue_>(0, scalar);
                if constexpr(check_only_) {
                    return val == 0;
                } else {
                    return val;
                }
            }
        }

        if constexpr(check_only_) {
            return false;
        } else {
            throw std::runtime_error("division by zero is not supported");
            return static_cast<decltype(delayed_arithmetic<op_, right_, InputValue_>(0, scalar))>(1);
        }

    } else {
        auto val = delayed_arithmetic<op_, right_, InputValue_>(0, scalar);
        if constexpr(check_only_) {
            return val == 0;
        } else {
            return val;
        }
    }
}

template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Scalar_>
bool delayed_arithmetic_actual_sparse(Scalar_ scalar) {
    return delayed_arithmetic_zeroish<true, op_, right_, InputValue_, Scalar_>(scalar);
}

template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Scalar_>
auto delayed_arithmetic_zero(Scalar_ scalar) {
    return delayed_arithmetic_zeroish<false, op_, right_, InputValue_, Scalar_>(scalar);
}
/**
 * @endcond
 */

/**
 * @brief Delayed unary isometric scalar arithmetic.
 *
 * This class applies the specified arithmetic operation to each element of a `Matrix` where the other operand is a scalar.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the scalar should be on the right hand side of the arithmetic operation.
 * Ignored for commutative operations, e.g., `ADD` and `MULTIPLY`.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
 * @tparam Scalar_ Type of the scalar value.
 */
template<ArithmeticOperation op_, bool right_, typename InputValue_, typename Scalar_>
class DelayedUnaryIsometricArithmeticScalar {
public:
    /**
     * @param scalar Scalar value to be used in the operation.
     */
    DelayedUnaryIsometricArithmeticScalar(Scalar_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_arithmetic_actual_sparse<op_, right_, InputValue_>(scalar);
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
    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, Index_, Index_ length, const InputValue_* input, OutputValue_* output) const {
        delayed_arithmetic_run_simple<op_, right_>(input, length, my_scalar, output);
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool, Index_, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        delayed_arithmetic_run_simple<op_, right_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool, Index_, Index_ number, const InputValue_* input_value, const Index_*, OutputValue_* output_value) const {
        delayed_arithmetic_run_simple<op_, right_>(input_value, number, my_scalar, output_value);
    }

    template<typename OutputValue_, typename, typename Index_>
    OutputValue_ fill(bool, Index_) const {
        // We perform the operation with the InputValue_ before casting it to
        // the OutputValue_, which is consistent with the behavior of all other
        // methods. See ../arithmetic_utils.hpp for some comments about the
        // safety of this cast when the value is known at compile time.
        return delayed_arithmetic_zero<op_, right_, InputValue_>(my_scalar);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed unary isometric vector arithmetic.
 *
 * This class applies the specified arithmetic operation to each element of a `Matrix` where the other operand is row/column-specific value.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the vector's values should be on the right hand side of the arithmetic operation.
 * Ignored for some `op_`.
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input, length, my_vector[idx], output);
        } else {
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_arithmetic<op_, right_>(val, my_vector[i + start]);
                } else {
                    output[i] = delayed_arithmetic<op_, right_>(input[i], my_vector[i + start]);
                }
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void dense(bool row, Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output[i];
                    val = delayed_arithmetic<op_, right_>(val, my_vector[indices[i]]);
                } else {
                    output[i] = delayed_arithmetic<op_, right_>(input[i], my_vector[indices[i]]);
                }
            }
        }
    }

    template<typename Index_, typename OutputValue_>
    void sparse(bool row, Index_ idx, Index_ number, const InputValue_* input_value, const Index_* indices, OutputValue_* output_value) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input_value, number, my_vector[idx], output_value);
        } else {
            for (Index_ i = 0; i < number; ++i) {
                if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                    auto& val = output_value[i];
                    val = delayed_arithmetic<op_, right_>(val, my_vector[indices[i]]);
                } else {
                    output_value[i] = delayed_arithmetic<op_, right_>(input_value[i], my_vector[indices[i]]);
                }
            }
        }
    }

    template<typename OutputValue_, typename, typename Index_>
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
 * @tparam InputValue_ Type of the matrix value to use in the operation.
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
