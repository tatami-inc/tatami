#ifndef TATAMI_ISOMETRIC_UNARY_ARITHMETIC_HELPERS_H
#define TATAMI_ISOMETRIC_UNARY_ARITHMETIC_HELPERS_H

#include "../arithmetic_utils.hpp"
#include "helper_interface.hpp"
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
void delayed_arithmetic_run_simple(const InputValue_* input, const Index_ length, const Scalar_ scalar, OutputValue_* const output) {
    if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
        input = output; // basically an assertion to the compiler to skip aliasing protection.
    }
    for (Index_ i = 0; i < length; ++i) {
        output[i] = delayed_arithmetic<op_, right_>(input[i], scalar);
    }
}

// The '*_actual_sparse' and '*_zero' functions should be mirrors of each other;
// we enforce this by putting their logic all in the same place.
template<bool check_only_, ArithmeticOperation op_, bool right_, typename OutputValue_, typename InputValue_, typename Scalar_>
auto delayed_arithmetic_zeroish(const Scalar_ scalar) {
    if constexpr(has_unsafe_divide_by_zero<op_, right_, InputValue_, Scalar_>()) {
        if constexpr(right_) {
            if (scalar) {
                const OutputValue_ val = delayed_arithmetic<op_, right_, InputValue_, Scalar_>(0, scalar);
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
            return 0;
        }

    } else {
        const OutputValue_ val = delayed_arithmetic<op_, right_, InputValue_, Scalar_>(0, scalar);
        if constexpr(check_only_) {
            return val == 0;
        } else {
            return val;
        }
    }
}

template<ArithmeticOperation op_, bool right_, typename OutputValue_, typename InputValue_, typename Scalar_>
bool delayed_arithmetic_actual_sparse(const Scalar_ scalar) {
    return delayed_arithmetic_zeroish<true, op_, right_, OutputValue_, InputValue_, Scalar_>(scalar);
}

template<ArithmeticOperation op_, bool right_, typename OutputValue_, typename InputValue_, typename Scalar_>
OutputValue_ delayed_arithmetic_zero(const Scalar_ scalar) {
    return delayed_arithmetic_zeroish<false, op_, right_, OutputValue_, InputValue_, Scalar_>(scalar);
}
/**
 * @endcond
 */

/**
 * @brief Helper for delayed unary isometric scalar arithmetic.
 *
 * This class applies the specified arithmetic operation to each element of a `Matrix` where the other operand is a scalar.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the scalar should be on the right hand side of the arithmetic operation.
 * Ignored for commutative operations, e.g., `op_ == ADD` and `MULTIPLY`.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<ArithmeticOperation op_, bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
class DelayedUnaryIsometricArithmeticScalarHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param scalar Scalar value to be used in the operation.
     */
    DelayedUnaryIsometricArithmeticScalarHelper(const Scalar_ scalar) : my_scalar(scalar) {
        my_sparse = delayed_arithmetic_actual_sparse<op_, right_, OutputValue_, InputValue_>(my_scalar);
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
    void dense(const bool, const Index_, const Index_, const Index_ length, const InputValue_* const input, OutputValue_* const output) const {
        delayed_arithmetic_run_simple<op_, right_>(input, length, my_scalar, output);
    }

    void dense(const bool, const Index_, const std::vector<Index_>& indices, const InputValue_* const input, OutputValue_* const output) const {
        delayed_arithmetic_run_simple<op_, right_>(input, static_cast<Index_>(indices.size()), my_scalar, output);
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(const bool, const Index_, const Index_ number, const InputValue_* const input_value, const Index_* const, OutputValue_* const output_value) const {
        delayed_arithmetic_run_simple<op_, right_>(input_value, number, my_scalar, output_value);
    }

    OutputValue_ fill(const bool, const Index_) const {
        // We perform the operation with the InputValue_ before casting it to
        // the OutputValue_, which is consistent with the behavior of all other
        // methods. See ../arithmetic_utils.hpp for some comments about the
        // safety of this cast when the value is known at compile time.
        return delayed_arithmetic_zero<op_, right_, OutputValue_, InputValue_>(my_scalar);
    }
};

/**
 * Convenient alias for the scalar addition helper.
 *
 * @tparam OutputValue_ Type of the result of the addition.
 * @tparam InputValue_ Type of the matrix value used in the addition.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricAddScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::ADD, true, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar subtraction helper.
 *
 * @tparam right_ Whether the scalar should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the subtraction.
 * @tparam InputValue_ Type of the matrix value used in the subtraction.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricSubtractScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::SUBTRACT, right_, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar multiplication helper.
 *
 * @tparam OutputValue_ Type of the result of the multiplication.
 * @tparam InputValue_ Type of the matrix value used in the multiplication.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricMultiplyScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::MULTIPLY, true, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar division helper.
 *
 * @tparam right_ Whether the scalar should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the division.
 * @tparam InputValue_ Type of the matrix value used in the division.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricDivideScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::DIVIDE, right_, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar power helper.
 *
 * @tparam right_ Whether the scalar should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the power operation.
 * @tparam InputValue_ Type of the matrix value used in the power operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricPowerScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::POWER, right_, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar modulo helper.
 *
 * @tparam right_ Whether the scalar should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the modulo operation.
 * @tparam InputValue_ Type of the matrix value used in the modulo operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricModuloScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::MODULO, right_, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * Convenient alias for the scalar modulo helper.
 *
 * @tparam right_ Whether the scalar should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the modulo operation.
 * @tparam InputValue_ Type of the matrix value used in the modulo operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Scalar_ Type of the scalar value.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Scalar_>
using DelayedUnaryIsometricIntegerDivideScalarHelper = DelayedUnaryIsometricArithmeticScalarHelper<ArithmeticOperation::INTEGER_DIVIDE, right_, OutputValue_, InputValue_, Index_, Scalar_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricAddScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricAddScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubtractScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricSubtractScalarHelper<right_, OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricMultiplyScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricMultiplyScalarHelper<OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricDivideScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricDivideScalarHelper<right_, OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricModuloScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricModuloScalarHelper<right_, OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricPowerScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricPowerScalarHelper<right_, OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Scalar_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricIntegerDivideScalar(Scalar_ scalar) {
    return std::make_shared<DelayedUnaryIsometricIntegerDivideScalarHelper<right_, OutputValue_, InputValue_, Index_, Scalar_> >(std::move(scalar));
}
/**
 * @endcond
 */

/**
 * @brief Helper for delayed unary isometric vector arithmetic.
 *
 * This class applies the specified arithmetic operation to each element of a `Matrix` where the other operand is row/column-specific value.
 * It should be used as the `Operation_` in the `DelayedUnaryIsometricOperation` class.
 * It may be used regardless of whether `InputValue_` and `OutputValue_` are equal (or not).
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the vector's values should be on the right hand side of the arithmetic operation.
 * Ignored for commutative operations, e.g., `op_ == ADD` and `MULTIPLY`.
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 * This should have a `[]` accessor method and a `size()` method.
 */
template<ArithmeticOperation op_, bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
class DelayedUnaryIsometricArithmeticVectorHelper final : public DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> {
public:
    /**
     * @param vector Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `by_row == true`, otherwise it should be of length equal to the number of columns.
     * @param by_row Whether `vector` corresponds to the rows.
     * If true, each element of the vector is assumed to correspond to a row, and that element is used as an operand with all entries in the same row of the matrix.
     * If false, each element of the vector is assumed to correspond to a column instead.
     */
    DelayedUnaryIsometricArithmeticVectorHelper(Vector_ vector, const bool by_row) : my_vector(std::move(vector)), my_by_row(by_row) {
        for (const auto x : my_vector) {
            if (!delayed_arithmetic_actual_sparse<op_, right_, OutputValue_, InputValue_>(x)) {
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
    void dense(const bool row, const Index_ idx, const Index_ start, const Index_ length, const InputValue_* input, OutputValue_* const output) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input, length, my_vector[idx], output);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to skip aliasing protection.
            }
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_arithmetic<op_, right_>(input[i], my_vector[i + start]);
            }
        }
    }

    void dense(const bool row, const Index_ idx, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* const output) const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input, static_cast<Index_>(indices.size()), my_vector[idx], output);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input = output; // basically an assertion to the compiler to skip aliasing protection.
            }
            Index_ length = indices.size();
            for (Index_ i = 0; i < length; ++i) {
                output[i] = delayed_arithmetic<op_, right_>(input[i], my_vector[indices[i]]);
            }
        }
    }

public:
    bool is_sparse() const {
        return my_sparse;
    }

    void sparse(
        const bool row,
        const Index_ idx,
        const Index_ number,
        const InputValue_* input_value,
        const Index_* const index,
        OutputValue_* const output_value)
    const {
        if (row == my_by_row) {
            delayed_arithmetic_run_simple<op_, right_>(input_value, number, my_vector[idx], output_value);
        } else {
            if constexpr(std::is_same<InputValue_, OutputValue_>::value) {
                input_value = output_value; // basically an assertion to the compiler to skip aliasing protection.
            }
            for (Index_ i = 0; i < number; ++i) {
                output_value[i] = delayed_arithmetic<op_, right_>(input_value[i], my_vector[index[i]]);
            }
        }
    }

    OutputValue_ fill(const bool row, const Index_ idx) const {
        if (row == my_by_row) {
            return delayed_arithmetic_zero<op_, right_, OutputValue_, InputValue_>(my_vector[idx]);
        } else {
            // We should only get to this point if it's sparse, otherwise no
            // single fill value would work across the length of my_vector.
            return 0;
        }
    }
};

/**
 * Convenient alias for the vector addition helper.
 *
 * @tparam OutputValue_ Type of the result of the addition.
 * @tparam InputValue_ Type of the matrix value used in the addition.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector. 
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricAddVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::ADD, true, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector subtraction helper.
 *
 * @tparam right_ Whether the vector should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the subtraction.
 * @tparam InputValue_ Type of the matrix value used in the subtraction.
 * @tparam Index_ Type of index.
 * @tparam Vector_ Type of the vector.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricSubtractVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::SUBTRACT, right_, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector multiplication helper.
 *
 * @tparam OutputValue_ Type of the result of the multiplication.
 * @tparam InputValue_ Type of the matrix value used in the multiplication.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricMultiplyVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::MULTIPLY, true, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector division helper.
 *
 * @tparam right_ Whether the vector should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the division.
 * @tparam InputValue_ Type of the matrix value used in the division.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricDivideVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::DIVIDE, right_, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector power helper.
 *
 * @tparam right_ Whether the vector should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the power operation.
 * @tparam InputValue_ Type of the matrix value used in the power operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricPowerVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::POWER, right_, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector modulo helper.
 *
 * @tparam right_ Whether the vector should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the modulo operation.
 * @tparam InputValue_ Type of the matrix value used in the modulo operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricModuloVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::MODULO, right_, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * Convenient alias for the vector modulo helper.
 *
 * @tparam right_ Whether the vector should be on the right hand side.
 * @tparam OutputValue_ Type of the result of the modulo operation.
 * @tparam InputValue_ Type of the matrix value used in the modulo operation.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Vector_ Type of the vector.
 */
template<bool right_, typename OutputValue_, typename InputValue_, typename Index_, typename Vector_>
using DelayedUnaryIsometricIntegerDivideVectorHelper = DelayedUnaryIsometricArithmeticVectorHelper<ArithmeticOperation::INTEGER_DIVIDE, right_, OutputValue_, InputValue_, Index_, Vector_>;

/**
 * @cond
 */
// Back-compatibility only.
template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricAddVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricAddVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricSubtractVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricSubtractVectorHelper<right_, OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricMultiplyVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricMultiplyVectorHelper<OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricDivideVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricDivideVectorHelper<right_, OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricModuloVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricModuloVectorHelper<right_, OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricPowerVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricPowerVectorHelper<right_, OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}

template<bool right_, typename OutputValue_ = double, typename InputValue_ = double, typename Index_ = int, typename Vector_>
std::shared_ptr<DelayedUnaryIsometricOperationHelper<OutputValue_, InputValue_, Index_> > make_DelayedUnaryIsometricIntegerDivideVector(Vector_ vector, bool by_row) {
    return std::make_shared<DelayedUnaryIsometricIntegerDivideVectorHelper<right_, OutputValue_, InputValue_, Index_, Vector_> >(std::move(vector), by_row);
}
/**
 * @endcond
 */

}

#endif
