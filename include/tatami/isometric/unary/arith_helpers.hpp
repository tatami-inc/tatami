#ifndef TATAMI_ARITH_HELPERS_H
#define TATAMI_ARITH_HELPERS_H

/**
 * @file arith_helpers.hpp
 *
 * @brief Helper classes for scalar arithmetic operations.
 * 
 * Classes defined here should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 */

namespace tatami {

/**
 * Type of the delayed arithmetic operation.
 */
enum class DelayedArithOp : char { 
    ADD, 
    SUBTRACT,
    MULTIPLY,
    DIVIDE
};

/**
 * @cond
 */
template<DelayedArithOp op_, typename Scalar_, typename Value_, typename Index_>
void delayed_arith_run_simple(Scalar_ scalar, Index_ length, Value_* buffer) {
    for (Index_ i = 0; i < length; ++i) {
        if constexpr(op_ == DelayedArithOp::ADD) {
            buffer[i] += scalar;
        } else if constexpr(op_ == DelayedArithOp::MULTIPLY) {
            buffer[i] *= scalar;
        } else if constexpr(op_ == DelayedArithOp::SUBTRACT) {
            if constexpr(right_) {
                buffer[i] -= scalar;
            } else {
                buffer[i] = scalar - buffer[i];
            }
        } else {
            // Assume IEEE behavior if divisor is zero.
            if constexpr(right_) {
                buffer[i] /= scalar;
            } else {
                buffer[i] = scalar / buffer[i];
            }
        }
    }
}
/**
 * @endcond
 */

/**
 * @brief Delayed scalar arithmetic.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the scalar should be on the right hand side of the arithmetic operation.
 * Ignored for some `op_`.
 * @tparam Scalar_ Type of the scalar value.
 */
template<DelayedArithOp op_, bool right_, typename Scalar_>
struct DelayedArithScalarHelper {
    /**
     * @param s Scalar value to be added.
     */
    DelayedArithScalarHelper(Scalar_ s) : scalar(s) {}

private:
    const Scalar_ scalar;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool potential_sparse = true;

    bool actual_sparse() const {
        if constexpr(op_ == DelayedArithOp::ADD || op_ == DelayedArithOp::SUBTRACT) {
            return scalar == 0;
        } else if (op_ == DelayedArithOp::MULTIPLY) {
            return true;
        } else {
            if constexpr(right_) {
                return scalar != 0;
            } else {
                return false;
            }
        }
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void dense(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_arith_run_simple(scalar, length, buffer);
    }

    template<bool, typename Value_, typename Index_>
    void sparse(Index_, Index_ number, Value_* buffer, const Index_*) const {
        delayed_arith_run_simple(scalar, range.number, buffer);
    }

    template<bool, typename Value_, typename Index_, typename ExtractType_>
    void expanded(Index_, ExtractType_, Index_ length, Value_* buffer) const {
        delayed_arith_run_simple(scalar, length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Delayed vector arithmetic.
 *
 * This should be used as the `OP` in the `DelayedUnaryIsometricOp` class.
 *
 * @tparam op_ The arithmetic operation.
 * @tparam right_ Whether the vector's values should be on the right hand side of the arithmetic operation.
 * Ignored for some `op_`.
 * @tparam margin_ Matrix dimension along which the operation is to occur.
 * If 0, each element of the vector is assumed to correspond to a row, and that value is subtracted from all entries in the same row of the matrix.
 * If 1, each element of the vector is assumed to correspond to a column instead.
 * @tparam Vector_ Type of the vector.
 */
template<DelayedArithOp op_, bool right_, int margin_, typename Vector_>
struct DelayedArithVectorHelper {
    /**
     * @param v Vector of values to use in the operation. 
     * This should be of length equal to the number of rows if `MARGIN = 0`, otherwise it should be of length equal to the number of columns.
     */
    DelayedSubtractVectorHelper(Vector_ v) : vec(std::move(v)) {
        if constexpr(op_ == DelayedArithOp::ADD || op_ == DelayedArithOp::SUBTRACT) {
            for (auto x : vec) {
                if (x != 0) {
                    still_sparse = false;
                    break;
                }
            }
        } else if constexpr(op_ == DelayedArithOp::DIVIDE) {
            if constexpr(right_) {
                // Division by any zero value renders it non-sparse in my eyes.
                for (auto x : vec) {
                    if (x == 0) {
                        still_sparse = false;
                        break;
                    }
                }
            }
        }
    }

private:
    const Vector_ vec;
    bool still_sparse = true;

public:
    /**
     * @cond
     */
    static constexpr bool needs_row = false;

    static constexpr bool needs_column = false;

    static constexpr bool potential_sparse = true;

    bool actual_sparse() const {
        if constexpr(op_ == DelayedArithOp::ADD || op_ == DelayedArithOp::SUBTRACT) {
            return still_sparse;
        } else if (op_ == DelayedArithOp::MULTIPLY) {
            return true;
        } else {
            if constexpr(right_) {
                return still_sparse;
            } else {
                return false;
            }
        }
    }
    /**
     * @endcond
     */

public:
    /**
     * @cond
     */
    template<bool accrow_, typename Value_, typename Index_, typename ExtractType>
    void dense(Index_ idx, ExtractType_ start, Index_ length, Value_* buffer) {
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_arith_run_simple(vec[idx], length, buffer);

        } else if constexpr(std::is_same<ExtractType_, Index_>::value) {
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(op_ == DelayedArithOp::ADD) {
                    buffer[i] += vec[i + start];
                } else if constexpr(op_ == DelayedArithOp::MULTIPLY) {
                    buffer[i] *= vec[i + start];
                } else if constexpr(op_ == DelayedArithOp::SUBTRACT) {
                    if constexpr(right_) {
                        buffer[i] -= vec[i + start];
                    } else {
                        buffer[i] = vec[i + start] - buffer[i];
                    }
                } else {
                    // Assume IEEE behavior if divisor is zero.
                    if constexpr(right_) {
                        buffer[i] /= vec[i + start];
                    } else {
                        buffer[i] = vec[i + start] / buffer[i];
                    }
                }
            }

        } else {
            for (Index_ i = 0; i < length; ++i) {
                auto scalar = vec[start[i]];
                if constexpr(op_ == DelayedArithOp::ADD) {
                    buffer[i] += scalar;
                } else if constexpr(op_ == DelayedArithOp::MULTIPLY) {
                    buffer[i] *= scalar;
                } else if constexpr(op_ == DelayedArithOp::SUBTRACT) {
                    if constexpr(right_) {
                        buffer[i] -= scalar;
                    } else {
                        buffer[i] = scalar - buffer[i];
                    }
                } else {
                    // Assume IEEE behavior if divisor is zero.
                    if constexpr(right_) {
                        buffer[i] /= scalar;
                    } else {
                        buffer[i] = scalar / buffer[i];
                    }
                }
            }
        }
    }

    template<bool, typename Value_, typename Index_>
    void run(Index_ idx, Index_ number, Value_* buffer, const Index_* indices);
        if constexpr(accrow_ == (margin_ == 0)) {
            delayed_arith_run_simple(vec[idx], range.number, buffer);

        } else {
            for (Index_ i = 0; i < length; ++i) {
                if constexpr(op_ == DelayedArithOp::ADD) {
                    values[indices[i]] += vec[indices[i]];
                } else if constexpr(op_ == DelayedArithOp::MULTIPLY) {
                    values[indices[i]] *= vec[indices[i]];
                } else if constexpr(op_ == DelayedArithOp::SUBTRACT) {
                    if constexpr(right_) {
                        values[indices[i]] -= vec[indices[i]];
                    } else {
                        values[indices[i]] = vec[indices[i]] - values[i];
                    }
                } else {
                    // Assume IEEE behavior if divisor is zero.
                    if constexpr(right_) {
                        values[indices[i]] /= vec[indices[i]];
                    } else {
                        values[indices[i]] = vec[indices[i]] / values[i];
                    }
                }
            }
        }
    }

    template<bool accrow_, typename Value_, typename Index_, typename ExtractType>
    void expanded(Index_ idx, ExtractType_&& start, Index_ length, Value_* buffer) {
        dense<accrow_>(idx, std::forward<ExtractType>(start), length, buffer);
    }
    /**
     * @endcond
     */
};

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be added.
 * @return A helper class for delayed scalar addition.
 */
template<typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::ADD, true, Scalar_> make_DelayedAddScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::ADD, true, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be subtracted.
 * @return A helper class for delayed scalar subtraction.
 */
template<bool right_, typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::SUBTRACT, right_, Scalar_> make_DelayedSubtractScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::SUBTRACT, right_, Scalar_>(std::move(s));
}

/**
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be multiplied.
 * @return A helper class for delayed scalar multiplication.
 */
template<typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::MULTIPLY, true, Scalar_> make_DelayedMultiplyScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::MULTIPLY, true, Scalar_>(std::move(s));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam Scalar_ Type of the scalar.
 * @param s Scalar value to be divided.
 * @return A helper class for delayed scalar division.
 */
template<bool right_, typename Scalar_>
DelayedArithScalarHelper<DelayedArithOp::DIVIDE, right_, Scalar_> make_DelayedDivideScalarHelper(Scalar_ s) {
    return DelayedArithScalarHelper<DelayedArithOp::DIVIDE, right_, Scalar_>(std::move(s));
}

/**
 * @tparam margin_ Matrix dimension along which the addition is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to be added to the rows/columns.
 * @return A helper class for delayed vector addition.
 */
template<int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::ADD, true, margin_, Vector_> make_DelayedAddVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::ADD, true, margin_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the subtraction.
 * @tparam margin_ Matrix dimension along which the subtraction is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to subtract from (or be subtracted by) the rows/columns.
 * @return A helper class for delayed vector subtraction.
 */
template<bool right_, int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::SUBTRACT, right_, margin_, Vector_> make_DelayedSubtractVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::SUBTRACT, right_, margin_, Vector_>(std::move(v));
}

/**
 * @tparam margin_ Matrix dimension along which the multiplication is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to multiply the rows/columns.
 * @return A helper class for delayed vector multiplication.
 */
template<int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::MULTIPLY, true, margin_, Vector_> make_DelayedMultiplyVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::MULTIPLY, true, margin_, Vector_>(std::move(v));
}

/**
 * @tparam right_ Whether the scalar should be on the right hand side of the division.
 * @tparam margin_ Matrix dimension along which the division is to occur, see `DelayedArithVectorHelper`.
 * @tparam Vector_ Type of the vector.
 *
 * @param v Vector to divide (or be divided by) the rows/columns.
 * @return A helper class for delayed vector division.
 */
template<bool right_, int margin_, typename Vector_>
DelayedArithVectorHelper<DelayedArithOp::DIVIDE, right_, margin_, Vector_> make_DelayedDivideVectorHelper(Vector_ v) {
    return DelayedArithVectorHelper<DelayedArithOp::DIVIDE, right_, margin_, Vector_>(std::move(v));
}

}

#endif
