#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_HELPER_INTERFACE_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_HELPER_INTERFACE_H

#include <vector>
#include "../../base/SparseRange.hpp"

/** 
 * @file mock_helpers.hpp
 * @brief Expectations for `tatami::DelayedUnaryIsometricOperation` helpers.
 */

namespace tatami {

/**
 * @brief Basic mock operation for a `DelayedUnaryIsometricOperation`. 
 *
 * This defines the basic expectations for an operation to be used as `Operation_` in `DelayedUnaryIsometricOperation`.
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
class DelayedUnaryIsometricMockBasic {
public:
    /**
     * This method accepts a contiguous block of an element of the target dimension from the underlying matrix (`input`),
     * applies the operation to each value, and stores the result in another array of a potentially different type (`output`).
     *
     * Implementions of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to use in the operation.
     * @tparam OutputValue_ Type of matrix value returned by the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `buffer` contains row `i`, otherwise it contains column `i`.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * Unlike `DelayedUnaryIsometricMockAdvanced::dense()`, this is always guaranteed to be available.
     * @param start Start of the contiguous block of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in] input Pointer to an array containing a contiguous block of a row/column extracted from the matrix.
     * This has `length` addressable elements.
     * @param[out] output Pointer to an array to store the results of the operation applied to elements of `input`. 
     * This has `length` addressable elements. 
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `input`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] Index_ start, 
        Index_ length, 
        [[maybe_unused]] const InputValue_* input,
        OutputValue_* output)
    const {
        std::fill_n(output, length, 0);
    }

    /**
     * This method accepts an indexed subset of an element of the target dimension from the underlying matrix (`input`),
     * applies the operation to each value, and stores the result in another array of different type (`output`).
     *
     * Implementions of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam InputValue_ Type of matrix value to use in the operation.
     * @tparam Index_ Type of index value.
     * @tparam OutputValue_ Type of matrix value returned by the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `buffer` contains row `i`, otherwise it contains column `i`.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * Unlike `DelayedUnaryIsometricMockAdvanced::dense()`, this is always guaranteed to be available.
     * @param indices Sorted and unique indices of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param[in] input Pointer to an array containing an indexed subset of a row/column extracted from the matrix.
     * This has `length` addressable elements.
     * @param[out] output Pointer to an array to store the results of the operation applied to elements of `input`. 
     * This has `length` addressable elements. 
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `input`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        const std::vector<Index_>& indices, 
        [[maybe_unused]] const InputValue_* input,
        OutputValue_* output)
    const {
        std::fill_n(output, indices.size(), 0);
    }

    /**
     * Whether this is a basic operation.
     * This should be true, otherwise an advanced operation interface is expected (see `DelayedUnaryIsometricMockAdvanced`).
     */
    static constexpr bool is_basic = true;
};

/**
 * @brief Advanced mock operation for `DelayedUnaryIsometricOperation`.
 *
 * This class defines the advanced expectations for an operation in `DelayedUnaryIsometricOperation`,
 * which improves efficiency by taking advantage of any sparsity in the underlying matrix.
 * Either the operation itself preserves sparsity, or any loss of sparsity is predictable,
 * i.e., zeros are transformed into a constant non-zero value that does not depend on its position in the `Matrix`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
class DelayedUnaryIsometricMockAdvanced {
public:
    /**
     * This method accepts a contiguous block of an element of the target dimension from the underlying matrix (`input`),
     * applies the operation to each value, and stores the result in another array of different type (`output`).
     *
     * Implementions of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to use in the operation.
     * @tparam OutputValue_ Type of matrix value returned by the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `buffer` contains row `i`, otherwise it contains column `i`.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     * @param start Start of the contiguous block of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in] input Pointer to an array containing a contiguous block of a row/column extracted from the matrix.
     * This has `length` addressable elements.
     * @param[out] output Pointer to an array to store the results of the operation applied to elements of `input`. 
     * This has `length` addressable elements. 
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `input`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] Index_ start, 
        Index_ length, 
        [[maybe_unused]] const InputValue_* input,
        OutputValue_* output)
    const {
        std::fill_n(output, length, 0);
    }

    /**
     * This method accepts an indexed subset of an element of the target dimension from the underlying matrix (`input`),
     * applies the operation to each value, and stores the result in another array of different type (`output`).
     *
     * Implementions of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam InputValue_ Type of matrix value to use in the operation.
     * @tparam Index_ Type of index value.
     * @tparam OutputValue_ Type of matrix value returned by the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `buffer` contains row `i`, otherwise it contains column `i`.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     * @param indices Sorted and unique indices of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param[in] input Pointer to an array containing an indexed subset of a row/column extracted from the matrix.
     * This has `length` addressable elements.
     * @param[out] output Pointer to an array to store the results of the operation applied to elements of `input`. 
     * This has `length` addressable elements. 
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `input`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        const std::vector<Index_>& indices, 
        [[maybe_unused]] const InputValue_* input,
        OutputValue_* output)
    const {
        std::fill_n(output, indices.size(), 0);
    }

    /**
     * This method is expected to iterate over `input_value`, apply the operation to each value, and store the result in `output_value`.
     * We assume that the operation only needs to be applied to the structural non-zeros;
     * structural zeros are either ignored for sparsity-preserving operations,
     * or the result of the operation on zeros will be populated by `fill()`.
     *
     * If `non_zero_depends_on_row() && !row` or `non_zero_depends_on_column() && row`, `index` is guaranteed to be non-NULL.
     * Otherwise, it may be NULL and should be ignored.
     * Even if non-NULL, indices are not guaranteed to be sorted.
     *
     * Implementations of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to use in the operation.
     * @tparam OutputValue_ Type of matrix value returned by the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `buffer` contains row `i`, otherwise it contains column `i`.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     * @param num Number of non-zero elements for row/column `i`.
     * @param[in] input_value Pointer to an array of values of the structural non-zero elements from the row/column of the matrix.
     * This is guaranteed to have `num` addressable elements.
     * @param[in] index Pointer to an array of column (if `row = true`) or row indices (otherwise) of the non-zero elements.
     * Alternatively NULL.
     * @param[out] output_value Pointer to an array in which to store the result of the operation on each element of `input_value`.
     * This is guaranteed to have `num` addressable elements.
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `input`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void sparse(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        Index_ num,
        [[maybe_unused]] const InputValue_* input_value,
        [[maybe_unused]] const Index_* index,
        OutputValue_* output_value)
    const {
        std::fill_n(output_value, num, 0);
    }

    /**
     * @tparam OutputValue_ Type of the result of the operation.
     * @tparam InputValue_ Type of the matrix value used in the operation.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `i` refers to the row or column index.
     * @param i The index of the row (if `row = true`) or column (otherwise) containing the zeros.
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     *
     * @return The result of the operation being applied on zeros from the `i`-th row/column of the matrix.
     *
     * This method will be called with the explicit `OutputValue_` and `InputValue_` template parameters.
     * Implementations of this method should either ensure that `Index_` is deducible or use a fixed integer type in the method signature.
     */
    template<typename OutputValue_, typename InputValue_, typename Index_>
    OutputValue_ fill([[maybe_unused]] bool row, [[maybe_unused]] Index_ i) const { 
        return 0;
    }

    /**
     * Whether this is a basic operation.
     * This should be false, otherwise a basic operation interface is expected (see `DelayedUnaryIsometricMockBasic`).
     */
    static constexpr bool is_basic = false;

    /**
     * @return Whether the operation will convert a structural zero to a non-zero value,
     * in a manner that depends on the identity of the column in which the structural zero occurs.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedUnaryIsometricOperation` will automatically recognize such operations as being row-independent.
     * 
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool zero_depends_on_row() const {
        return false;
    }

    /**
     * @return Whether the operation will convert a structural zero to a non-zero value,
     * in a manner that depends on the identity of the column in which the structural zero occurs.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedUnaryIsometricOperation` will automatically recognize such operations as being row-independent.
     *
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool zero_depends_on_column() const {
        return false;
    }

    /**
     * @return Whether the result of the operation on a non-zero operand depends on the identity of the row containing the operand.
     *
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool non_zero_depends_on_row() const {
        return false;
    }

    /**
     * @return Whether the result of the operation on a non-zero operand depends on the identity of the column containing the operand.
     *
     * This method may also omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool non_zero_depends_on_column() const {
        return false;
    }

    /** 
     * @return Does this operation preserve sparsity?
     * This may return false.
     */
    bool is_sparse() const {
        return true;
    }
};

}

#endif
