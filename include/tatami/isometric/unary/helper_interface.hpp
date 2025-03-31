#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OPERATION_HELPER_INTERFACE_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OPERATION_HELPER_INTERFACE_H

#include <vector>
#include <optional>
#include "../../base/SparseRange.hpp"

/** 
 * @file helper_interface.hpp
 * @brief Interface for `tatami::DelayedUnaryIsometricOperation` helpers.
 */

namespace tatami {

/**
 * @brief Helper operation interface for `DelayedUnaryIsometricOperation`.
 * @brief Advanced mock operation for `DelayedUnaryIsometricOperation`.
 *
 * This class defines the interface for an operation helper in `DelayedUnaryIsometricOperation`,
 * Operations should generally inherit from this class, though it is possible for developers to define their own classes with the same signatures for compile-time polymorphism.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedUnaryIsometricOperationHelper {
public:
    /**
     * @cond
     */
    DelayedUnaryIsometricOperationHelper() = default;
    DelayedUnaryIsometricOperationHelper(const DelayedUnaryIsometricOperationHelper&) = default;
    DelayedUnaryIsometricOperationHelper& operator=(const DelayedUnaryIsometricOperationHelper&) = default;
    DelayedUnaryIsometricOperationHelper(DelayedUnaryIsometricOperationHelper&&) = default;
    DelayedUnaryIsometricOperationHelper& operator=(DelayedUnaryIsometricOperationHelper&&) = default;
    virtual ~DelayedUnaryIsometricOperationHelper() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method accepts a contiguous block of an element of the target dimension from the underlying matrix (`input`),
     * applies the operation to each value, and stores the result in another array of different type (`output`).
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
    virtual void dense(bool row, Index_ i, Index_ start, Index_ length, const InputValue_* input, OutputValue_* output) const = 0;

    /**
     * This method accepts an indexed subset of an element of the target dimension from the underlying matrix (`input`),
     * applies the operation to each value, and stores the result in another array of different type (`output`).
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
    virtual void dense(bool row, Index_ i, const std::vector<Index_>& indices, const InputValue_* input, OutputValue_* output) const = 0;

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
    virtual void sparse(bool row, Index_ i, Index_ num, const InputValue_* input_value, const Index_* index, OutputValue_* output_value) const = 0;

    /**
     * @param row Whether `i` refers to the row or column index.
     * @param i The index of the row (if `row = true`) or column (otherwise) containing the zeros.
     * This argument should be ignored if the operation does not depend on the row/column,
     * i.e., when all of `zero_depends_on_row()` and friends return false.
     *
     * @return The result of the operation being applied on zeros from the `i`-th row/column of the matrix.
     *
     * This function will never be called by `DelayedUnaryIsometricOperation` if the operation depends on the dimension that is not specified by `row`,
     * i.e., when `row = true && zero_depends_on_column()` or `row = false && zero_depends_on_row()`.
     * In such cases, no single fill value would exist.
     */
    virtual OutputValue_ fill(bool row, Index_ i) const = 0; 

    /**
     * @return Whether the operation will convert a structural zero to a non-zero value,
     * in a manner that depends on the identity of the column in which the structural zero occurs.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedUnaryIsometricOperation` will automatically recognize such operations as being row-independent.
     */
    virtual bool zero_depends_on_row() const = 0;

    /**
     * @return Whether the operation will convert a structural zero to a non-zero value,
     * in a manner that depends on the identity of the column in which the structural zero occurs.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedUnaryIsometricOperation` will automatically recognize such operations as being row-independent.
     */
    virtual bool zero_depends_on_column() const = 0;

    /**
     * @return Whether the result of the operation on a non-zero operand depends on the identity of the row containing the operand.
     */
    virtual bool non_zero_depends_on_row() const = 0;

    /**
     * @return Whether the result of the operation on a non-zero operand depends on the identity of the column containing the operand.
     */
    virtual bool non_zero_depends_on_column() const = 0;

    /** 
     * @return Does this operation preserve sparsity?
     * This may return false.
     */
    virtual bool is_sparse() const = 0;

    /**
     * @return Expected number of rows in the matrix to which this operation is to be applied (i.e., the underlying matrix in the `DelayedUnaryIsometricOperation` constructor).
     * If no value is returned, the matrix may have any number of rows.
     */
    virtual std::optional<Index_> nrow() const = 0; 

    /**
     * @return Expected number of columns in the matrix to which this operation is to be applied (i.e., the underlying matrix in the `DelayedUnaryIsometricOperation` constructor).
     * If no value is returned, the matrix may have any number of columns.
     */
    virtual std::optional<Index_> ncol() const = 0; 
};

}

#endif
