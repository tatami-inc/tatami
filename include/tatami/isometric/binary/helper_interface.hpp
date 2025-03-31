#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OPERATION_HELPER_INTERFACE_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OPERATION_HELPER_INTERFACE_H

#include <vector>
#include <optional>
#include "../../base/SparseRange.hpp"

/** 
 * @file helper_interface.hpp
 * @brief Interface for `tatami::DelayedBinaryIsometricOperation` helpers.
 */

namespace tatami {

/**
 * @brief Helper operation interface for `DelayedBinaryIsometricOperation`.
 *
 * This class defines the interface for an operation helper in `DelayedBinaryIsometricOperation`.
 * Operations should generally inherit from this class, though it is possible for developers to define their own classes with the same signatures for compile-time polymorphism.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the matrix value used in the operation.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename OutputValue_, typename InputValue_, typename Index_>
class DelayedBinaryIsometricOperationHelper {
public:
    /**
     * @cond
     */
    DelayedBinaryIsometricOperationHelper() = default;
    DelayedBinaryIsometricOperationHelper(const DelayedBinaryIsometricOperationHelper&) = default;
    DelayedBinaryIsometricOperationHelper& operator=(const DelayedBinaryIsometricOperationHelper&) = default;
    DelayedBinaryIsometricOperationHelper(DelayedBinaryIsometricOperationHelper&&) = default;
    DelayedBinaryIsometricOperationHelper& operator=(DelayedBinaryIsometricOperationHelper&&) = default;
    virtual ~DelayedBinaryIsometricOperationHelper() = default;
    /**
     * @endcond
     */

public:
    /**
     * @param row Whether `i` refers to the row or column index.
     * @param i The index of the row (if `row = true`) or column (otherwise) containing the zeros.
     * This argument should be ignored if the operation does not depend on the row/column, 
     * i.e., when `row = true && !zero_depends_on_row()` or `row = false && !zero_depends_on_column()`.
     *
     * @return The result of `OP(lz, rz)` where `OP` is the operation,
     * `lz` is a structural zero from the `i`-th row/column of the left matrix,
     * and `rz` is a structural zero from the `i`-th row/column of the left matrix.
     *
     * This function will never be called by `DelayedBinaryIsometricOperation` if the operation depends on the dimension that is not specified by `row`,
     * i.e., when `row = true && zero_depends_on_column()` or `row = false && zero_depends_on_row()`.
     * In such cases, no single fill value would exist.
     */
    virtual OutputValue_ fill(bool row, Index_ i) const = 0;

    /**
     * @return Whether applying the operation to a pair of structural zeros (one from each matrix)
     * yields a value that depends on the identity of the row containing those zeros.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedBinaryIsometricOperation` will automatically recognize such operations as being row-independent.
     */
    virtual bool zero_depends_on_row() const = 0;

    /**
     * @return Whether applying the operation to a pair of structural zeros (one from each matrix)
     * yields a value that depends on the identity of the column containing those zeros.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedBinaryIsometricOperation` will automatically recognize such operations as being row-independent.
     */
    virtual bool zero_depends_on_column() const = 0;

    /**
     * @return Whether the result of the operation depends on the identity of the row containing the operands,
     * where at least one of the operands is non-zero.
     */
    virtual bool non_zero_depends_on_row() const = 0;

    /**
     * @return Whether the result of the operation depends on the identity of the column containing the operands,
     * where at least one of the operands is non-zero.
     */
    virtual bool non_zero_depends_on_column() const = 0;

    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
     * These buffers represent the same element of the target dimension from the left and right matrices, respectively, in dense form.
     * Each buffer holds values from a contiguous block of the non-target dimension.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `left_buffer` and `right_buffer` hold the contents of the `i`-th row from both matrices;
     * otherwise, they hold the contents of the `i`-th column.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     * @param start Start of the contiguous block of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in,out] left_buffer Pointer to an array containing the row/column extracted from the left matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     * @param[in] right_buffer Pointer to an array containing the row/column extracted from the right matrix.
     * This has `length` addressable elements.
     * @param[out] output_buffer Pointer to an array in which to store the result of the operation.
     * This has `length` addressable elements.
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `left_buffer`.
     */
    virtual void dense(bool row, Index_ i, Index_ start, Index_ length, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const = 0;

    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
     * These buffers represent the same element of the target dimension from the left and right matrices, respectively, in dense form.
     * Each buffer holds values from an indexed subset of the non-target dimension.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `left_buffer` and `right_buffer` hold the contents of the `i`-th row from both matrices;
     * otherwise, they hold the contents of the `i`-th column.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     * @param indices Sorted and unique indices of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param[in,out] left_buffer Pointer to an array containing the row/column extracted from the left matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     * @param[in] right_buffer Pointer to an array containing the row/column extracted from the right matrix.
     * This has `length` addressable elements.
     * @param[out] output_buffer Pointer to an array in which to store the result of the operation.
     * This has `length` addressable elements.
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `left_buffer`.
     */
    virtual void dense(bool row, Index_ i, const std::vector<Index_>& indices, const InputValue_* left_buffer, const InputValue_* right_buffer, OutputValue_* output_buffer) const = 0;

    /**
     * This method applies the operation to the sparse ranges in `left` and `right`, 
     * containing values from the same element of the target dimension from the left and right matrices, respectively.
     * Specifically, the operation only needs to be applied to the structural non-zeros,
     * and results of the operation should be stored in the `output_*` buffers.
     * Structural zeros are either ignored for sparsity-preserving operations,
     * or the result of the operation on zeros will be populated by `fill()`.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `left_buffer` and `right_buffer` hold the contents of the `i`-th row from both matrices;
     * otherwise, they hold the contents of the `i`-th column.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     * @param left Contents of row/column `i` extracted from the left matrix.
     * @param right Contents of row/column `i` extracted from the right matrix.
     * @param[out] output_value Pointer to an array for storing output values of the operation.
     * This is guaranteed to have enough space for the union of indices in `left` and `right`.
     * @param[out] output_index Pointer to an array for storing output indices.
     * This is guaranteed to have enough space for the union of indices in `left` and `right`.
     * @param report_value Whether to return the values in `output_value`.
     * @param report_index Whether to return the indices in `output_index`.
     *
     * @return Number of structural non-zero elements in the `output_*` buffers.
     *
     * If `report_value = true`, `left.value` and `right.value` and `output_value` are all guaranteed to be non-NULL.
     * Otherwise, any of these pointers may be NULL and should be ignored.
     *
     * If `report_index = true`, `output_index` is guaranteed to be non-NULL; otherwise, `output_index` should be ignored.
     * `left.index` and `right.index` will always return be non-NULL regardless of `report_index`.
     * Indices in `left.index` and `right.index` are also guaranteed to be in ascending order.
     *
     * It is expected that the results of the operation are sorted in ascending order,
     * i.e., indices in `output_index` should be increasing.
     *
     * The settings of `report_index` and `report_value` should not change the number or ordering of the results.
     * That is, `output_index` should have the same indices regardless of `report_value`,
     * and `output_value` should have the same values regardless of `report_index`.
     * This implies that implementations should not omit structural non-zeros even if the actual value is zero,
     * as the computation of the actual value requires `report_value = true`.
     */
    virtual Index_ sparse(
        bool row, 
        Index_ i, 
        const SparseRange<InputValue_, Index_>& left, 
        const SparseRange<InputValue_, Index_>& right, 
        OutputValue_* output_value, 
        Index_* output_index, 
        bool report_value, 
        bool report_index)
    const = 0;

    /** 
     * @return Whether this operation preserves sparsity.
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
