#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OP_HELPER_INTERFACE_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OP_HELPER_INTERFACE_H

#include <vector>
#include "../../base/SparseRange.hpp"

namespace tatami {

/**
 * @brief Mock helper for dense operations in `DelayedBinaryIsometricOp`.
 *
 * This defines the expectations for operations that discard sparsity in a variable manner
 * (i.e., zeros are transformed into non-zeros of different value depending on their position in the `Matrix`).
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryMockVariableDenseHelper {
    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
     * representing the dense contents of the same row/column from the left and right matrices respectively.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `left_buffer` and `right_buffer` contain the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param start Start of the contiguous block of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in,out] left_buffer Contents of the row/column extracted from the left matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     * @param[int] right_buffer Contents of the row/column extracted from the right matrix.
     * This has `length` addressable elements.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] Index_ start, 
        Index_ length, 
        Value_* left_buffer, 
        [[maybe_unused]] const Value_* right_buffer)
    const {
        // Just filling it with something as a mock.
        std::fill_n(left_buffer, length, 0);
    }

    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
     * representing the dense contents of the same row/column from the left and right matrices respectively.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `left_buffer` and `right_buffer` contain the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param indices Sorted and unique indices of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in,out] left_buffer Contents of the row/column extracted from the left matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     * @param[int] right_buffer Contents of the row/column extracted from the right matrix.
     * This has `length` addressable elements.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        const std::vector<Index_>& indices, 
        Value_* left_buffer, 
        [[maybe_unused]] const Value_* right_buffer) 
    const {
        std::fill_n(left_buffer, indices.size(), 0);
    }
 
    /**
     * This operation does not preserve sparsity.
     */
    static constexpr bool is_sparse = false;

    /**
     * Conversion of zeros to non-zero values is dependent on rows.
     */
    static constexpr bool zero_depends_on_row = true;

    /**
     * Conversion of zeros to non-zero values is dependent on columns.
     */
    static constexpr bool zero_depends_on_column = true;

    /**
     * @cond
     */
    virtual ~DelayedBinaryMockVariableDenseHelper() = default;
    /**
     * @endcond
     */
};

/**
 * @brief Mock helper for constant dense operations in `DelayedBinaryIsometricOp`.
 *
 * This defines the expectations for operations that discard sparsity in a variable manner
 * (i.e., zeros are transformed into non-zeros of different value depending on their position in the `Matrix`).
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryMockConstantDenseHelper : public DelayedBinaryMockVariableDenseHelper {
    /**
     * This method applies the operation to the sparse ranges in `left` and `right`,
     * representing the contents of the same row/column from the left and right matrices respectively.
     * Specifically, the operation only needs to be applied to the structural non-zeros, as the zeros are populated by `fill()`.
     * The results of the operation should then be stored in the `output_*` buffers.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `left_buffer` and `right_buffer` contain the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param left Contents of row/column `i` extracted from the left matrix.
     * @param right Contents of row/column `i` extracted from the right matrix.
     * @param[out] output_value Pointer to an array for storing output values of the operation.
     * This is guaranteed to have enough space for the union of indices in `left` and `right`.
     * @param[out] output_index Pointer to an array for storing output indices of the operation.
     * This is guaranteed to have enough space for the union of indices in `left` and `right`.
     *
     * @return Number of structural non-zero elements in the `output_*` buffers.
     *
     * Both `left` and `right` are guaranteed to have non-NULL `value` pointers and `index` pointers.
     * Indices in `index` are also guaranteed to be in ascending order.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    Index_ sparse(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] const SparseRange<Value_, Index_>& left, 
        [[maybe_unused]] const SparseRange<Value_, Index_>& right, 
        [[maybe_unused]] Value_* output_value,
        [[maybe_unused]] Index_* output_index)
    const {
        return 0;
    }

    /**
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param i The index of the row containing the zero, if `zero_variance` is `DelayedBinaryDenseZeroVariance::ROW`;
     * the index of the column containing the zero, if `zero_variance` is `DelayedBinaryDenseZeroVariance::COLUMN`;
     * or ignored, if `zero_variance` is `DelayedBinaryDenseZeroVariance::NONE`.
     *
     * @return The result of the operation being applied on zeros from both the left and right matrices.
     * This should be constant for all elements in the row/column/matrix, depending on the interpretation of `i`.
     *
     * This method will be called with the `Value_` template parameter.
     */
    template<typename Value_, typename Index_>
    Value_ fill([[maybe_unused]] Index_ i) const { 
        return 0;
    }

    /**
     * This operation does not preserve sparsity.
     */
    static constexpr bool is_sparse = false;

    /**
     * Conversion of zeros to non-zero values is not dependent on rows.
     * This may also be `true` provided that `zero_depends_on_column = false`.
     */
    static constexpr bool zero_depends_on_row = false;

    /**
     * Conversion of zeros to non-zero values is not dependent on columns.
     * This may also be `true` provided that `zero_depends_on_row = false`.
     */
    static constexpr bool zero_depends_on_column = false;
};

/**
 * @brief Interface for sparse operations in `DelayedBinaryIsometricOp`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryMockSparseHelper : public DelayedBinaryMockVariableDenseHelper {
    /**
     * This method applies the operation to the sparse ranges in `left` and `right`, 
     * representing the contents of the same row/column from the left and right matrices respectively.
     * It should only be called for sparsity-preserving operations, enabling optimizations by skipping structural zeros.
     * The results of the operation should then be stored in the `output_*` buffers.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `left_buffer` and `right_buffer` contain the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param left Contents of row/column `i` extracted from the left matrix.
     * @param right Contents of row/column `i` extracted from the right matrix.
     * @param[out] output_value Pointer to an array for storing output values of the operation.
     * This is guaranteed to have enough space for the union of indices in `left` and `right`.
     * @param[out] output_index Pointer to an array for storing output indices.
     * This is guaranteed to have enough space for the union of indices in `left` and `right`.
     * @param needs_value Whether to return the values in `output_value`.
     * @param needs_index Whether to return the indices in `output_index`.
     *
     * @return Number of structural non-zero elements in the `output_*` buffers.
     *
     * If `needs_value = true`, both `left` and `right` will have non-NULL `value` pointers, and `output_value` is guaranteed to be non-NULL.
     * Otherwise, `left` and `right` will have NULL `value` pointers, and `output_value` may be NULL and should be ignored.
     *
     * If `needs_index = true`, `output_index` is guaranteed to be non-NULL; otherwise, `output_index` should be ignored.
     * `left` and `right` will always return contain non-NULL `index` pointers regardless of `needs_index`, and these are always in ascending order.
     *
     * It is expected that the contents of `output_index` are also sorted in ascending order.
     *
     * As described in `Options`, the settings of `needs_index` and `needs_value` should not change the number or order of reported structural non-zero values.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    Index_ sparse(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] const SparseRange<Value_, Index_>& left, 
        [[maybe_unused]] const SparseRange<Value_, Index_>& right, 
        [[maybe_unused]] Value_* output_value, 
        [[maybe_unused]] Index_* output_index, 
        [[maybe_unused]] bool needs_value, 
        [[maybe_unused]] bool needs_index)
    const {
        return 0;
    }

    /**
     * This operation preserves sparsity.
     */
    static constexpr bool is_sparse = true;
};

}

#endif
