#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OP_HELPER_INTERFACE_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OP_HELPER_INTERFACE_H

#include <vector>
#include "../../base/SparseRange.hpp"

namespace tatami {

/**
 * @brief Mock helper for dense operations in `DelayedBinaryIsometricOp`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryMockDenseHelper {
    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
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
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
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
     * This operation is not sparsity preserving. 
     */
    static constexpr bool is_sparse = false;

    /**
     * @cond
     */
    virtual ~DelayedBinaryMockDenseHelper() = default;
    /**
     * @endcond
     */
};

/**
 * @brief Interface for sparse operations in `DelayedBinaryIsometricOp`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryMockSparseHelper : public DelayedBinaryMockDenseHelper {
    /**
     * This method applies the operation to the sparse ranges in `left` and `right`, storing results in the `output_*` buffers.
     * It should only be called for sparsity-preserving operations, enabling optimizations by skipping structural zeros.
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
     * This operation is sparsity preserving. 
     */
    static constexpr bool is_sparse = true;
};

}

#endif
