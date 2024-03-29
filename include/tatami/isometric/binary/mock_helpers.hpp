#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OP_HELPER_INTERFACE_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OP_HELPER_INTERFACE_H

#include <vector>
#include "../../base/SparseRange.hpp"

/** 
 * @file mock_helpers.hpp
 * @brief Expectations for `tatami::DelayedBinaryIsometricOp` helpers.
 */

namespace tatami {

/**
 * @brief Basic mock operation for `DelayedBinaryIsometricOp`.
 *
 * This class defines the basic expectations for an operation in `DelayedBinaryIsometricOp`.
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryBasicMockHelper {
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
     * @param[in] right_buffer Contents of the row/column extracted from the right matrix.
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
     * @param[in,out] left_buffer Contents of the row/column extracted from the left matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     * @param[in] right_buffer Contents of the row/column extracted from the right matrix.
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
     * Conversion of zeros to non-zero values is dependent on rows.
     * This should be `true`, otherwise a `DelayedBinaryAdvancedMockHelper` is expected.
     */
    static constexpr bool zero_depends_on_row = true;

    /**
     * Conversion of zeros to non-zero values is dependent on columns.
     * This should be `true`, otherwise a `DelayedBinaryAdvancedMockHelper` is expected.
     */
    static constexpr bool zero_depends_on_column = true;

    /**
     * @cond
     */
    virtual ~DelayedBinaryBasicMockHelper() = default;
    /**
     * @endcond
     */
};

/**
 * @brief Advanced mock operation for `DelayedBinaryIsometricOp`.
 *
 * This class defines the advanced expectations for an operation in `DelayedBinaryIsometricOp`,
 * which improves efficiency by taking advantage of any sparsity in the underlying matrices.
 * Either the operation itself preserves sparsity, or any loss of sparsity is predictable,
 * i.e., zeros are transformed into a constant non-zero value that does not depend on its position in the `Matrix`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedBinaryAdvancedMockHelper : public DelayedBinaryBasicMockHelper {
    /**
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param i The index of the row containing the zero, if `zero_depends_on_row = true`;
     * the index of the column containing the zero, if `zero_depends_on_column = true`;
     * or ignored, if neither are true.
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
     * Conversion of zeros to non-zero values is not dependent on rows.
     * This may also be `true` provided that `zero_depends_on_column = false`.
     */
    static constexpr bool zero_depends_on_row = false;

    /**
     * Conversion of zeros to non-zero values is not dependent on columns.
     * This may also be `true` provided that `zero_depends_on_row = false`.
     */
    static constexpr bool zero_depends_on_column = false;

    /**
     * This method applies the operation to the sparse ranges in `left` and `right`, 
     * representing the contents of the same row/column from the left and right matrices respectively.
     * Specifically, the operation only needs to be applied to the structural non-zeros,
     * and results of the operation should be stored in the `output_*` buffers.
     * Structural zeros are either ignored for sparsity-preserving operations,
     * or the result of the operation on zeros will be populated by `fill()`.
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
        [[maybe_unused]] bool report_value, 
        [[maybe_unused]] bool report_index)
    const {
        return 0;
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
