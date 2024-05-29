#ifndef TATAMI_ISOMETRIC_BINARY_MOCK_HELPERS_H
#define TATAMI_ISOMETRIC_BINARY_MOCK_HELPERS_H

#include <vector>
#include "../../base/SparseRange.hpp"

/** 
 * @file mock_helpers.hpp
 * @brief Expectations for `tatami::DelayedBinaryIsometricOperation` helpers.
 */

namespace tatami {

/**
 * @brief Basic mock operation for `DelayedBinaryIsometricOperation`.
 *
 * This class defines the basic expectations for an operation in `DelayedBinaryIsometricOperation`.
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
class DelayedBinaryIsometricMockBasic {
public:
    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
     * These buffers represent the same element of the target dimension from the left and right matrices, respectively, in dense form.
     * Each buffer holds values from a contiguous block of the non-target dimension.
     *
     * Implementations of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to be used in the operation.
     * @tparam OutputValue_ Type of the result of the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `left_buffer` and `right_buffer` hold the contents of the `i`-th row from both matrices;
     * otherwise, they hold the contents of the `i`-th column.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * Unlike `DelayedBinaryIsometricMockAdvanced::dense()`, this is always guaranteed to be the actual index and not a placeholder.
     * @param start Start of the contiguous block of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in] left_buffer Pointer to an array containing the row/column extracted from the left matrix.
     * This has `length` addressable elements. 
     * @param[in] right_buffer Pointer to an array containing the row/column extracted from the right matrix.
     * This has `length` addressable elements.
     * @param[out] output_buffer Pointer to an array in which to store the result of the operation.
     * This has `length` addressable elements.
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `left_buffer`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] Index_ start, 
        Index_ length, 
        [[maybe_unused]] const InputValue_* left_buffer, 
        [[maybe_unused]] const InputValue_* right_buffer,
        OutputValue_* output_buffer)
    const {
        // Just filling it with something as a mock.
        std::fill_n(output_buffer, length, 0);
    }

    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
     * These buffers represent the same element of the target dimension from the left and right matrices, respectively, in dense form.
     * Each buffer holds values from an indexed subset of the non-target dimension.
     * 
     * Implementations of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to be used in the operation.
     * @tparam OutputValue_ Type of the result of the operation.
     *
     * @param row Whether the rows are the target dimension.
     * If true, `left_buffer` and `right_buffer` hold the contents of the `i`-th row from both matrices;
     * otherwise, they hold the contents of the `i`-th column.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * Unlike `DelayedBinaryIsometricMockAdvanced::dense()`, this is always guaranteed to be the actual index and not a placeholder.
     * @param indices Sorted and unique indices of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param[in] left_buffer Pointer to an array containing the row/column extracted from the left matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     * @param[in] right_buffer Pointer to an array containing the row/column extracted from the right matrix.
     * This has `length` addressable elements.
     * @param[out] output_buffer Pointer to an array in which to store the result of the operation.
     * This has `length` addressable elements.
     * If `InputValue_ == OutputValue_`, this is guaranteed to be the same as `left_buffer`.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        const std::vector<Index_>& indices, 
        [[maybe_unused]] const InputValue_* left_buffer, 
        [[maybe_unused]] const InputValue_* right_buffer, 
        OutputValue_* output_buffer) 
    const {
        std::fill_n(output_buffer, indices.size(), 0);
    }

    /**
     * Whether this is a basic operation.
     * This should be true, otherwise an advanced operation is expected (see `DelayedBinaryIsometricMockAdvanced`).
     */
    static constexpr bool is_basic = true;
};

/**
 * @brief Advanced mock operation for `DelayedBinaryIsometricOperation`.
 *
 * This class defines the advanced expectations for an operation in `DelayedBinaryIsometricOperation`,
 * which improves efficiency by taking advantage of any sparsity in the underlying matrices.
 * Either the operation itself preserves sparsity, or any loss of sparsity is predictable,
 * i.e., zeros are transformed into a constant non-zero value that does not depend on its position in the `Matrix`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
class DelayedBinaryIsometricMockAdvanced {
public:
    /**
     * This method will be called with explicit `OutputValue_` and `InputValue_` template parameters.
     * Implementations of this method should either ensure that `Index_` is deducible or use a fixed integer type in the method signature.
     *
     * @tparam OutputValue_ Type of the result of the operation.
     * @tparam InputValue_ Type of the matrix value used in the operation.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `i` refers to the row or column index.
     * @param i The index of the row (if `row = true`) or column (otherwise) containing the zeros.
     * This argument should be ignored if the operation does not depend on the row/column (i.e., when all of `zero_depends_on_row()` and friends return false),
     * in which case an arbitrary placeholder may be supplied. 
     *
     * @return The result of `OP(lz, rz)` where `OP` is the operation,
     * `lz` is a structural zero from the `i`-th row/column of the left matrix,
     * and `rz` is a structural zero from the `i`-th row/column of the left matrix,
     */
    template<typename OutputValue_, typename InputValue_, typename Index_>
    OutputValue_ fill([[maybe_unused]] bool row, [[maybe_unused]] Index_ i) const { 
        return 0;
    }

    /**
     * Whether this is a basic operation.
     * This should be false, otherwise a basic operation interface is expected (see `DelayedBinaryIsometricMockBasic`).
     */
    static constexpr bool is_basic = false;

    /**
     * @return Whether applying the operation to a pair of structural zeros (one from each matrix)
     * yields a value that depends on the identity of the row containing those zeros.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedBinaryIsometricOperation` will automatically recognize such operations as being row-independent.
     *
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool zero_depends_on_row() const {
        return false;
    }

    /**
     * @return Whether applying the operation to a pair of structural zeros (one from each matrix)
     * yields a value that depends on the identity of the column containing those zeros.
     *
     * This method is only called when `is_sparse()` returns false.
     * It is not necessary to explicitly return `false` here for sparsity-preserving operations,
     * as `DelayedBinaryIsometricOperation` will automatically recognize such operations as being row-independent.
     *
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool zero_depends_on_column() const {
        return false;
    }

    /**
     * @return Whether the result of the operation depends on the identity of the row containing the operands,
     * where at least one of the operands is non-zero.
     * 
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool non_zero_depends_on_row() const {
        return false;
    }

    /**
     * @return Whether the result of the operation depends on the identity of the column containing the operands,
     * where at least one of the operands is non-zero.
     *
     * This method may be omitted from the class definition, in which case it is assumed to always return false. 
     */
    bool non_zero_depends_on_column() const {
        return false;
    }

    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
     * These buffers represent the same element of the target dimension from the left and right matrices, respectively, in dense form.
     * Each buffer holds values from a contiguous block of the non-target dimension.
     *
     * Implementations of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to be used in the operation.
     * @tparam OutputValue_ Type of the result of the operation.
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
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] Index_ start, 
        Index_ length, 
        [[maybe_unused]] const InputValue_* left_buffer, 
        [[maybe_unused]] const InputValue_* right_buffer,
        OutputValue_* output_buffer)
    const {
        // Just filling it with something as a mock.
        std::fill_n(output_buffer, length, 0);
    }

    /**
     * This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
     * These buffers represent the same element of the target dimension from the left and right matrices, respectively, in dense form.
     * Each buffer holds values from an indexed subset of the non-target dimension.
     *
     * Implementations of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to be used in the operation.
     * @tparam OutputValue_ Type of the result of the operation.
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
    template<typename Index_, typename InputValue_, typename OutputValue_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        const std::vector<Index_>& indices, 
        [[maybe_unused]] const InputValue_* left_buffer, 
        [[maybe_unused]] const InputValue_* right_buffer,
        OutputValue_* output_buffer) 
    const {
        std::fill_n(output_buffer, indices.size(), 0);
    }

    /**
     * This method applies the operation to the sparse ranges in `left` and `right`, 
     * containing values from the same element of the target dimension from the left and right matrices, respectively.
     * Specifically, the operation only needs to be applied to the structural non-zeros,
     * and results of the operation should be stored in the `output_*` buffers.
     * Structural zeros are either ignored for sparsity-preserving operations,
     * or the result of the operation on zeros will be populated by `fill()`.
     *
     * @tparam Index_ Type of index value.
     * @tparam InputValue_ Type of matrix value to be used in the operation.
     * @tparam OutputValue_ Type of the result of the operation.
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
     *
     * Note that implementations of this method do not necessarily need to have the same template arguments as shown here.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Index_, typename InputValue_, typename OutputValue_>
    Index_ sparse(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        [[maybe_unused]] const SparseRange<InputValue_, Index_>& left, 
        [[maybe_unused]] const SparseRange<InputValue_, Index_>& right, 
        [[maybe_unused]] OutputValue_* output_value, 
        [[maybe_unused]] Index_* output_index, 
        [[maybe_unused]] bool report_value, 
        [[maybe_unused]] bool report_index)
    const {
        return 0;
    }

    /** 
     * @return Whether this operation preserves sparsity.
     */
    bool is_sparse() const {
        return true;
    }
};

}

#endif
