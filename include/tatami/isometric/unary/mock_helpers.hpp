#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_HELPER_INTERFACE_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_HELPER_INTERFACE_H

#include <vector>
#include "../../base/SparseRange.hpp"

namespace tatami {

/**
 * @brief Mock helper for dense operations in `DelayedUnaryIsometricOp`.
 *
 * This defines the expectations for operations that discard sparsity in a variable manner
 * (i.e., zeros are transformed into non-zeros of different value depending on their position in the `Matrix`).
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedUnaryMockVariableDenseHelper {
    /**
     * This method should apply the operation to values in `buffer`,
     * representing a contiguous block of values from a row/column.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `buffer` contains the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param start Start of the contiguous block of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in,out] buffer Contents of the row/column extracted from the matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
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
        Value_* buffer)
    const {
        // Just filling it with something as a mock.
        std::fill_n(buffer, length, 0);
    }

    /**
     * This method should apply the operation to values in `buffer`,
     * representing an indexed subset of values from a row/column.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `buffer` contains the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param indices Sorted and unique indices of columns (if `row = true`) or rows (otherwise) extracted from `i`.
     * @param length Length of the contiguous block.
     * @param[in,out] buffer Contents of the row/column extracted from the matrix.
     * This has `length` addressable elements, and the result of the operation should be stored here.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    void dense(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        const std::vector<Index_>& indices, 
        Value_* buffer)
    const {
        std::fill_n(buffer, indices.size(), 0);
    }
 
    /**
     * This operation does not preserve sparsity.
     */
    static constexpr bool is_sparse = false;

    /**
     * Conversion of zeros to non-zero values is dependent on the row of origin.
     * This should be true, otherwise a `DelayeUnaryMockConstantDenseHelper` is expected.
     */
    static constexpr bool zero_depends_on_row = true;

    /**
     * Conversion of zeros to non-zero values is dependent on the column of origin.
     * This should be true, otherwise a `DelayeUnaryMockConstantDenseHelper` is expected.
     */
    static constexpr bool zero_depends_on_column = true;

    /**
     * @cond
     */
    virtual ~DelayedUnaryMockVariableDenseHelper() = default;
    /**
     * @endcond
     */
};

/**
 * @brief Mock helper for constant dense operations in `DelayedUnaryIsometricOp`.
 *
 * This defines the expectations for operations that discard sparsity in a variable manner
 * (i.e., zeros are transformed into non-zeros of different value depending on their position in the `Matrix`).
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedUnaryMockConstantDenseHelper : public DelayedUnaryMockVariableDenseHelper {
    /**
     * This method applies the operation to a sparse range representing the contents of a row/column from the underyling matrix.
     * Specifically, the operation only needs to be applied to the structural non-zeros, as the zeros are populated by `fill()`.
     * The results of the operation should then be stored in the `output_*` buffers.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `buffer` contains the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param num Number of non-zero elements for row/column `i`.
     * @param[in,out] value Pointer to an array of values of the non-zero elements.
     * This is guaranteed to have `num` addressable elements.
     * @param[in] index Pointer to an array of column/row indices of the non-zero elements.
     * This is guaranteed to have `num` addressable elements.
     * Note that indices are not guaranteed to be sorted.
     *
     * This method is expected to iterate over `value` and modify it in place,
     * i.e., replace each value with the result of the operation on that value.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    void sparse(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        Index_ num,
        Value_* value,
        [[maybe_unused]] const Index_* index)
    const {
        std::fill(value, value + num, 0);
    }

    /**
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param i The index of the row containing the zero, if `zero_depends_on_row = true`;
     * the index of the column containing the zero, if `zero_depends_on_column = true`;
     * or ignored, if both `zero_depends_on_row` and `zero_depends_on_column` are both false.
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
     * Conversion of zeros to non-zero values is not dependent on the row of origin.
     * This may also be `true` provided that `zero_depends_on_column = false`.
     */
    static constexpr bool zero_depends_on_row = false;

    /**
     * Conversion of zeros to non-zero values is not dependent on the column of origin.
     * This may also be `true` provided that `zero_depends_on_row = false`.
     */
    static constexpr bool zero_depends_on_column = false;
};

/**
 * @brief Interface for sparse operations in `DelayedUnaryIsometricOp`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedUnaryMockSparseHelper : public DelayedUnaryMockVariableDenseHelper {
    /**
     * This method applies the operation to a sparse range representing the contents of a row/column from the underlying matrix.
     * It should only be called for sparsity-preserving operations, enabling optimizations by skipping structural zeros.
     * The results of the operation should then be stored in the `output_*` buffers.
     *
     * @tparam Value_ Type of matrix value.
     * @tparam Index_ Type of index value.
     *
     * @param row Whether `buffer` contains the row contents.
     * @param i Index of the extracted row (if `row = true`) or column (otherwise).
     * @param[in,out] value Pointer to an array of values of the non-zero elements.
     * This is guaranteed to have `num` addressable elements.
     * @param[in] index Pointer to an array of column/row indices of the non-zero elements.
     * Alternatively this may be NULL.
     *
     * This method is expected to iterate over `value` and modify it in place,
     * i.e., replace each value with the result of the operation on that value.
     *
     * If `depends_on_row && !row` or `depends_on_column && row`, `index` is guaranteed to be non-NULL.
     * Otherwise, it may be NULL and should be ignored.
     * If supplied, indices are not guaranteed to be sorted.
     *
     * Note that this method does not necessarily need to have the same template arguments.
     * It will be called without any explicit template arguments so anything can be used as long as type deduction works.
     */
    template<typename Value_, typename Index_>
    void sparse(
        [[maybe_unused]] bool row, 
        [[maybe_unused]] Index_ i, 
        Index_ num,
        Value_* value,
        const Index_* index)
    const {
        std::fill(value, value + num, 0);
    }

    /**
     * This operation preserves sparsity.
     */
    static constexpr bool is_sparse = true;

    /**
     * Whether the operation requires the identity of the row of origin.
     */
    static constexpr bool depends_on_row = false;

    /**
     * Whether the operation requires the identity of the column of origin.
     */
    static constexpr bool depends_on_column = false;
};

}

#endif
