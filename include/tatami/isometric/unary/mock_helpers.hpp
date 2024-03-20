#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_HELPER_INTERFACE_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_HELPER_INTERFACE_H

#include <vector>
#include "../../base/SparseRange.hpp"

/** 
 * @file mock_helpers.hpp
 * @brief Expectations for `tatami::DelayedUnaryIsometricOp` helpers.
 */

namespace tatami {

/**
 * @brief Basic mock operation for a `DelayedUnaryIsometricOp`. 
 *
 * This defines the basic expectations for an operation to use in `DelayedUnaryIsometricOp`.
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedUnaryBasicMockHelper {
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
    virtual ~DelayedUnaryBasicMockHelper() = default;
    /**
     * @endcond
     */
};

/**
 * @brief Advanced mock operation for `DelayedUnaryIsometricOp`.
 *
 * This class defines the advanced expectations for an operation in `DelayedUnaryIsometricOp`,
 * which improves efficiency by taking advantage of any sparsity in the underlying matrix.
 * Either the operation itself preserves sparsity, or any loss of sparsity is predictable,
 * i.e., zeros are transformed into a constant non-zero value that does not depend on its position in the `Matrix`.
 *
 * Actual operations aren't expected to inherit from this class;
 * this is only provided for documentation purposes.
 * Operations only need to implement methods with the same signatures for compile-time polymorphism.
 */
struct DelayedUnaryAdvancedMockHelper : public DelayedUnaryBasicMockHelper {
    /**
     * This method applies the operation to a sparse range representing the contents of a row/column from the underyling matrix.
     * Specifically, the operation only needs to be applied to the structural non-zeros;
     * structural zeros are either ignored for sparsity-preserving operations,
     * or the result of the operation on zeros will be populated by `fill()`.
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
     * Alternatively NULL.
     *
     * This method is expected to iterate over `value` and modify it in place,
     * i.e., replace each value with the result of the operation on that value.
     *
     * If `non_zero_depends_on_row && !row` or `non_zero_depends_on_column && row`, `index` is guaranteed to be non-NULL.
     * Otherwise, it may be NULL and should be ignored.
     * Even if non-NULL, indices are not guaranteed to be sorted.
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
     * Conversion of zeros to non-zero values is not dependent on the row of origin.
     * This may also be `true` provided that `zero_depends_on_column = false`.
     */
    static constexpr bool zero_depends_on_row = false;

    /**
     * Conversion of zeros to non-zero values is not dependent on the column of origin.
     * This may also be `true` provided that `zero_depends_on_row = false`.
     */
    static constexpr bool zero_depends_on_column = false;

    /**
     * Whether the operation requires the identity of the row of origin.
     * This only affects the presence of `index` in `sparse()`.
     * May be true or false.
     */
    static constexpr bool non_zero_depends_on_row = false;

    /**
     * Whether the operation requires the identity of the column of origin.
     * This only affects the presence of `index` in `sparse()`.
     * May be true or false.
     */
    static constexpr bool non_zero_depends_on_column = false;

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
