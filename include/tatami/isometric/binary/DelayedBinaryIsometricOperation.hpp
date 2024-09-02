#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OPERATION_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OPERATION_H

#include "../../base/Matrix.hpp"
#include "../../utils/new_extractor.hpp"
#include "../../utils/copy.hpp"
#include "../../dense/SparsifiedWrapper.hpp"
#include "../depends_utils.hpp"

#include <memory>
#include <vector>
#include <type_traits>

/**
 * @file DelayedBinaryIsometricOperation.hpp
 *
 * @brief Delayed binary isometric operations.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedBinaryIsometricOperation_internal {

/********************
 *** Dense simple ***
 ********************/

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseSimpleFull : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseSimpleFull(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        my_left_ext = new_extractor<false, oracle_>(left, my_row, oracle, opt);
        my_right_ext = new_extractor<false, oracle_>(right, my_row, std::move(oracle), opt);
        my_extent = my_row ? left->ncol() : left->nrow();

        my_right_holding_buffer.resize(my_extent);
        if constexpr(!same_value) {
            my_left_holding_buffer.resize(my_extent);
        }
    }

    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto rptr = my_right_ext->fetch(i, my_right_holding_buffer.data());

        if constexpr(same_value) {
            auto lptr = my_left_ext->fetch(i, buffer);
            copy_n(lptr, my_extent, buffer);
            my_operation.dense(my_row, my_oracle.get(i), static_cast<Index_>(0), my_extent, buffer, rptr, buffer);
        } else {
            auto lptr = my_left_ext->fetch(i, my_left_holding_buffer.data());
            my_operation.dense(my_row, my_oracle.get(i), static_cast<Index_>(0), my_extent, lptr, rptr, buffer);
        }

        return buffer;
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    std::unique_ptr<DenseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;
    Index_ my_extent;
    std::vector<InputValue_> my_right_holding_buffer;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_left_holding_buffer;
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseSimpleBlock : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseSimpleBlock(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_block_start(block_start),
        my_block_length(block_length)
    {
        my_left_ext = new_extractor<false, oracle_>(left, my_row, oracle, my_block_start, my_block_length, opt);
        my_right_ext = new_extractor<false, oracle_>(right, my_row, std::move(oracle), my_block_start, my_block_length, opt);

        my_right_holding_buffer.resize(my_block_length);
        if constexpr(!same_value) {
            my_left_holding_buffer.resize(my_block_length);
        }
    }

    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto rptr = my_right_ext->fetch(i, my_right_holding_buffer.data());

        if constexpr(same_value) {
            auto lptr = my_left_ext->fetch(i, buffer);
            copy_n(lptr, my_block_length, buffer);
            my_operation.dense(my_row, my_oracle.get(i), my_block_start, my_block_length, buffer, rptr, buffer);
        } else {
            auto lptr = my_left_ext->fetch(i, my_left_holding_buffer.data());
            my_operation.dense(my_row, my_oracle.get(i), my_block_start, my_block_length, lptr, rptr, buffer);
        }

        return buffer;
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    Index_ my_block_start, my_block_length;
    std::unique_ptr<DenseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;
    std::vector<InputValue_> my_right_holding_buffer;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_left_holding_buffer;
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseSimpleIndex : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseSimpleIndex(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_indices_ptr(std::move(indices_ptr))
    {
        my_left_ext = new_extractor<false, oracle_>(left, my_row, oracle, my_indices_ptr, opt);
        my_right_ext = new_extractor<false, oracle_>(right, my_row, std::move(oracle), my_indices_ptr, opt);

        my_right_holding_buffer.resize(my_indices_ptr->size());
        if constexpr(!same_value) {
            my_left_holding_buffer.resize(my_indices_ptr->size());
        }
    }

    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto rptr = my_right_ext->fetch(i, my_right_holding_buffer.data());
        const auto& indices = *my_indices_ptr;

        if constexpr(same_value) {
            auto lptr = my_left_ext->fetch(i, buffer);
            copy_n(lptr, indices.size(), buffer);
            my_operation.dense(my_row, my_oracle.get(i), indices, buffer, rptr, buffer);
        } else {
            auto lptr = my_left_ext->fetch(i, my_left_holding_buffer.data());
            my_operation.dense(my_row, my_oracle.get(i), indices, lptr, rptr, buffer);
        }

        return buffer;
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    VectorPtr<Index_> my_indices_ptr;
    std::unique_ptr<DenseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;
    std::vector<InputValue_> my_right_holding_buffer;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_left_holding_buffer;
};

/**********************
 *** Dense expanded ***
 **********************/

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseExpandedFull : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseExpandedFull(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Options opt) :
        my_operation(op),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        opt.sparse_ordered_index = true;
        my_left_ext = new_extractor<true, oracle_>(left, my_row, oracle, opt);
        my_right_ext = new_extractor<true, oracle_>(right, my_row, std::move(oracle), opt);

        my_extent = my_row ? left->ncol() : left->nrow();
        my_left_vbuffer.resize(my_extent);
        my_right_vbuffer.resize(my_extent);
        my_output_vbuffer.resize(my_extent);
        my_left_ibuffer.resize(my_extent);
        my_right_ibuffer.resize(my_extent);
        my_output_ibuffer.resize(my_extent);
    }

    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto lres = my_left_ext->fetch(i, my_left_vbuffer.data(), my_left_ibuffer.data());
        auto rres = my_right_ext->fetch(i, my_right_vbuffer.data(), my_right_ibuffer.data());

        i = my_oracle.get(i);
        auto num = my_operation.sparse(my_row, i, lres, rres, my_output_vbuffer.data(), my_output_ibuffer.data(), true, true);

        // Avoid calling my_operation.fill() if possible, as this might throw
        // zero-related errors with non-IEEE-float types.
        if (num < my_extent) { 
            std::fill_n(buffer, my_extent, my_operation.template fill<OutputValue_, InputValue_>(my_row, i));
        }

        for (Index_ j = 0; j < num; ++j) {
            buffer[my_output_ibuffer[j]] = my_output_vbuffer[j];
        }
        return buffer;
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;
    Index_ my_extent;

    std::vector<InputValue_> my_left_vbuffer, my_right_vbuffer;
    std::vector<OutputValue_> my_output_vbuffer;
    std::vector<Index_> my_left_ibuffer, my_right_ibuffer, my_output_ibuffer;
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseExpandedBlock : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseExpandedBlock(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_block_start(block_start),
        my_block_length(block_length)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        opt.sparse_ordered_index = true;
        my_left_ext = new_extractor<true, oracle_>(left, my_row, oracle, my_block_start, my_block_length, opt);
        my_right_ext = new_extractor<true, oracle_>(right, my_row, std::move(oracle), my_block_start, my_block_length, opt);

        my_left_vbuffer.resize(my_block_length);
        my_right_vbuffer.resize(my_block_length);
        my_output_vbuffer.resize(my_block_length);
        my_left_ibuffer.resize(my_block_length);
        my_right_ibuffer.resize(my_block_length);
        my_output_ibuffer.resize(my_block_length);
    }

    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto lres = my_left_ext->fetch(i, my_left_vbuffer.data(), my_left_ibuffer.data());
        auto rres = my_right_ext->fetch(i, my_right_vbuffer.data(), my_right_ibuffer.data());

        i = my_oracle.get(i);
        auto num = my_operation.sparse(my_row, i, lres, rres, my_output_vbuffer.data(), my_output_ibuffer.data(), true, true);

        // Avoid calling my_operation.fill() if possible, as this might throw
        // zero-related errors with non-IEEE-float types.
        if (num < my_block_length) { 
            std::fill_n(buffer, my_block_length, my_operation.template fill<OutputValue_, InputValue_>(my_row, i));
        }

        for (Index_ j = 0; j < num; ++j) {
            buffer[my_output_ibuffer[j] - my_block_start] = my_output_vbuffer[j];
        }
        return buffer;
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_block_start, my_block_length;

    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;

    std::vector<InputValue_> my_left_vbuffer, my_right_vbuffer;
    std::vector<OutputValue_> my_output_vbuffer;
    std::vector<Index_> my_left_ibuffer, my_right_ibuffer, my_output_ibuffer;
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseExpandedIndex : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseExpandedIndex(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_extent(indices_ptr->size())
    {
        // Create a remapping vector to map the extracted indices back to the
        // dense buffer. We use the 'remapping_offset' to avoid allocating the
        // full extent of the dimension.
        const auto& indices = *indices_ptr;
        if (my_extent) {
            my_remapping_offset = indices.front();
            my_remapping.resize(indices.back() - my_remapping_offset + 1);

            for (Index_ i = 0; i < my_extent; ++i) {
                my_remapping[indices[i] - my_remapping_offset] = i;
            }
        }

        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        opt.sparse_ordered_index = true;
        my_left_ext = new_extractor<true, oracle_>(left, my_row, oracle, indices_ptr, opt);
        my_right_ext = new_extractor<true, oracle_>(right, my_row, std::move(oracle), std::move(indices_ptr), opt);

        my_left_vbuffer.resize(my_extent);
        my_right_vbuffer.resize(my_extent);
        my_output_vbuffer.resize(my_extent);
        my_left_ibuffer.resize(my_extent);
        my_right_ibuffer.resize(my_extent);
        my_output_ibuffer.resize(my_extent);
    }

    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto lres = my_left_ext->fetch(i, my_left_vbuffer.data(), my_left_ibuffer.data());
        auto rres = my_right_ext->fetch(i, my_right_vbuffer.data(), my_right_ibuffer.data());

        i = my_oracle.get(i);
        auto num = my_operation.sparse(my_row, i, lres, rres, my_output_vbuffer.data(), my_output_ibuffer.data(), true, true);

        // Avoid calling my_operation.fill() if possible, as this might throw
        // zero-related errors with non-IEEE-float types.
        if (num < my_extent) { 
            std::fill_n(buffer, my_extent, my_operation.template fill<OutputValue_, InputValue_>(my_row, i));
        }

        for (Index_ j = 0; j < num; ++j) {
            buffer[my_remapping[my_output_ibuffer[j] - my_remapping_offset]] = my_output_vbuffer[j];
        }
        return buffer;
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_extent;

    std::vector<Index_> my_remapping;
    Index_ my_remapping_offset = 0;

    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;

    std::vector<InputValue_> my_left_vbuffer, my_right_vbuffer;
    std::vector<OutputValue_> my_output_vbuffer;
    std::vector<Index_> my_left_ibuffer, my_right_ibuffer, my_output_ibuffer;
};

/**************
 *** Sparse ***
 **************/

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class Sparse : public SparseExtractor<oracle_, OutputValue_, Index_> {
public:
    Sparse(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        initialize(my_row ? left->ncol() : left->nrow(), opt);
        my_left_ext = new_extractor<true, oracle_>(left, my_row, oracle, opt);
        my_right_ext = new_extractor<true, oracle_>(right, my_row, std::move(oracle), opt);
    }

    Sparse(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        initialize(block_length, opt);
        my_left_ext = new_extractor<true, oracle_>(left, my_row, oracle, block_start, block_length, opt);
        my_right_ext = new_extractor<true, oracle_>(right, my_row, std::move(oracle), block_start, block_length, opt);
    }

    Sparse(
        const Matrix<InputValue_, Index_>* left,
        const Matrix<InputValue_, Index_>* right,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        initialize(indices_ptr->size(), opt); // do this before the move.
        my_left_ext = new_extractor<true, oracle_>(left, my_row, oracle, indices_ptr, opt);
        my_right_ext = new_extractor<true, oracle_>(right, my_row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    void initialize(size_t extent, Options& opt) {
        my_report_value = opt.sparse_extract_value;
        my_report_index = opt.sparse_extract_index;

        my_left_ibuffer.resize(extent);
        my_right_ibuffer.resize(extent);
        if (my_report_value) {
            my_left_vbuffer.resize(extent);
            my_right_vbuffer.resize(extent);
        }

        opt.sparse_ordered_index = true;
        opt.sparse_extract_index = true;
    }

public:
    SparseRange<OutputValue_, Index_> fetch(Index_ i, OutputValue_* value_buffer, Index_* index_buffer) {
        auto left_ranges = my_left_ext->fetch(i, my_left_vbuffer.data(), my_left_ibuffer.data());
        auto right_ranges = my_right_ext->fetch(i, my_right_vbuffer.data(), my_right_ibuffer.data());
        auto num = my_operation.sparse(
            my_row, 
            my_oracle.get(i), 
            left_ranges, 
            right_ranges, 
            value_buffer,
            index_buffer,
            my_report_value,
            my_report_index
        );

        return SparseRange(
            num, 
            (my_report_value ? value_buffer: NULL),
            (my_report_index ? index_buffer: NULL)
        );
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_left_ext, my_right_ext;
    std::vector<InputValue_> my_left_vbuffer, my_right_vbuffer;
    std::vector<Index_> my_left_ibuffer, my_right_ibuffer;

    bool my_report_value = false;
    bool my_report_index = false;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed isometric operations on two matrices
 *
 * Implements any operation that takes two matrices of the same shape and returns another matrix of that shape.
 * Each entry of the output matrix is a function of the corresponding values in the two input matrices.
 * This operation is "delayed" in that it is only evaluated during data extraction, e.g., with `MyopicDenseExtractor::fetch()` or friends.
 *
 * This class is inspired by the `DelayedNaryIsoOp` class in the **DelayedArray** Bioconductor package.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * This is the type of the value of the output matrix.
 * @tparam InputValue_ Type of the value of the input matrices, to use in the operation.
 * This may or may not be the same as `OutputValue_`, depending on the methods available in `Operation_`.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class implementing the operation.
 * This should implement the same methods as `DelayedBinaryIsometricMockBasic` or `DelayedBinaryIsometricMockAdvanced`,
 * depending on whether it can take advantage of matrix sparsity.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DelayedBinaryIsometricOperation : public Matrix<OutputValue_, Index_> {
public:
    /**
     * @param left Pointer to the left matrix.
     * @param right Pointer to the right matrix.
     * @param operation Instance of the functor class.
     */
    DelayedBinaryIsometricOperation(std::shared_ptr<const Matrix<InputValue_, Index_> > left, std::shared_ptr<const Matrix<InputValue_, Index_> > right, Operation_ operation) : 
        my_left(std::move(left)), my_right(std::move(right)), my_operation(std::move(operation)) 
    {
        if (my_left->nrow() != my_right->nrow() || my_left->ncol() != my_right->ncol()) {
            throw std::runtime_error("shape of the left and right matrices should be the same");
        }

        my_prefer_rows_proportion = (my_left->prefer_rows_proportion() + my_right->prefer_rows_proportion()) / 2;

        if constexpr(!Operation_::is_basic) {
            if (my_operation.is_sparse()) {
                my_is_sparse = my_left->is_sparse() && my_right->is_sparse();

                // Well, better than nothing, I guess.
                my_is_sparse_proportion = (my_left->is_sparse_proportion() + my_right->is_sparse_proportion())/2;
            }
        }
    }

private:
    std::shared_ptr<const Matrix<InputValue_, Index_> > my_left, my_right;
    Operation_ my_operation;

    double my_prefer_rows_proportion;
    double my_is_sparse_proportion = 0;
    bool my_is_sparse = false;

public:
    Index_ nrow() const {
        return my_left->nrow();
    }

    Index_ ncol() const {
        return my_left->ncol();
    }

    bool is_sparse() const {
        return my_is_sparse;
    }

    double is_sparse_proportion() const {
        return my_is_sparse_proportion;
    }

    bool prefer_rows() const { 
        return my_prefer_rows_proportion > 0.5;
    }

    double prefer_rows_proportion() const { 
        return my_prefer_rows_proportion;
    }

    bool uses_oracle(bool row) const {
        return my_left->uses_oracle(row) || my_right->uses_oracle(row);
    }

    using Matrix<OutputValue_, Index_>::dense;

    using Matrix<OutputValue_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_simple_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOperation_internal::DenseSimpleFull<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_left.get(),
            my_right.get(),
            my_operation,
            row, 
            std::move(oracle),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_simple_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOperation_internal::DenseSimpleBlock<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_left.get(),
            my_right.get(),
            my_operation,
            row, 
            std::move(oracle),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_simple_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOperation_internal::DenseSimpleIndex<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_left.get(),
            my_right.get(),
            my_operation,
            row, 
            std::move(oracle),
            std::move(indices_ptr),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOperation_internal::DenseExpandedFull<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_left.get(),
            my_right.get(),
            my_operation,
            row, 
            std::move(oracle),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOperation_internal::DenseExpandedBlock<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_left.get(),
            my_right.get(),
            my_operation,
            row, 
            std::move(oracle),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOperation_internal::DenseExpandedIndex<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_left.get(),
            my_right.get(),
            my_operation,
            row, 
            std::move(oracle),
            std::move(indices_ptr),
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_internal(bool row, Args_&& ... args) const {
        if constexpr(!Operation_::is_basic) {
            if (my_left->is_sparse() && my_right->is_sparse()) {
                if (DelayedIsometricOperation_internal::can_dense_expand(my_operation, row)) {
                    return dense_expanded_internal<oracle_>(row, std::forward<Args_>(args)...);
                }
            }
        } 

        return dense_simple_internal<oracle_>(row, std::forward<Args_>(args)...);
    }

public:
    std::unique_ptr<MyopicDenseExtractor<OutputValue_, Index_> > dense(bool row, const Options& opt) const {
        return dense_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<OutputValue_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<OutputValue_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        if constexpr(!Operation_::is_basic) {
            if (my_is_sparse) {
                return std::make_unique<DelayedBinaryIsometricOperation_internal::Sparse<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
                    my_left.get(),
                    my_right.get(),
                    my_operation,
                    row, 
                    std::move(oracle),
                    opt
                );
            }
        } 

        return std::make_unique<FullSparsifiedWrapper<oracle_, OutputValue_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), opt),
            row ? my_left->ncol() : my_left->nrow(),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        if constexpr(!Operation_::is_basic) {
            if (my_is_sparse) {
                return std::make_unique<DelayedBinaryIsometricOperation_internal::Sparse<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
                    my_left.get(),
                    my_right.get(),
                    my_operation,
                    row, 
                    std::move(oracle),
                    block_start,
                    block_length,
                    opt
                );
            }
        }

        return std::make_unique<BlockSparsifiedWrapper<oracle_, OutputValue_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), block_start, block_length, opt),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if constexpr(!Operation_::is_basic) {
            if (my_is_sparse) {
                return std::make_unique<DelayedBinaryIsometricOperation_internal::Sparse<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
                    my_left.get(),
                    my_right.get(),
                    my_operation,
                    row, 
                    std::move(oracle),
                    std::move(indices_ptr),
                    opt
                );
            }
        }

        return std::make_unique<IndexSparsifiedWrapper<oracle_, OutputValue_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), indices_ptr, opt),
            indices_ptr,
            opt
        );
    }

public:
    std::unique_ptr<MyopicSparseExtractor<OutputValue_, Index_> > sparse(bool row, const Options& opt) const {
        return sparse_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<OutputValue_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<OutputValue_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<OracularDenseExtractor<OutputValue_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<OutputValue_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<OutputValue_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
public:
    std::unique_ptr<OracularSparseExtractor<OutputValue_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<OutputValue_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<OutputValue_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam OutputValue_ Type of the result of the operation.
 * @tparam InputValue_ Type of the value of the input matrices.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class defining the operation.
 *
 * @param left Pointer to a (possibly `const`) `Matrix`.
 * @param right Pointer to a (possibly `const`) `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedBinaryIsometricOperation` clas.
 */
template<typename OutputValue_ = double, typename InputValue_, typename Index_, class Operation_>
std::shared_ptr<Matrix<OutputValue_, Index_> > make_DelayedBinaryIsometricOperation(std::shared_ptr<const Matrix<InputValue_, Index_> > left, std::shared_ptr<const Matrix<InputValue_, Index_> > right, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Operation2_;
    return std::shared_ptr<Matrix<OutputValue_, Index_> >(new DelayedBinaryIsometricOperation<OutputValue_, InputValue_, Index_, Operation2_>(std::move(left), std::move(right), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename OutputValue_ = double, typename InputValue_, typename Index_, class Operation_>
std::shared_ptr<Matrix<OutputValue_, Index_> > make_DelayedBinaryIsometricOperation(std::shared_ptr<Matrix<InputValue_, Index_> > left, std::shared_ptr<Matrix<InputValue_, Index_> > right, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Operation2_;
    return std::shared_ptr<Matrix<OutputValue_, Index_> >(new DelayedBinaryIsometricOperation<OutputValue_, InputValue_, Index_, Operation2_>(std::move(left), std::move(right), std::move(op)));
}
/**
 * @endcond
 */

}

#endif
