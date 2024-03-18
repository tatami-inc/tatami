#ifndef TATAMI_DELAYED_BINARY_ISOMETRIC_OP_H
#define TATAMI_DELAYED_BINARY_ISOMETRIC_OP_H

#include "../../base/Matrix.hpp"
#include "../../utils/new_extractor.hpp"
#include "../../utils/copy.hpp"
#include "../../dense/SparsifiedWrapper.hpp"

#include <memory>
#include <vector>

/**
 * @file DelayedBinaryIsometricOp.hpp
 *
 * @brief Delayed binary isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedBinaryIsometricOp_internal {

/********************
 *** Dense simple ***
 ********************/

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseSimpleFull : public DenseExtractor<oracle_, Value_, Index_> {
    DenseSimpleFull(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt) :
        operation(op),
        row(row)
    {
        extent = row ? lmat->ncol() : lmat->nrow();
        left = new_extractor<false, oracle_>(lmat, row, oracle, opt);
        right = new_extractor<false, oracle_>(rmat, row, std::move(oracle), opt);
        holding_buffer.resize(extent);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto lptr = left->fetch(i, buffer);
        copy_n(lptr, extent, buffer);
        auto rptr = right->fetch(i, holding_buffer.data());
        operation.dense(row, i, static_cast<Index_>(0), extent, buffer, rptr);
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    Index_ extent;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> holding_buffer;
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseSimpleBlock : public DenseExtractor<oracle_, Value_, Index_> {
    DenseSimpleBlock(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        const Options& opt) :
        operation(op),
        row(row),
        block_start(block_start),
        block_length(block_length)
    {
        left = new_extractor<false, oracle_>(lmat, row, oracle, block_start, block_length, opt);
        right = new_extractor<false, oracle_>(rmat, row, std::move(oracle), block_start, block_length, opt);
        holding_buffer.resize(block_length);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto lptr = left->fetch(i, buffer);
        copy_n(lptr, block_length, buffer);
        auto rptr = right->fetch(i, holding_buffer.data());
        operation.dense(row, i, block_start, block_length, buffer, rptr);
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    Index_ block_start, block_length;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> holding_buffer;
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseSimpleIndex : public DenseExtractor<oracle_, Value_, Index_> {
    DenseSimpleIndex(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> idx_ptr,
        const Options& opt) :
        operation(op),
        row(row),
        indices_ptr(std::move(idx_ptr))
    {
        left = new_extractor<false, oracle_>(lmat, row, oracle, indices_ptr, opt);
        right = new_extractor<false, oracle_>(rmat, row, std::move(oracle), indices_ptr, opt);
        holding_buffer.resize(indices_ptr->size());
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto lptr = left->fetch(i, buffer);
        copy_n(lptr, holding_buffer.size(), buffer);
        auto rptr = right->fetch(i, holding_buffer.data());
        operation.dense(row, i, *indices_ptr, buffer, rptr);
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    VectorPtr<Index_> indices_ptr;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> holding_buffer;
};

/**********************
 *** Dense expanded ***
 **********************/

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseExpandedFull : public DenseExtractor<oracle_, Value_, Index_> {
    DenseExpandedFull(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        Operation_ op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Options opt) :
        operation(op),
        row(row)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        opt.sparse_ordered_index = true;
        left = new_extractor<true, oracle_>(lmat, row, oracle, opt);
        right = new_extractor<true, oracle_>(rmat, row, std::move(oracle), opt);

        extent = row ? lmat->ncol() : lmat->nrow();
        left_vbuffer.resize(extent);
        right_vbuffer.resize(extent);
        output_vbuffer.resize(extent);
        left_ibuffer.resize(extent);
        right_ibuffer.resize(extent);
        output_ibuffer.resize(extent);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto lres = left->fetch(i, left_vbuffer.data(), left_ibuffer.data());
        auto rres = right->fetch(i, right_vbuffer.data(), right_ibuffer.data());
        auto num = operation.sparse(row, i, lres, rres, output_vbuffer.data(), output_ibuffer.data(), true, true);

        std::fill(buffer, buffer + extent, operation.template fill<Value_>(i));
        for (Index_ j = 0; j < num; ++j) {
            buffer[output_ibuffer[j]] = output_vbuffer[j];
        }
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    Index_ extent;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> left_vbuffer, right_vbuffer, output_vbuffer;
    std::vector<Index_> left_ibuffer, right_ibuffer, output_ibuffer;
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseExpandedBlock : public DenseExtractor<oracle_, Value_, Index_> {
    DenseExpandedBlock(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        Options opt) :
        operation(op),
        row(row),
        block_start(block_start),
        block_length(block_length)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        opt.sparse_ordered_index = true;
        left = new_extractor<true, oracle_>(lmat, row, oracle, block_start, block_length, opt);
        right = new_extractor<true, oracle_>(rmat, row, std::move(oracle), block_start, block_length, opt);

        left_vbuffer.resize(block_length);
        right_vbuffer.resize(block_length);
        output_vbuffer.resize(block_length);
        left_ibuffer.resize(block_length);
        right_ibuffer.resize(block_length);
        output_ibuffer.resize(block_length);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto lres = left->fetch(i, left_vbuffer.data(), left_ibuffer.data());
        auto rres = right->fetch(i, right_vbuffer.data(), right_ibuffer.data());
        auto num = operation.sparse(row, i, lres, rres, output_vbuffer.data(), output_ibuffer.data(), true, true);

        std::fill(buffer, buffer + block_length, operation.template fill<Value_>(i));
        for (Index_ j = 0; j < num; ++j) {
            buffer[output_ibuffer[j] - block_start] = output_vbuffer[j];
        }
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    Index_ block_start, block_length;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> left_vbuffer, right_vbuffer, output_vbuffer;
    std::vector<Index_> left_ibuffer, right_ibuffer, output_ibuffer;
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseExpandedIndex : public DenseExtractor<oracle_, Value_, Index_> {
    DenseExpandedIndex(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        operation(op),
        row(row)
    {
        // Create a remapping vector to map the extracted indices back to the
        // dense buffer. We use the 'remapping_offset' to avoid allocating the
        // full extent of the dimension.
        const auto& indices = *indices_ptr;
        extent = indices.size();
        if (extent) {
            remapping_offset = indices.front();
            remapping.resize(indices.back() - remapping_offset + 1);
            for (Index_ i = 0; i < extent; ++i) {
                remapping[indices[i] - remapping_offset] = i;
            }
        }

        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        opt.sparse_ordered_index = true;
        left = new_extractor<true, oracle_>(lmat, row, oracle, indices_ptr, opt);
        right = new_extractor<true, oracle_>(rmat, row, std::move(oracle), std::move(indices_ptr), opt);

        left_vbuffer.resize(extent);
        right_vbuffer.resize(extent);
        output_vbuffer.resize(extent);
        left_ibuffer.resize(extent);
        right_ibuffer.resize(extent);
        output_ibuffer.resize(extent);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto lres = left->fetch(i, left_vbuffer.data(), left_ibuffer.data());
        auto rres = right->fetch(i, right_vbuffer.data(), right_ibuffer.data());
        auto num = operation.sparse(row, i, lres, rres, output_vbuffer.data(), output_ibuffer.data(), true, true);

        std::fill(buffer, buffer + extent, operation.template fill<Value_>(i));
        for (Index_ j = 0; j < num; ++j) {
            buffer[remapping[output_ibuffer[j] - remapping_offset]] = output_vbuffer[j];
        }
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    Index_ extent;
    std::vector<Index_> remapping;
    Index_ remapping_offset = 0;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> left_vbuffer, right_vbuffer, output_vbuffer;
    std::vector<Index_> left_ibuffer, right_ibuffer, output_ibuffer;
};

/**************
 *** Sparse ***
 **************/

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct Sparse : public SparseExtractor<oracle_, Value_, Index_> {
    Sparse(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Options opt) :
        operation(op),
        row(row)
    {
        initialize(row ? lmat->ncol() : lmat->nrow(), opt);
        left = new_extractor<true, oracle_>(lmat, row, oracle, opt);
        right = new_extractor<true, oracle_>(rmat, row, std::move(oracle), opt);
    }

    Sparse(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        Options opt) :
        operation(op),
        row(row)
    {
        initialize(block_length, opt);
        left = new_extractor<true, oracle_>(lmat, row, oracle, block_start, block_length, opt);
        right = new_extractor<true, oracle_>(rmat, row, std::move(oracle), block_start, block_length, opt);
    }

    Sparse(
        const Matrix<Value_, Index_>* lmat,
        const Matrix<Value_, Index_>* rmat,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        operation(op)
    {
        initialize(indices_ptr->size(), opt); // do this before the move.
        left = new_extractor<true, oracle_>(lmat, row, oracle, indices_ptr, opt);
        right = new_extractor<true, oracle_>(rmat, row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    void initialize(size_t extent, Options& opt) {
        report_value = opt.sparse_extract_value;
        report_index = opt.sparse_extract_index;

        left_internal_ibuffer.resize(extent);
        right_internal_ibuffer.resize(extent);
        if (report_value) {
            left_internal_vbuffer.resize(extent);
            right_internal_vbuffer.resize(extent);
        }

        opt.sparse_ordered_index = true;
        opt.sparse_extract_index = true;
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto left_ranges = left->fetch(i, left_internal_vbuffer.data(), left_internal_ibuffer.data());
        auto right_ranges = right->fetch(i, right_internal_vbuffer.data(), right_internal_ibuffer.data());
        auto num = operation.sparse(
            row, 
            i, 
            left_ranges, 
            right_ranges, 
            vbuffer,
            ibuffer,
            report_value,
            report_index
        );
        return SparseRange(num, (report_value ? vbuffer : NULL), (report_index ? ibuffer : NULL));
    }

private:
    const Operation_& operation;
    bool row;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> left_internal_vbuffer, right_internal_vbuffer;
    std::vector<Index_> left_internal_ibuffer, right_internal_ibuffer;
    bool report_value = false;
    bool report_index = false;
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
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `DenseExtractor::fetch()` or friends.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
 * This should implement the same methods as `DelayedBinaryBasicMockHelper` or `DelayedBinaryAdvancedMockHelper`,
 * depending on whether it can take advantage of matrix sparsity.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedBinaryIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param l Pointer to the left matrix.
     * @param r Pointer to the right matrix.
     * @param op Instance of the functor class.
     */
    DelayedBinaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > l, std::shared_ptr<const Matrix<Value_, Index_> > r, Operation_ op) : 
        left(std::move(l)), right(std::move(r)), operation(std::move(op)) 
    {
        if (left->nrow() != right->nrow() || left->ncol() != right->ncol()) {
            throw std::runtime_error("shape of the left and right matrices should be the same");
        }

        prefer_rows_proportion_internal = (left->prefer_rows_proportion() + right->prefer_rows_proportion()) / 2;

        if constexpr(is_advanced) {
            if (operation.is_sparse()) {
                is_sparse_internal = left->sparse() && right->sparse();

                // Well, better than nothing, I guess.
                sparse_proportion_internal = (left->sparse_proportion() + right->sparse_proportion())/2;
            }
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > left, right;
    Operation_ operation;
    double prefer_rows_proportion_internal;
    double sparse_proportion_internal = 0;
    bool is_sparse_internal = false;

    static constexpr bool is_advanced = (!Operation_::zero_depends_on_row || !Operation_::zero_depends_on_column);

public:
    Index_ nrow() const {
        return left->nrow();
    }

    Index_ ncol() const {
        return left->ncol();
    }

    bool sparse() const {
        return is_sparse_internal;
    }

    double sparse_proportion() const {
        return sparse_proportion_internal;
    }

    bool prefer_rows() const { 
        return prefer_rows_proportion_internal > 0.5;
    }

    double prefer_rows_proportion() const { 
        return prefer_rows_proportion_internal;
    }

    bool uses_oracle(bool row) const {
        return left->uses_oracle(row) || right->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_simple_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOp_internal::DenseSimpleFull<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_simple_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOp_internal::DenseSimpleBlock<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_simple_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOp_internal::DenseSimpleIndex<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            std::move(indices_ptr),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOp_internal::DenseExpandedFull<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOp_internal::DenseExpandedBlock<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedBinaryIsometricOp_internal::DenseExpandedIndex<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            std::move(indices_ptr),
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, Args_&& ... args) const {
        if constexpr(is_advanced) {
            if (left->sparse() && right->sparse()) {
                // If we don't depend on the rows, then we don't need row indices when 'row = false'.
                // Similarly, if we don't depend on columns, then we don't column row indices when 'row = true'.
                if ((!Operation_::zero_depends_on_row && !row) || (!Operation_::zero_depends_on_column && row)) {
                    return dense_expanded_internal<oracle_>(row, std::forward<Args_>(args)...);
                }
            }
        } 

        return dense_simple_internal<oracle_>(row, std::forward<Args_>(args)...);
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options& opt) const {
        return dense_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        if constexpr(is_advanced) {
            if (is_sparse_internal) {
                return std::make_unique<DelayedBinaryIsometricOp_internal::Sparse<oracle_, Value_, Index_, Operation_> >(
                    left.get(),
                    right.get(),
                    operation,
                    row, 
                    std::move(oracle),
                    opt
                );
            }
        } 

        return std::make_unique<FullSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), opt),
            row ? left->ncol() : left->nrow(),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        if constexpr(is_advanced) {
            if (is_sparse_internal) {
                return std::make_unique<DelayedBinaryIsometricOp_internal::Sparse<oracle_, Value_, Index_, Operation_> >(
                    left.get(),
                    right.get(),
                    operation,
                    row, 
                    std::move(oracle),
                    block_start,
                    block_length,
                    opt
                );
            }
        }

        return std::make_unique<BlockSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), block_start, block_length, opt),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if constexpr(is_advanced) {
            if (is_sparse_internal) {
                return std::make_unique<DelayedBinaryIsometricOp_internal::Sparse<oracle_, Value_, Index_, Operation_> >(
                    left.get(),
                    right.get(),
                    operation,
                    row, 
                    std::move(oracle),
                    std::move(indices_ptr),
                    opt
                );
            }
        }

        return std::make_unique<IndexSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), indices_ptr, opt),
            indices_ptr,
            opt
        );
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        return sparse_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return dense_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return sparse_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class defining the operation.
 *
 * @param left Pointer to a (possibly `const`) `Matrix`.
 * @param right Pointer to a (possibly `const`) `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedBinaryIsometricOp` clas.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBinaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > left, std::shared_ptr<const Matrix<Value_, Index_> > right, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBinaryIsometricOp<Value_, Index_, Op_>(std::move(left), std::move(right), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBinaryIsometricOp(std::shared_ptr<Matrix<Value_, Index_> > left, std::shared_ptr<Matrix<Value_, Index_> > right, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBinaryIsometricOp<Value_, Index_, Op_>(std::move(left), std::move(right), std::move(op)));
}
/**
 * @endcond
 */

}

#include "arith_helpers.hpp"

#include "compare_helpers.hpp"

#include "boolean_helpers.hpp"

#include "mock_helpers.hpp"

#endif
