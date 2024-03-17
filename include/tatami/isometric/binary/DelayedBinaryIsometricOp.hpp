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
namespace DelayedIsometricBinaryOp_internal {

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseFull : public DenseExtractor<oracle_, Value_, Index_> {
    DenseFull(
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
        operation.template dense<Value_, Index_>(row, i, 0, extent, buffer, rptr);
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
struct DenseBlock : public DenseExtractor<oracle_, Value_, Index_> {
    DenseBlock(
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
        operation.template dense<Value_, Index_>(row, i, block_start, block_length, buffer, rptr);
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
struct DenseIndex : public DenseExtractor<oracle_, Value_, Index_> {
    DenseIndex(
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
        operation.template dense<Value_, Index_>(row, i, *indices_ptr, buffer, rptr);
        return buffer;
    }

private:
    const Operation_& operation;
    bool row;
    VectorPtr<Index_> indices_ptr;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > left, right;
    std::vector<Value_> holding_buffer;
};

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
        return operation.template sparse<Value_, Index_>(
            row, 
            i, 
            left_ranges, 
            right_ranges, 
            vbuffer,
            ibuffer,
            report_value,
            report_index
        );
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
 * The `Operation_` class is expected to provide the following static `constexpr` member variables:
 *
 * - `always_sparse`: whether the operation can be optimized to return a sparse result if both input matrices are sparse.
 * 
 * The class should implement the following method:
 *
 * - `void dense(bool row, Index_ i, Index_ start, Index_ length, Value_* left_buffer, const Value_* right_buffer) const`: 
 *   This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
 *   When `row = true`, these buffers contain a contiguous block of elements from row `i` of the left and right matrices, 
 *   where the block starts at column `start` and is of length `length`.
 *   If `row = false`, `i` is instead a column and the block starts at row `start`.
 *   The result of the operation should be stored in `left_buffer`.
 * - `void dense(bool row, Index_ i, const std::vector<Index_>& indices, Value_* buffer1, const Value_* buffer2) const`: 
 *   This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`.
 *   When `row = true`, these buffers contain a subset of elements from row `i` of the left and right matrices, 
 *   where the subset is defined by column indices in the `indices` array of length `length`.
 *   If `row_ = false`, `i` is instead a column and `indices` contains row indices.
 *   The result of the operation should be stored in `left_buffer`.
 * 
 * If `always_sparse = true`, the class should implement:
 *
 * - `Index_ sparse(bool row, Index_ i, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer, bool needs_value, bool needs_index) const`:
 *   This method should apply the operation to the sparse ranges in `left` and `right`.
 *   These ranges consist of the contents of row `i` (if `row = true)` or column `i` (otherwise) from the left and right matrices.
 *   `left` and `right` will always return indices regardless of `needs_index`, and these are always sorted by increasing index;
 *   however, values will only be returned if `needs_value = true`.
 *   All non-zero values resulting from the operation should be stored in `value_buffer` if `needs_value = true`, otherwise `value_buffer = NULL` and should be ignored.
 *   The corresponding indices of those values should be stored in `index_buffer` if `needs_index = true`, otherwise `index_buffer = NULL` and should be ignored.
 *   The return value should be the number of structural non-zero elements in the output buffers.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
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
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > left, right;
    Operation_ operation;
    double prefer_rows_proportion_internal;

public:
    Index_ nrow() const {
        return left->nrow();
    }

    Index_ ncol() const {
        return left->ncol();
    }

    /**
     * @return `true` if both underlying (pre-operation) matrices are sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(Operation_::always_sparse) {
            return left->sparse() && right->sparse();
        }
        return false;
    }

    double sparse_proportion() const {
        if constexpr(Operation_::always_sparse) {
            // Well, better than nothing.
            return (left->sparse_proportion() + right->sparse_proportion())/2;
        }
        return 0;
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
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedIsometricBinaryOp_internal::DenseFull<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedIsometricBinaryOp_internal::DenseBlock<oracle_, Value_, Index_, Operation_> >(
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
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedIsometricBinaryOp_internal::DenseIndex<oracle_, Value_, Index_, Operation_> >(
            left.get(),
            right.get(),
            operation,
            row, 
            std::move(oracle),
            std::move(indices_ptr),
            opt
        );
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
        if constexpr(Operation_::always_sparse) {
            return std::make_unique<DelayedIsometricBinaryOp_internal::Sparse<oracle_, Value_, Index_, Operation_> >(
                left.get(),
                right.get(),
                operation,
                row, 
                std::move(oracle),
                opt
            );
        } else {
            return std::make_unique<FullSparsifiedWrapper<oracle_, Value_, Index_> >(
                dense_internal(row, std::move(oracle), opt),
                row ? left->ncol() : left->nrow(),
                opt
            );
        }
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        if constexpr(Operation_::always_sparse) {
            return std::make_unique<DelayedIsometricBinaryOp_internal::Sparse<oracle_, Value_, Index_, Operation_> >(
                left.get(),
                right.get(),
                operation,
                row, 
                std::move(oracle),
                block_start,
                block_length,
                opt
            );
        } else {
            return std::make_unique<BlockSparsifiedWrapper<oracle_, Value_, Index_> >(
                dense_internal(row, std::move(oracle), block_start, block_length, opt),
                block_start,
                block_length,
                opt
            );
        }
    }

    template<bool oracle_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if constexpr(Operation_::always_sparse) {
            return std::make_unique<DelayedIsometricBinaryOp_internal::Sparse<oracle_, Value_, Index_, Operation_> >(
                left.get(),
                right.get(),
                operation,
                row, 
                std::move(oracle),
                std::move(indices_ptr),
                opt
            );
        } else {
            return std::make_unique<IndexSparsifiedWrapper<oracle_, Value_, Index_> >(
                dense_internal(row, std::move(oracle), indices_ptr, opt),
                indices_ptr,
                opt
            );
        }
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

#endif
