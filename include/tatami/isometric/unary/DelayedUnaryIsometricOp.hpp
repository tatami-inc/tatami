#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H

#include "../../base/Matrix.hpp"
#include "../../utils/copy.hpp"
#include "../../utils/new_extractor.hpp"

#include <memory>
#include <algorithm>
#include <vector>

/**
 * @file DelayedUnaryIsometricOp.hpp
 *
 * @brief Delayed unary isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedUnaryIsometricOp_internal {

/**
 * DenseBasic is used if:
 *
 * - the underlying matrix is dense.
 *
 * OR
 *
 * - the underlying matrix is sparse
 * - the operation discards sparsity in a variable manner.
 */
template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseBasicFull : public DenseExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    DenseBasicFull(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        const Options& opt) :
        operation(op),
        row(row),
        extent(row ? p->ncol() : p->nrow()),
        internal(new_extractor<false, oracle_>(p, row, std::move(oracle), opt))
    {}

private:
    const Operation_& op;
    bool row;
    Index_ extent;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto ptr = internal->fetch(i, buffer);
        copy_n(ptr, extent, buffer);
        operation.dense(row, i, 0, extent, buffer);
        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseBasicBlock : public DenseExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    DenseBasicBlock(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        Index_ bs,
        Index_ bl,
        const Options& opt) :
        operation(op),
        row(row),
        block_start(bs),
        block_length(bl),
        internal(new_extractor<false, oracle_>(p, row, std::move(oracle), bs, bl, opt))
    {}

private:
    const Operation_& op;
    bool row;
    Index_ block_start, block_length;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto ptr = internal->fetch(i, buffer);
        copy_n(ptr, extent, buffer);
        operation.dense(row, i, block_start, block_length, buffer);
        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct DenseBasicIndex : public DenseExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    DenseBasicIndex(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> idx_ptr,
        const Options& opt) :
        operation(op),
        row(row),
        indices_ptr(idx_ptr),
        internal(new_extractor<false, oracle_>(p, row, std::move(oracle), std::move(idx_ptr), opt))
    {}

private:
    const Operation_& op;
    bool row;
    VectorPtr<Index_> indices_ptr;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto ptr = internal->fetch(i, buffer);
        copy_n(ptr, extent, buffer);
        operation.dense(row, i, *indices_ptr, buffer);
        return buffer;
    }
};

/**
 * DenseFromSparse is used if:
 *
 * - the underlying matrix is sparse
 * - the operation preserves sparsity
 * 
 * OR
 *
 * - the underlying matrix is sparse
 * - the operation discards sparsity in a constant manner.
 */
template<bool oracle_, typename Value_, typename Index_, class Operation_> 
struct DenseFromSparseFull : public DenseIsometricExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    DenseFromSparseFull(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op,
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Options opt) :
        operation(op),
        row(row),
        extent(row ? p->ncol() : p->nrow()),
        internal_vbuffer(extent),
        internal_ibuffer(extent)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        internal = new_extractor<true, oracle_>(p, row, std::move(oracle), opt);
    }

private:
    const Operation_& op;
    bool row;
    Index_ extent;
    std::vector<Value_> internal_vbuffer;
    std::vector<Index_> internal_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto vbuffer = internal_vbuffer.data();
        auto range = this->internal->fetch(i, vbuffer, internal_ibuffer.data());
        operation.sparse(row, i, range.number, vbuffer, range.index);

        // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
        if (range.number < extent) { 
            std::fill_n(buffer, extent, operation.zero(i));
        }

        for (Index_ i = 0; i < range.number; ++i) {
            buffer[range.index[i]] = vbuffer[i];
        }

        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_> 
struct DenseFromSparseBlock : public DenseIsometricExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    DenseFromSparseBlock(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op,
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Index_ bs,
        Index_ bl,
        Options opt) :
        operation(op),
        row(row),
        block_start(bs),
        block_length(bl),
        internal_vbuffer(block_length),
        internal_ibuffer(block_length)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        internal = new_extractor<true, oracle_>(p, row, std::move(oracle), block_start, block_length, opt);
    }

private:
    const Operation_& op;
    bool row;
    Index_ block_start, block_length;
    std::vector<Value_> internal_vbuffer;
    std::vector<Index_> internal_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto vbuffer = internal_vbuffer.data();
        auto range = this->internal->fetch(i, vbuffer, internal_ibuffer.data());
        operation.sparse(row, i, range.number, vbuffer, range.index);

        // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
        if (range.number < block_length) { 
            std::fill_n(buffer, extent, operation.zero(i));
        }

        for (Index_ i = 0; i < range.number; ++i) {
            buffer[range.index[i] - block_start] = vbuffer[i];
        }

        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_> 
struct DenseFromSparseIndex : public DenseIsometricExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    DenseFromSparseIndex(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op,
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        operation(op),
        row(row)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;

        const auto& indices = *indices_ptr;
        extent = indices.size();
        internal_vbuffer.resize(extent);
        internal_ibuffer.resize(extent);

        if (extent) {
            index_mapping.resize(row ? p->ncol() : p->nrow());
            for (Index_ i = 0; i < extent; ++i) {
                index_mapping[indices[i]] = i;
            }
        }

        internal = new_extractor<true, oracle_>(p, row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    const Operation_& op;
    bool row;
    Index_ extent;
    std::vector<Value_> internal_vbuffer;
    std::vector<Index_> internal_ibuffer;
    std::vector<Index_> index_mapping;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto vbuffer = internal_vbuffer.data();
        auto range = this->internal->fetch(i, vbuffer, internal_ibuffer.data());
        operation.sparse(row, i, range.number, vbuffer, range.index);

        // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
        if (range.number < extent) { 
            std::fill_n(buffer, extent, operation.zero(i));
        }

        for (Index_ i = 0; i < range.number; ++i) {
            buffer[index_mapping[range.index[i]]] = vbuffer[i];
        }

        return buffer;
    }
};

/**
 * SparseSimple is used if:
 *
 * - the underlying matrix is sparse
 * - the operation preserves sparsity
 * - indices are not necessary to perform the operation 
 */
template<bool oracle_, typename Value_, typename Index_, class Operation_>
struct SparseSimple : public SparseExtractor<oracle_, Value_, Index_> {
    SparseSimple(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        const Options& opt) :
        operation(op),
        row(row),
        internal(new_extractor<true, oracle_>(p, row, std::move(oracle), opt)) 
    {}

    SparseSimple(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        Index_ bs,
        Index_ bl,
        const Options& opt) :
        operation(op),
        row(row),
        internal(new_extractor<true, oracle_>(p, row, std::move(oracle), bs, bl, opt)) 
    {}

    SparseSimple(
        const Matrix<Value_, Index_>* p, 
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> indices_ptr,
        const Options& opt) :
        operation(op),
        row(row),
        internal(new_extractor<true, oracle_>(p, row, std::move(oracle), std::move(indices_ptr), opt)) 
    {}

private:
    const Operation_& op;
    bool row;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;

public:
    SparseRange<Value_, Index_> fetch_internal(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        auto raw = internal->fetch(i, vbuffer, ibuffer);
        if (raw.value) {
            copy_n(raw.value, raw.number, vbuffer);
            operation.sparse(row, i, raw.number, vbuffer, raw.index);
            raw.value = vbuffer;
        }
        return raw;
    }
};

/**
 * SparseNeedsIndices is used if:
 *
 * - the underlying matrix is sparse
 * - the operation preserves sparsity
 * - indices are necessary to perform the operation 
 */
template<bool oracle_, typename Value_, typename Index_>
struct SparseNeedsIndices : public SparseExtractor<oracle_, Value_, Index_> {
    SparseNeedsIndices(
        const Matrix<Value_, Index_>* p,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        Options opt) :
        operation(op),
        row(row)
    {
        initialize(opt, row ? p->ncol() : p->nrow());
        internal = new_extractor<true, oracle_>(p, row, std::move(oracle), opt);
    }

    SparseNeedsIndices(
        const Matrix<Value_, Index_>* p,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        Index_ bs,
        Index_ bl,
        Options opt) :
        operation(op),
        row(row)
    {
        initialize(opt, bl);
        internal = new_extractor<true, oracle_>(p, row, std::move(oracle), bs, bl, opt);
    }

    SparseNeedsIndices(
        const Matrix<Value_, Index_>* p,
        const Operation_& op, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        operation(op),
        row(row)
    {
        initialize(opt, indices_ptr->size());
        internal = new_extractor<true, oracle_>(p, row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    void initialize(Options& opt, size_t extent) {
        report_value = opt.sparse_extract_value;
        report_index = opt.sparse_extract_index;

        // The index is not required if we don't even want the values,
        // in which case Operation_::sparse() isn't even called.
        if (report_value) {
            opt.sparse_extract_index = true;

            // We only need an internal ibuffer if the user wants the
            // values but didn't provide enough space to store the indices
            // (which we need to pass to Operation_::sparse()).
            if (!report_index) {
                internal_ibuffer.resize(extent);
            }
        }
    }

private:
    const Operation_& op;
    bool row;
    bool report_value, report_index;
    std::vector<Index_> internal_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto raw = internal->fetch(i, vbuffer, report_index ? ibuffer : internal_ibuffer.data());

        if (report_value) {
            copy_n(raw.value, raw.number, vbuffer);
            operation.sparse(row, i, raw.number, vbuffer, raw.index);
            raw.value = vbuffer;
        }

        if (!report_index) {
            raw.index = NULL;
        } 

        return raw;
    }
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed isometric operations on a single matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `DenseExtractor::fetch()` or friends.
 * We only consider "unary" operations that involve a single `Matrix` - see `DelayedBinaryIsometricOp` for operations between two `Matrix` instances.
 *
 * The `Operation_` class is expected to provide the following static `constexpr` member variables:
 *
 * - `needs_row`: whether the operation needs the row index, e.g., if the operation varies by row.
 * - `needs_column`: whether the operation needs the column index, e.g., if the operation varies by column.
 * 
 * The class should implement the following methods:
 *
 * - `void dense(bool row, Index_ i, Index_ start, Index_ length, Value_* buffer) const`: 
 *   This method should apply the operation in-place to all values in `buffer`, which contains a contiguous block of columns from row `i` (if `row = true`) or a block of rows from column `i` (otherwise).
 *   This block spans indices from `[start, start + length)`.
 * - `void dense(bool row, Index_ i, const vector<Index_>& subset, Value_* buffer) const`: 
 *   This method should apply the operation in-place to all values in `buffer`, which contains a subset of columns from row `i` (if `row = true`) or a subset of rows from column `i` (otherwise).
 *   `subset` contains the sorted and unique indices of columns/rows in the subset.
 * - `void sparse(bool row, Index_ i, Index_ number, Value_* values, const Index_* indices) const`:
 *   This method should apply the operation to all values in `values`, which contains `number` structural non-zero elements from row `i` (if `row = true`) or column `i` (otherwise).
 *   The column/row index of each element is reported in `indices` - unless `needs_column` (if `row = true`) or `needs_row` (otherwise) is false, in which case `indices` may be `NULL` and should be ignored.
 *   (Note that, in this method, the operation only needs to be applied to the non-zero values in `values`, even if the operation yields a non-zero result when applied to zero.)
 * - `bool is_sparse() const`: whether this particular instance of the operation yields a sparse result when applied on sparse data of type `Value_`.
 *   For example, an addition operation remains sparse if the added value is zero.
 * - `Value_ zero(Index_ i) const`:
 *   This method should return the result of applying the operation on a zero input for row `i` (if `needs_column = false`) or column `i` (if `needs_row = false`).
 *   If both `needs_column` and `needs_row` are false, `i` should not be used as the result of the operation should be constant for all rows/columns.
 *   If `needs_column` and `needs_row` are both true, this method will never be called so any placeholder value may be returned.
 * 
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedUnaryIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedUnaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) : mat(std::move(p)), operation(std::move(op)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    Operation_ operation;
    static_assert(!Operation_::always_dense || !Operation_::always_sparse);

public:
    Index_ nrow() const {
        return mat->nrow();
    }
    
    Index_ ncol() const {
        return mat->ncol();
    }

    /**
     * @return `true` if both the underlying (pre-operation) matrix is sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        return mat->sparse() && operation.is_sparse();
    }

    double sparse_proportion() const {
        if (operation.is_sparse()) {
            return mat->sparse_proportion();
        } else {
            return 0;
        }
    }

    bool prefer_rows() const { 
        return mat->prefer_rows();
    }

    double prefer_rows_proportion() const { 
        return mat->prefer_rows_proportion();
    }

    bool uses_oracle(bool row) const {
        return mat->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, Args_&& ... args) const {
        if (!(mat->sparse())) {
            // Not sparse so we just extract it as dense.
            return std::make_unique<DelayedIsometricUnaryOp_internal::DenseBasic<oracle_, Value_, Index_> >(mat.get(), operation, row, std::forward<Args_>(args)...);
        }

        if (operation.is_sparse()) {
            // Operation preserves sparsity, no questions asked.
            return std::make_unique<DelayedIsometricUnaryOp_internal::DenseFromSparse<oracle_, Value_, Index_> >(mat.get(), operation, row, std::forward<Args_>(args)...);
        }

        if constexpr(!Operation_::needs_column) {
            if (row) {
                // Constant sparsity breaking for each row.
                return std::make_unique<DelayedIsometricUnaryOp_internal::DenseFromSparse<oracle_, Value_, Index_> >(mat.get(), operation, row, std::forward<Args_>(args)...);
            }
        }

        if constexpr(!Operation_::needs_row) {
            if (!row) {
                // Constant sparsity breaking for each column.
                return std::make_unique<DelayedIsometricUnaryOp_internal::DenseFromSparse<oracle_, Value_, Index_> >(mat.get(), operation, row, std::forward<Args_>(args)...);
            }
        }

        // Otherwise, variable sparsity breaking
        return std::make_unique<DelayedIsometricUnaryOp_internal::DenseBasic<oracle_, Value_, Index_> >(mat.get(), operation, row, std::forward<Args_>(args)...);
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
public:
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Option& opt) const {
        return std::make_unique<FullSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal(row, std::move(oracle), opt),
            (row ? mat->ncol() : mat->nrow()),
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Option& opt) const {
        return std::make_unique<BlockSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal(row, std::move(oracle), block_start, block_length, opt),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Option& opt) const {
        return std::make_unique<IndexSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal(row, std::move(oracle), indices_ptr, opt),
            indices_ptr,
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) const {
        if (!preserves_sparsity()) {
            return sparse_to_dense_internal(row, std::move(oracle), std::forward<Args_>(args)...);
        }

        if constexpr(!Operation_::needs_column) {
            if (row) {
                // no need for column indices
                return std::make_unique<DelayedIsometricUnaryOp_internal::SparseSimple<oracle_, Value_, Index_> >(mat.get(), row, std::move(oracle), std::forward<Args_>(args)...);
            }
        }

        if constexpr(!Operation_::needs_row) {
            if (!row) {
                // no need for row indices
                return std::make_unique<DelayedIsometricUnaryOp_internal::SparseSimple<oracle_, Value_, Index_> >(mat.get(), row, std::move(oracle), std::forward<Args_>(args)...);
            }
        }

        return std::make_unique<DelayedIsometricUnaryOp_internal::SparseNeedsIndices<oracle_, Value_, Index_>(mat.get(), row, std::move(oracle), std::forward<Args_>(args)...);
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
 * @param p Pointer to a (possibly `const`) `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedUnaryIsometricOp` clas.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOp(std::shared_ptr<Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

// For back-compatibility.
template<typename ... Args_>
auto make_DelayedUnaryIsometricOp(Args_&&... args) {
    return make_DelayedUnaryIsometricOp(std::forward<Args_>(args)...);
}

template<typename Value_, typename Index_, class Operation_>
using DelayedIsometricOp = DelayedUnaryIsometricOp<Value_, Index_, Operation_>;
/**
 * @endcond
 */

}

#include "arith_helpers.hpp"

#include "math_helpers.hpp"

#include "compare_helpers.hpp"

#include "boolean_helpers.hpp"

#endif
