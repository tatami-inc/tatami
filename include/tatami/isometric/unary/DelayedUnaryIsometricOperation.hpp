#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OPERATION_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OPERATION_H

#include "../../base/Matrix.hpp"
#include "../../utils/copy.hpp"
#include "../../utils/new_extractor.hpp"
#include "../depends_utils.hpp"

#include <memory>
#include <algorithm>
#include <vector>

/**
 * @file DelayedUnaryIsometricOperation.hpp
 *
 * @brief Delayed unary isometric operations.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedUnaryIsometricOperation_internal {

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
class DenseBasicFull : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseBasicFull(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_extent(row ? matrix->ncol() : matrix->nrow()),
        my_ext(new_extractor<false, oracle_>(matrix, row, std::move(oracle), opt))
    {}

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_extent;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto ptr = my_ext->fetch(i, buffer);
        copy_n(ptr, my_extent, buffer);
        my_operation.dense(my_row, my_oracle.get(i), static_cast<Index_>(0), my_extent, buffer);
        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
class DenseBasicBlock : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseBasicBlock(
        const Matrix<Value_, Index_>* matrix, 
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
        my_block_length(block_length),
        my_ext(new_extractor<false, oracle_>(matrix, row, std::move(oracle), block_start, block_length, opt))
    {}

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_block_start, my_block_length;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto ptr = my_ext->fetch(i, buffer);
        copy_n(ptr, my_block_length, buffer);
        my_operation.dense(my_row, my_oracle.get(i), my_block_start, my_block_length, buffer);
        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_>
class DenseBasicIndex : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseBasicIndex(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> indices_ptr,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_indices_ptr(indices_ptr),
        my_ext(new_extractor<false, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt))
    {}

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    VectorPtr<Index_> my_indices_ptr;
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto ptr = my_ext->fetch(i, buffer);
        const auto& indices = *my_indices_ptr;
        copy_n(ptr, indices.size(), buffer);
        my_operation.dense(my_row, my_oracle.get(i), indices, buffer);
        return buffer;
    }
};

/**
 * DenseExpanded is used if:
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
class DenseExpandedFull : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseExpandedFull(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation,
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_extent(row ? matrix->ncol() : matrix->nrow()),
        my_vbuffer(my_extent),
        my_ibuffer(my_extent)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        my_ext = new_extractor<true, oracle_>(matrix, my_row, std::move(oracle), opt);
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    Index_ my_extent;
    std::vector<Value_> my_vbuffer;
    std::vector<Index_> my_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto vbuffer = my_vbuffer.data();
        auto range = my_ext->fetch(i, vbuffer, my_ibuffer.data());
        copy_n(range.value, range.number, vbuffer);

        i = my_oracle.get(i);
        my_operation.sparse(my_row, i, range.number, vbuffer, range.index);

        // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
        if (range.number < my_extent) { 
            std::fill_n(buffer, my_extent, my_operation.template fill<Value_>(my_row, i));
        }

        for (Index_ i = 0; i < range.number; ++i) {
            buffer[range.index[i]] = vbuffer[i];
        }

        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_> 
class DenseExpandedBlock : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseExpandedBlock(
        const Matrix<Value_, Index_>* matrix, 
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
        my_block_length(block_length),
        my_vbuffer(block_length),
        my_ibuffer(block_length)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), block_start, block_length, opt);
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    Index_ my_block_start, my_block_length;
    std::vector<Value_> my_vbuffer;
    std::vector<Index_> my_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto vbuffer = my_vbuffer.data();
        auto range = my_ext->fetch(i, vbuffer, my_ibuffer.data());
        copy_n(range.value, range.number, vbuffer);

        i = my_oracle.get(i);
        my_operation.sparse(my_row, i, range.number, vbuffer, range.index);

        // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
        if (range.number < my_block_length) { 
            std::fill_n(buffer, my_block_length, my_operation.template fill<Value_>(my_row, i));
        }

        for (Index_ i = 0; i < range.number; ++i) {
            buffer[range.index[i] - my_block_start] = vbuffer[i];
        }

        return buffer;
    }
};

template<bool oracle_, typename Value_, typename Index_, class Operation_> 
class DenseExpandedIndex : public DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseExpandedIndex(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation,
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        opt.sparse_extract_value = true;
        opt.sparse_extract_index = true;

        const auto& indices = *indices_ptr;
        my_extent = indices.size();
        my_vbuffer.resize(my_extent);
        my_ibuffer.resize(my_extent);

        // Create a remapping vector to map the extracted indices back to the
        // dense buffer. We use the 'remapping_offset' to avoid allocating the
        // full extent of the dimension.
        if (my_extent) {
            my_remapping_offset = indices.front();
            my_remapping.resize(indices.back() - my_remapping_offset + 1);
            for (Index_ i = 0; i < my_extent; ++i) {
                my_remapping[indices[i] - my_remapping_offset] = i;
            }
        }

        my_ext = new_extractor<true, oracle_>(matrix, my_row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    Index_ my_extent;
    std::vector<Value_> my_vbuffer;
    std::vector<Index_> my_ibuffer;

    std::vector<Index_> my_remapping;
    Index_ my_remapping_offset = 0;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto vbuffer = my_vbuffer.data();
        auto range = my_ext->fetch(i, vbuffer, my_ibuffer.data());
        copy_n(range.value, range.number, vbuffer);

        i = my_oracle.get(i);
        my_operation.sparse(my_row, i, range.number, vbuffer, range.index);

        // avoid calling zero() if possible, as this might throw zero-related errors in non-IEEE platforms.
        if (range.number < my_extent) { 
            std::fill_n(buffer, my_extent, my_operation.template fill<Value_>(my_row, i));
        }

        for (Index_ i = 0; i < range.number; ++i) {
            buffer[my_remapping[range.index[i] - my_remapping_offset]] = vbuffer[i];
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
class SparseSimple : public SparseExtractor<oracle_, Value_, Index_> {
public:
    SparseSimple(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), opt)) 
    {}

    SparseSimple(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start,
        Index_ block_length,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), block_start, block_length, opt)) 
    {}

    SparseSimple(
        const Matrix<Value_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> indices_ptr,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt)) 
    {}

private:
    const Operation_& my_operation;
    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto raw = my_ext->fetch(i, value_buffer, index_buffer);
        if (raw.value) {
            copy_n(raw.value, raw.number, value_buffer);
            my_operation.sparse(my_row, my_oracle.get(i), raw.number, value_buffer, raw.index);
            raw.value = value_buffer;
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
template<bool oracle_, typename Value_, typename Index_, class Operation_>
class SparseNeedsIndices : public SparseExtractor<oracle_, Value_, Index_> {
public:
    SparseNeedsIndices(
        const Matrix<Value_, Index_>* matrix,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        initialize(opt, row ? matrix->ncol() : matrix->nrow());
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), opt);
    }

    SparseNeedsIndices(
        const Matrix<Value_, Index_>* matrix,
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
        initialize(opt, block_length);
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), block_start, block_length, opt);
    }

    SparseNeedsIndices(
        const Matrix<Value_, Index_>* matrix,
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> indices_ptr,
        Options opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        initialize(opt, indices_ptr->size());
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    void initialize(Options& opt, size_t extent) {
        my_report_value = opt.sparse_extract_value;
        my_report_index = opt.sparse_extract_index;

        // The index is not required if we don't even want the values,
        // in which case Operation_::is_sparse() isn't even called.
        if (my_report_value) {
            opt.sparse_extract_index = true;

            // We only need an internal ibuffer if the user wants the
            // values but didn't provide enough space to store the indices
            // (which we need to pass to Operation_::is_sparse()).
            if (!my_report_index) {
                my_ibuffer.resize(extent);
            }
        }
    }

private:
    const Operation_& my_operation;
    bool my_row;
    bool my_report_value, my_report_index;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    std::vector<Index_> my_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto raw = my_ext->fetch(i, value_buffer, my_report_index ? index_buffer : my_ibuffer.data());

        if (my_report_value) {
            copy_n(raw.value, raw.number, value_buffer);
            my_operation.sparse(my_row, my_oracle.get(i), raw.number, value_buffer, raw.index);
            raw.value = value_buffer;
        }

        if (!my_report_index) {
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
 * @brief Delayed isometric operation on a single matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated during data extraction, e.g., with `MyopicDenseExtractor::fetch()` or friends.
 * We only consider "unary" operations that involve a single `Matrix` - see `DelayedBinaryIsometricOperation` for operations between two `Matrix` instances.
 *
 * This class is inspired by the `DelayedUnaryIsoOp` classes from the **DelayedArray** Bioconductor package.
 * 
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class implementing the operation.
 * This should implement the same methods as `DelayedUnaryIsometricMockBasic` or `DelayedUnaryIsometricMockAdvanced`,
 * depending on whether it can take advantage of matrix sparsity.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedUnaryIsometricOperation : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying matrix.
     * @param operation Instance of the functor class.
     */
    DelayedUnaryIsometricOperation(std::shared_ptr<const Matrix<Value_, Index_> > matrix, Operation_ operation) : my_matrix(std::move(matrix)), my_operation(std::move(operation)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > my_matrix;
    Operation_ my_operation;

public:
    Index_ nrow() const {
        return my_matrix->nrow();
    }
    
    Index_ ncol() const {
        return my_matrix->ncol();
    }

    bool is_sparse() const {
        if constexpr(!Operation_::is_basic) {
            if (my_operation.is_sparse()) {
                return my_matrix->is_sparse();
            }
        }
        return false;
    }

    double is_sparse_proportion() const {
        if constexpr(!Operation_::is_basic) {
            if (my_operation.is_sparse()) {
                return my_matrix->is_sparse_proportion();
            }
        }
        return 0;
    }

    bool prefer_rows() const { 
        return my_matrix->prefer_rows();
    }

    double prefer_rows_proportion() const { 
        return my_matrix->prefer_rows_proportion();
    }

    bool uses_oracle(bool row) const {
        return my_matrix->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense;

    using Matrix<Value_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_basic_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseBasicFull<oracle_, Value_, Index_, Operation_> >(my_matrix.get(), my_operation, row, std::move(oracle), opt);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_basic_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseBasicBlock<oracle_, Value_, Index_, Operation_> >(my_matrix.get(), my_operation, row, std::move(oracle), block_start, block_length, opt);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_basic_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseBasicIndex<oracle_, Value_, Index_, Operation_> >(my_matrix.get(), my_operation, row, std::move(oracle), std::move(indices_ptr), opt);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseExpandedFull<oracle_, Value_, Index_, Operation_> >(my_matrix.get(), my_operation, row, std::move(oracle), opt);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseExpandedBlock<oracle_, Value_, Index_, Operation_> >(my_matrix.get(), my_operation, row, std::move(oracle), block_start, block_length, opt);
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseExpandedIndex<oracle_, Value_, Index_, Operation_> >(my_matrix.get(), my_operation, row, std::move(oracle), std::move(indices_ptr), opt);
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, Args_&& ... args) const {
        if constexpr(!Operation_::is_basic) {
            if (my_matrix->is_sparse()) {
                if (DelayedIsometricOperation_internal::can_dense_expand(my_operation, row)) {
                    return dense_expanded_internal<oracle_>(row, std::forward<Args_>(args)...);
                }
            }
        }

        return dense_basic_internal<oracle_>(row, std::forward<Args_>(args)...);
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
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<FullSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), opt),
            (row ? my_matrix->ncol() : my_matrix->nrow()),
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<BlockSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), block_start, block_length, opt),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<IndexSparsifiedWrapper<oracle_, Value_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), indices_ptr, opt),
            indices_ptr,
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) const {
        if constexpr(!Operation_::is_basic) {
            if (my_operation.is_sparse() && my_matrix->is_sparse()) { 
                if (DelayedIsometricOperation_internal::needs_sparse_indices(my_operation, row)) {
                    return std::make_unique<DelayedUnaryIsometricOperation_internal::SparseNeedsIndices<oracle_, Value_, Index_, Operation_> >(
                        my_matrix.get(),
                        my_operation,
                        row, 
                        std::move(oracle),
                        std::forward<Args_>(args)...
                    );

                } else {
                    return std::make_unique<DelayedUnaryIsometricOperation_internal::SparseSimple<oracle_, Value_, Index_, Operation_> >(
                        my_matrix.get(),
                        my_operation, 
                        row, 
                        std::move(oracle), 
                        std::forward<Args_>(args)...
                    );
                }
            }
        }

        return sparse_to_dense_internal<oracle_>(row, std::move(oracle), std::forward<Args_>(args)...);
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
 * @tparam Operation_ Helper class implementing the operation.
 *
 * @param matrix Pointer to a (possibly `const`) `Matrix`.
 * @param operation Instance of the operation helper class.
 *
 * @return Instance of a `DelayedUnaryIsometricOperation` class.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOperation(std::shared_ptr<const Matrix<Value_, Index_> > matrix, Operation_ operation) {
    typedef typename std::remove_reference<Operation_>::type Operation2_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOperation<Value_, Index_, Operation2_>(std::move(matrix), std::move(operation)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOperation(std::shared_ptr<Matrix<Value_, Index_> > matrix, Operation_ operation) {
    typedef typename std::remove_reference<Operation_>::type Operation2_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOperation<Value_, Index_, Operation2_>(std::move(matrix), std::move(operation)));
}

// For back-compatibility.
template<typename ... Args_>
auto make_DelayedIsometricOperation(Args_&&... args) {
    return make_DelayedUnaryIsometricOperation(std::forward<Args_>(args)...);
}

template<typename Value_, typename Index_, class Operation_>
using DelayedIsometricOperation = DelayedUnaryIsometricOperation<Value_, Index_, Operation_>;
/**
 * @endcond
 */

}

#include "arithmetic_helpers.hpp"

#include "math_helpers.hpp"

#include "compare_helpers.hpp"

#include "boolean_helpers.hpp"

#include "mock_helpers.hpp"

#endif
