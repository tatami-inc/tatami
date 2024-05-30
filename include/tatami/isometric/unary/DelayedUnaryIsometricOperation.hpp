#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OPERATION_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OPERATION_H

#include "../../base/Matrix.hpp"
#include "../../utils/copy.hpp"
#include "../../utils/new_extractor.hpp"
#include "../depends_utils.hpp"

#include <memory>
#include <algorithm>
#include <vector>
#include <type_traits>

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
template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseBasicFull : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseBasicFull(
        const Matrix<InputValue_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_extent(row ? matrix->ncol() : matrix->nrow()),
        my_ext(new_extractor<false, oracle_>(matrix, row, std::move(oracle), opt))
    {
        if constexpr(!same_value) {
            my_holding_buffer.resize(my_extent);
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_extent;

    std::unique_ptr<DenseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_holding_buffer;

public:
    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        if constexpr(same_value) {
            auto ptr = my_ext->fetch(i, buffer);
            copy_n(ptr, my_extent, buffer);
            my_operation.dense(my_row, my_oracle.get(i), static_cast<Index_>(0), my_extent, buffer, buffer);
        } else {
            auto ptr = my_ext->fetch(i, my_holding_buffer.data());
            my_operation.dense(my_row, my_oracle.get(i), static_cast<Index_>(0), my_extent, ptr, buffer);
        }
        return buffer;
    }
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseBasicBlock : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseBasicBlock(
        const Matrix<InputValue_, Index_>* matrix, 
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
    {
        if constexpr(!same_value) {
            my_holding_buffer.resize(my_block_length);
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_block_start, my_block_length;

    std::unique_ptr<DenseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_holding_buffer;

public:
    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        if constexpr(same_value) {
            auto ptr = my_ext->fetch(i, buffer);
            copy_n(ptr, my_block_length, buffer);
            my_operation.dense(my_row, my_oracle.get(i), my_block_start, my_block_length, buffer, buffer);
        } else {
            auto ptr = my_ext->fetch(i, my_holding_buffer.data());
            my_operation.dense(my_row, my_oracle.get(i), my_block_start, my_block_length, ptr, buffer);
        }
        return buffer;
    }
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseBasicIndex : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseBasicIndex(
        const Matrix<InputValue_, Index_>* matrix, 
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
    {
        if constexpr(!same_value) {
            my_holding_buffer.resize(my_indices_ptr->size());
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    VectorPtr<Index_> my_indices_ptr;

    std::unique_ptr<DenseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_holding_buffer;

public:
    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        const auto& indices = *my_indices_ptr;
        if constexpr(same_value) {
            auto ptr = my_ext->fetch(i, buffer);
            copy_n(ptr, indices.size(), buffer);
            my_operation.dense(my_row, my_oracle.get(i), indices, buffer, buffer);
        } else {
            auto ptr = my_ext->fetch(i, my_holding_buffer.data());
            my_operation.dense(my_row, my_oracle.get(i), indices, ptr, buffer);
        }
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
template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DenseExpandedFull : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseExpandedFull(
        const Matrix<InputValue_, Index_>* matrix, 
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

        if constexpr(!same_value) {
            my_result_vbuffer.resize(my_extent);
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_extent;

    std::vector<InputValue_> my_vbuffer;
    std::vector<Index_> my_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<OutputValue_>, bool>::type my_result_vbuffer;

public:
    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto vbuffer = my_vbuffer.data();
        auto range = my_ext->fetch(i, vbuffer, my_ibuffer.data());

        i = my_oracle.get(i);
        if constexpr(same_value) {
            copy_n(range.value, range.number, vbuffer);
            my_operation.sparse(my_row, i, range.number, vbuffer, range.index, vbuffer);
        } else {
            my_operation.sparse(my_row, i, range.number, range.value, range.index, my_result_vbuffer.data());
        }

        // avoid calling fill() if possible, as this might throw zero-related errors with non-IEEE-float types.
        if (range.number < my_extent) { 
            std::fill_n(buffer, my_extent, my_operation.template fill<OutputValue_, InputValue_>(my_row, i));
        }

        if constexpr(same_value) {
            for (Index_ i = 0; i < range.number; ++i) {
                buffer[range.index[i]] = vbuffer[i];
            }
        } else {
            for (Index_ i = 0; i < range.number; ++i) {
                buffer[range.index[i]] = my_result_vbuffer[i];
            }
        }

        return buffer;
    }
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_> 
class DenseExpandedBlock : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseExpandedBlock(
        const Matrix<InputValue_, Index_>* matrix, 
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

        if constexpr(!same_value) {
            my_result_vbuffer.resize(my_block_length);
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_block_start, my_block_length;

    std::vector<InputValue_> my_vbuffer;
    std::vector<Index_> my_ibuffer;
    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<OutputValue_>, bool>::type my_result_vbuffer;

public:
    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto vbuffer = my_vbuffer.data();
        auto range = my_ext->fetch(i, vbuffer, my_ibuffer.data());

        i = my_oracle.get(i);
        if constexpr(same_value) {
            copy_n(range.value, range.number, vbuffer);
            my_operation.sparse(my_row, i, range.number, vbuffer, range.index, vbuffer);
        } else {
            my_operation.sparse(my_row, i, range.number, range.value, range.index, my_result_vbuffer.data());
        }

        // avoid calling fill() if possible, as this might throw zero-related errors with non-IEEE-float types.
        if (range.number < my_block_length) { 
            std::fill_n(buffer, my_block_length, my_operation.template fill<OutputValue_, InputValue_>(my_row, i));
        }

        if constexpr(same_value) {
            for (Index_ i = 0; i < range.number; ++i) {
                buffer[range.index[i] - my_block_start] = vbuffer[i];
            }
        } else {
            for (Index_ i = 0; i < range.number; ++i) {
                buffer[range.index[i] - my_block_start] = my_result_vbuffer[i];
            }
        }

        return buffer;
    }
};

template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_> 
class DenseExpandedIndex : public DenseExtractor<oracle_, OutputValue_, Index_> {
public:
    DenseExpandedIndex(
        const Matrix<InputValue_, Index_>* matrix, 
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

        if constexpr(!same_value) {
            my_result_vbuffer.resize(my_extent);
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;
    Index_ my_extent;

    std::vector<InputValue_> my_vbuffer;
    std::vector<Index_> my_ibuffer;

    std::vector<Index_> my_remapping;
    Index_ my_remapping_offset = 0;
    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<OutputValue_>, bool>::type my_result_vbuffer;

public:
    const OutputValue_* fetch(Index_ i, OutputValue_* buffer) {
        auto vbuffer = my_vbuffer.data();
        auto range = my_ext->fetch(i, vbuffer, my_ibuffer.data());

        i = my_oracle.get(i);
        if constexpr(same_value) {
            copy_n(range.value, range.number, vbuffer);
            my_operation.sparse(my_row, i, range.number, vbuffer, range.index, vbuffer);
        } else {
            my_operation.sparse(my_row, i, range.number, range.value, range.index, my_result_vbuffer.data());
        }

        // avoid calling fill() if possible, as this might throw zero-related errors with non-IEEE-float types.
        if (range.number < my_extent) { 
            std::fill_n(buffer, my_extent, my_operation.template fill<OutputValue_, InputValue_>(my_row, i));
        }

        if constexpr(same_value) {
            for (Index_ i = 0; i < range.number; ++i) {
                buffer[my_remapping[range.index[i] - my_remapping_offset]] = vbuffer[i];
            }
        } else {
            for (Index_ i = 0; i < range.number; ++i) {
                buffer[my_remapping[range.index[i] - my_remapping_offset]] = my_result_vbuffer[i];
            }
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
template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class SparseSimple : public SparseExtractor<oracle_, OutputValue_, Index_> {
public:
    SparseSimple(
        const Matrix<InputValue_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row),
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), opt)) 
    {
        initialize(opt, row ? matrix->ncol() : matrix->nrow());
    }

    SparseSimple(
        const Matrix<InputValue_, Index_>* matrix, 
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
    {
        initialize(opt, block_length);
    }

    SparseSimple(
        const Matrix<InputValue_, Index_>* matrix, 
        const Operation_& operation, 
        bool row, 
        MaybeOracle<oracle_, Index_> oracle, 
        VectorPtr<Index_> indices_ptr,
        const Options& opt) :
        my_operation(operation),
        my_row(row),
        my_oracle(oracle, my_operation, row)
    {
        initialize(opt, indices_ptr->size()); 

        // Need to construct this here so that indices_ptr isn't moved until after calling initialize().
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt);
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_holding_vbuffer;

    void initialize(const Options& opt, size_t extent) {
        if constexpr(!same_value) {
            if (opt.sparse_extract_value) {
                my_holding_vbuffer.resize(extent);
            }
        }
    }

public:
    SparseRange<OutputValue_, Index_> fetch(Index_ i, OutputValue_* value_buffer, Index_* index_buffer) {
        if constexpr(same_value) {
            auto raw = my_ext->fetch(i, value_buffer, index_buffer);
            if (raw.value) {
                copy_n(raw.value, raw.number, value_buffer);
                my_operation.sparse(my_row, my_oracle.get(i), raw.number, value_buffer, raw.index, value_buffer);
                raw.value = value_buffer;
            }
            return raw;
        } else {
            auto raw = my_ext->fetch(i, my_holding_vbuffer.data(), index_buffer);

            SparseRange<OutputValue_, Index_> output(raw.number);
            output.index = raw.index;
            if (raw.value) {
                my_operation.sparse(my_row, my_oracle.get(i), raw.number, raw.value, raw.index, value_buffer);
                output.value = value_buffer;
            }

            return output;
        }
    }
};

/**
 * SparseNeedsIndices is used if:
 *
 * - the underlying matrix is sparse
 * - the operation preserves sparsity
 * - indices are necessary to perform the operation 
 */
template<bool oracle_, typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class SparseNeedsIndices : public SparseExtractor<oracle_, OutputValue_, Index_> {
public:
    SparseNeedsIndices(
        const Matrix<InputValue_, Index_>* matrix,
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
        const Matrix<InputValue_, Index_>* matrix,
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
        const Matrix<InputValue_, Index_>* matrix,
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

        if constexpr(!same_value) {
            if (my_report_value) {
                my_holding_vbuffer.resize(extent);
            }
        }
    }

private:
    const Operation_& my_operation;

    bool my_row;
    DelayedIsometricOperation_internal::MaybeOracleDepends<oracle_, Operation_, Index_> my_oracle;

    bool my_report_value, my_report_index;
    std::vector<Index_> my_ibuffer;

    std::unique_ptr<SparseExtractor<oracle_, InputValue_, Index_> > my_ext;

    static constexpr bool same_value = std::is_same<OutputValue_, InputValue_>::value;
    typename std::conditional<!same_value, std::vector<InputValue_>, bool>::type my_holding_vbuffer;

public:
    SparseRange<OutputValue_, Index_> fetch(Index_ i, OutputValue_* value_buffer, Index_* index_buffer) {
        auto iptr = my_report_index ? index_buffer : my_ibuffer.data();

        if constexpr(same_value) {
            auto raw = my_ext->fetch(i, value_buffer, iptr);
            if (my_report_value) {
                copy_n(raw.value, raw.number, value_buffer);
                my_operation.sparse(my_row, my_oracle.get(i), raw.number, value_buffer, raw.index, value_buffer);
                raw.value = value_buffer;
            }
            if (!my_report_index) {
                raw.index = NULL;
            } 
            return raw;

        } else {
            auto raw = my_ext->fetch(i, my_holding_vbuffer.data(), iptr);
            SparseRange<OutputValue_, Index_> output(raw.number, NULL, (my_report_index ? raw.index : NULL));
            if (my_report_value) {
                my_operation.sparse(my_row, my_oracle.get(i), raw.number, raw.value, raw.index, value_buffer);
                output.value = value_buffer;
            }
            return output;
        }
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
 * @tparam OutputValue_ Type of the result of the operation.
 * This is the type of the value of the output matrix.
 * @tparam InputValue_ Type of the value of the input matrix, to use in the operation.
 * This may or may not be the same as `OutputValue_`, depending on the methods available in `Operation_`.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class representing the operation of interest.
 * This should implement the same methods as `DelayedUnaryIsometricMockBasic` or `DelayedUnaryIsometricMockAdvanced`,
 * depending on whether it can take advantage of matrix sparsity.
 */
template<typename OutputValue_, typename InputValue_, typename Index_, class Operation_>
class DelayedUnaryIsometricOperation : public Matrix<OutputValue_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying matrix.
     * @param operation Instance of the functor class.
     */
    DelayedUnaryIsometricOperation(std::shared_ptr<const Matrix<InputValue_, Index_> > matrix, Operation_ operation) : 
        my_matrix(std::move(matrix)), my_operation(std::move(operation)) {}

private:
    std::shared_ptr<const Matrix<InputValue_, Index_> > my_matrix;
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

    using Matrix<OutputValue_, Index_>::dense;

    using Matrix<OutputValue_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_basic_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseBasicFull<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_matrix.get(), my_operation, row, std::move(oracle), opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_basic_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseBasicBlock<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_matrix.get(), my_operation, row, std::move(oracle), block_start, block_length, opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_basic_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseBasicIndex<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_matrix.get(), my_operation, row, std::move(oracle), std::move(indices_ptr), opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseExpandedFull<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_matrix.get(), my_operation, row, std::move(oracle), opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseExpandedBlock<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_matrix.get(), my_operation, row, std::move(oracle), block_start, block_length, opt
        );
    }

    template<bool oracle_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_expanded_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<DelayedUnaryIsometricOperation_internal::DenseExpandedIndex<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
            my_matrix.get(), my_operation, row, std::move(oracle), std::move(indices_ptr), opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, OutputValue_, Index_> > dense_internal(bool row, Args_&& ... args) const {
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
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) const {
        return std::make_unique<FullSparsifiedWrapper<oracle_, OutputValue_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), opt),
            (row ? my_matrix->ncol() : my_matrix->nrow()),
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return std::make_unique<BlockSparsifiedWrapper<oracle_, OutputValue_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), block_start, block_length, opt),
            block_start,
            block_length,
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_to_dense_internal(bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<IndexSparsifiedWrapper<oracle_, OutputValue_, Index_> >(
            dense_internal<oracle_>(row, std::move(oracle), indices_ptr, opt),
            indices_ptr,
            opt
        );
    }

    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, OutputValue_, Index_> > sparse_internal(bool row, MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) const {
        if constexpr(!Operation_::is_basic) {
            if (my_operation.is_sparse() && my_matrix->is_sparse()) { 
                if (DelayedIsometricOperation_internal::needs_sparse_indices(my_operation, row)) {
                    return std::make_unique<DelayedUnaryIsometricOperation_internal::SparseNeedsIndices<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
                        my_matrix.get(),
                        my_operation,
                        row, 
                        std::move(oracle),
                        std::forward<Args_>(args)...
                    );

                } else {
                    return std::make_unique<DelayedUnaryIsometricOperation_internal::SparseSimple<oracle_, OutputValue_, InputValue_, Index_, Operation_> >(
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
 * @tparam OutputValue_ Type of matrix value after the operation.
 * @tparam InputValue_ Type of matrix value before the operation.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class implementing the operation.
 *
 * @param matrix Pointer to a (possibly `const`) `Matrix`.
 * @param operation Instance of the operation helper class.
 *
 * @return Instance of a `DelayedUnaryIsometricOperation` class.
 */
template<typename OutputValue_ = double, typename InputValue_, typename Index_, class Operation_>
std::shared_ptr<Matrix<OutputValue_, Index_> > make_DelayedUnaryIsometricOperation(std::shared_ptr<const Matrix<InputValue_, Index_> > matrix, Operation_ operation) {
    typedef typename std::remove_reference<Operation_>::type Operation2_;
    return std::shared_ptr<Matrix<OutputValue_, Index_> >(new DelayedUnaryIsometricOperation<OutputValue_, InputValue_, Index_, Operation2_>(std::move(matrix), std::move(operation)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename OutputValue_ = double, typename InputValue_, typename Index_, class Operation_> 
std::shared_ptr<Matrix<OutputValue_, Index_> > make_DelayedUnaryIsometricOperation(std::shared_ptr<Matrix<InputValue_, Index_> > matrix, Operation_ operation) {
    typedef typename std::remove_reference<Operation_>::type Operation2_;
    return std::shared_ptr<Matrix<OutputValue_, Index_> >(new DelayedUnaryIsometricOperation<OutputValue_, InputValue_, Index_, Operation2_>(std::move(matrix), std::move(operation)));
}

// For back-compatibility.
template<typename ... Args_>
auto make_DelayedIsometricOperation(Args_&&... args) {
    return make_DelayedUnaryIsometricOperation(std::forward<Args_>(args)...);
}

template<typename OutputValue_, typename Index_, class Operation_>
using DelayedIsometricOperation = DelayedUnaryIsometricOperation<OutputValue_, Index_, Operation_, OutputValue_>;
/**
 * @endcond
 */

}

#endif
