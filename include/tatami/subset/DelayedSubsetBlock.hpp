#ifndef TATAMI_DELAYED_SUBSET_BLOCK
#define TATAMI_DELAYED_SUBSET_BLOCK

#include "../base/Matrix.hpp"
#include "../utils/new_extractor.hpp"

#include <vector>
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetBlock.hpp
 *
 * @brief Delayed subsetting to a single contiguous block.
 *
 * This is a specialized implementation that is more efficient than the `tatami::DelayedSubset` class.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedSubsetBlock_internal {

template<typename Index_>
void bump_indices(VectorPtr<Index_>& indices_ptr, const Index_ subset_start) {
    if (subset_start) {
        const auto ptr2 = new std::vector<Index_>(*indices_ptr);
        indices_ptr.reset(ptr2);
        for (auto& i : *ptr2) {
            i += subset_start;
        }
    }
}

template<bool oracle_, typename Value_, typename Index_>
class AlongDense final : public DenseExtractor<oracle_, Value_, Index_> {
public:
    AlongDense(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const Index_ subset_length,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) :
        my_ext(new_extractor<false, oracle_>(matrix, row, std::move(oracle), subset_start, subset_length, opt))
    {}

    AlongDense(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const Index_ /* for consistency */,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) :
        my_ext(new_extractor<false, oracle_>(matrix, row, std::move(oracle), subset_start + block_start, block_length, opt))
    {}

    AlongDense(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const Index_ /* for consistency */,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) {
        bump_indices(indices_ptr, subset_start); 
        my_ext = new_extractor<false, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt);
    }

    const Value_* fetch(const Index_ i, Value_* const buffer) {
        return my_ext->fetch(i, buffer);
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;
};

template<bool oracle_, typename Value_, typename Index_>
class AlongSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    AlongSparse(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const Index_ subset_length,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) :
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), subset_start, subset_length, opt)),
        my_shift(subset_start)
    {}

    AlongSparse(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const Index_ /* for consistency */,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) :
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), subset_start + block_start, block_length, opt)),
        my_shift(subset_start)
    {}

    AlongSparse(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const Index_ /* for consistency */,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr, 
        const Options& opt
    ) : 
        my_shift(subset_start) 
    {
        bump_indices(indices_ptr, subset_start); 
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::move(indices_ptr), opt);
    }

    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        auto output = my_ext->fetch(i, value_buffer, index_buffer);
        if (output.index && my_shift) {
            for (Index_ i = 0; i < output.number; ++i) {
                index_buffer[i] = output.index[i] - my_shift;
            }
            output.index = index_buffer;
        }
        return output;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;
    Index_ my_shift;
};

template<typename Index_>
class SubsetOracle final : public Oracle<Index_> {
public:
    SubsetOracle(
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ shift
    ) :
        my_oracle(std::move(oracle)),
        my_shift(shift)
    {}

    PredictionIndex total() const {
        return my_oracle->total();
    }

    Index_ get(const PredictionIndex i) const {
        return my_oracle->get(i) + my_shift;
    }

private:
    std::shared_ptr<const Oracle<Index_> > my_oracle;
    Index_ my_shift;
};

template<bool oracle_, typename Value_, typename Index_>
class AcrossDense final : public DenseExtractor<oracle_, Value_, Index_> {
public:
    template<typename ... Args_>
    AcrossDense(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Args_&& ... args
    ) :
        my_shift(subset_start)
    {
        if constexpr(oracle_) {
            oracle.reset(new SubsetOracle(std::move(oracle), my_shift));
        } 
        my_ext = new_extractor<false, oracle_>(matrix, row, std::move(oracle), std::forward<Args_>(args)...);
    }

    const Value_* fetch(const Index_ i, Value_* const buffer) {
        return my_ext->fetch(i + my_shift, buffer);
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;
    Index_ my_shift;
};

template<bool oracle_, typename Value_, typename Index_>
class AcrossSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    template<typename ... Args_>
    AcrossSparse(
        const Matrix<Value_, Index_>& matrix,
        const Index_ subset_start,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Args_&& ... args
    ) :
        my_shift(subset_start)
    {
        if constexpr(oracle_) {
            oracle.reset(new SubsetOracle(std::move(oracle), my_shift));
        }
        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::forward<Args_>(args)...);
    }

    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        return my_ext->fetch(i + my_shift, value_buffer, index_buffer);
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;
    Index_ my_shift;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting to a contiguous block.
 *
 * Implements delayed subsetting (i.e., slicing) of a matrix to a single contiguous block of rows or columns.
 * This is a specialized implementation that is more efficient than the `tatami::DelayedSubset` class.
 * This operation is "delayed" in that it is only evaluated when data is extracted from the matrix.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<typename Value_, typename Index_>
class DelayedSubsetBlock final : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying (pre-subset) matrix.
     * @param subset_start Index of the start of the block. This should be a row index if `by_row = true` and a column index otherwise.
     * @param subset_length Length of the block, in terms of the number of rows (if `by_row = true`) or columns (otherwise).
     * @param by_row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     */
    DelayedSubsetBlock(
        std::shared_ptr<const Matrix<Value_, Index_> > matrix,
        const Index_ subset_start,
        const Index_ subset_length,
        const bool by_row
    ) : 
        my_matrix(std::move(matrix)),
        my_subset_start(subset_start),
        my_subset_length(subset_length),
        my_by_row(by_row)
    {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > my_matrix;
    Index_ my_subset_start, my_subset_length;
    bool my_by_row;

public:
    Index_ nrow() const {
        if (my_by_row) {
            return my_subset_length;
        } else {
            return my_matrix->nrow();
        }
    }

    Index_ ncol() const {
        if (my_by_row) {
            return my_matrix->ncol();
        } else {
            return my_subset_length;
        }
    }

    bool is_sparse() const {
        return my_matrix->is_sparse();
    }

    double is_sparse_proportion() const {
        return my_matrix->is_sparse_proportion();
    }

    bool prefer_rows() const {
        return my_matrix->prefer_rows();
    }

    double prefer_rows_proportion() const {
        return my_matrix->prefer_rows_proportion();
    }

    bool uses_oracle(const bool row) const {
        return my_matrix->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::sparse_column;

    using Matrix<Value_, Index_>::sparse_row;

    /************************
     ***** Myopic dense *****
     ************************/
private:
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        const bool row,
        Args_&&... args
    ) const {
        if (row != my_by_row) {
            return std::make_unique<DelayedSubsetBlock_internal::AlongDense<oracle_, Value_, Index_> >(
                *my_matrix,
                my_subset_start,
                my_subset_length,
                row,
                std::forward<Args_>(args)...
            );
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::AcrossDense<oracle_, Value_, Index_> >(
                *my_matrix,
                my_subset_start,
                row,
                std::forward<Args_>(args)...
            );
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        const Options& opt
    ) const {
        return dense_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return dense_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return dense_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /*************************
     ***** Myopic sparse *****
     *************************/
private:
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(
        const bool row,
        Args_&&... args
    ) const {
        if (row != my_by_row) {
            return std::make_unique<DelayedSubsetBlock_internal::AlongSparse<oracle_, Value_, Index_> >(
                *my_matrix,
                my_subset_start,
                my_subset_length,
                row,
                std::forward<Args_>(args)...
            );
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::AcrossSparse<oracle_, Value_, Index_> >(
                *my_matrix,
                my_subset_start,
                row,
                std::forward<Args_>(args)...
            );
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        const Options& opt
    ) const {
        return sparse_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return sparse_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return sparse_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**************************
     ***** Oracular dense *****
     **************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        return dense_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle, 
        const Index_ block_start,
        const Index_ block_length, 
        const Options& opt
    ) const {
        return dense_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return dense_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /***************************
     ***** Oracular sparse *****
     ***************************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        return sparse_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return sparse_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return sparse_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

/**
 * @cond
 */
template<typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<const Matrix<Value_, Index_> > matrix, const Index_ subset_start, const Index_ subset_length, bool by_row) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedSubsetBlock<Value_, Index_>(std::move(matrix), subset_start, subset_length, by_row));
}

template<typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<Matrix<Value_, Index_> > matrix, Index_ subset_start, Index_ subset_length, bool by_row) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedSubsetBlock<Value_, Index_>(std::move(matrix), subset_start, subset_length, by_row));
}

template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<const Matrix<Value_, Index_> > matrix, Index_ subset_start, Index_ subset_length) {
    return make_DelayedSubsetBlock(std::move(matrix), subset_start, subset_length, margin_ == 0);
}

template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<Matrix<Value_, Index_> > matrix, Index_ subset_start, Index_ subset_length) {
    return make_DelayedSubsetBlock(std::move(matrix), subset_start, subset_length, margin_ == 0);
}
/**
 * @endcond
 */

}

#endif
