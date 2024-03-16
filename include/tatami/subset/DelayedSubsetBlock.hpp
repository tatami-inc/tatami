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
void bump_indices(VectorPtr<Index_>& indices_ptr, Index_ subset_start) {
    if (subset_start) {
        auto ptr2 = new std::vector<Index_>(*indices_ptr);
        indices_ptr.reset(ptr2);
        for (auto& i : *ptr2) {
            i += subset_start;
        }
    }
}

template<typename Value_, typename Index_>
void debump_indices(SparseRange<Value_, Index_>& range, Index_* buffer, Index_ subset_start) {
    if (range.index && subset_start) {
        for (Index_ i = 0; i < range.number; ++i) {
            buffer[i] = range.index[i] - subset_start;
        }
        range.index = buffer;
    }
}

template<bool oracle_, typename Value_, typename Index_>
struct AlongDense : public DenseExtractor<oracle_, Value_, Index_> {
    AlongDense(const Matrix<Value_, Index_>* mat, Index_ subset_start, Index_ subset_length, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) :
        internal(new_extractor<false, oracle_>(mat, row, std::move(oracle), subset_start, subset_length, opt)) {}

    AlongDense(const Matrix<Value_, Index_>* mat, Index_ subset_start, [[maybe_unused]] Index_ subset_length, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<false, oracle_>(mat, row, std::move(oracle), subset_start + block_start, block_length, opt)) {}

    AlongDense(const Matrix<Value_, Index_>* mat, Index_ subset_start, [[maybe_unused]] Index_ subset_length, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) {
        bump_indices(indices_ptr, subset_start); 
        internal = new_extractor<false, oracle_>(mat, row, std::move(oracle), std::move(indices_ptr), opt);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i, buffer);
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;
};

template<bool oracle_, typename Value_, typename Index_>
struct AlongSparse : public SparseExtractor<oracle_, Value_, Index_> {
    AlongSparse(const Matrix<Value_, Index_>* mat, Index_ subset_start, Index_ subset_length, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) :
        internal(new_extractor<true, oracle_>(mat, row, std::move(oracle), subset_start, subset_length, opt)), shift(subset_start) {}

    AlongSparse(const Matrix<Value_, Index_>* mat, Index_ subset_start, [[maybe_unused]] Index_ subset_length, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<true, oracle_>(mat, row, std::move(oracle), subset_start + block_start, block_length, opt)), shift(subset_start) {}

    AlongSparse(const Matrix<Value_, Index_>* mat, Index_ subset_start, [[maybe_unused]] Index_ subset_length, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) : 
        shift(subset_start) 
    {
        bump_indices(indices_ptr, subset_start); 
        internal = new_extractor<true, oracle_>(mat, row, std::move(oracle), std::move(indices_ptr), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto output = internal->fetch(i, vbuffer, ibuffer);
        debump_indices(output, ibuffer, shift);
        return output;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;
    Index_ shift;
};

template<typename Index_>
struct SubsetOracle : public Oracle<Index_> {
    SubsetOracle(std::shared_ptr<const Oracle<Index_> > input, Index_ shift) : input(std::move(input)), shift(shift) {}

    size_t total() const {
        return input->total();
    }

    Index_ get(size_t i) const {
        return input->get(i) + shift;
    }

private:
    std::shared_ptr<const Oracle<Index_> > input;
    Index_ shift;
};

template<bool oracle_, typename Value_, typename Index_>
struct AcrossDense : public DenseExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    AcrossDense(const Matrix<Value_, Index_>* mat, Index_ subset_start, bool row, MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) : shift(subset_start) {
        if constexpr(oracle_) {
            auto ptr = new SubsetOracle(std::move(oracle), shift);
            oracle.reset(ptr);
        } 
        internal = new_extractor<false, oracle_>(mat, row, std::move(oracle), std::forward<Args_>(args)...);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i + shift, buffer);
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;
    Index_ shift;
};

template<bool oracle_, typename Value_, typename Index_>
struct AcrossSparse : public SparseExtractor<oracle_, Value_, Index_> {
    template<typename ... Args_>
    AcrossSparse(const Matrix<Value_, Index_>* mat, Index_ subset_start, bool row, MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) : shift(subset_start) {
        if constexpr(oracle_) {
            auto ptr = new SubsetOracle(std::move(oracle), shift);
            oracle.reset(ptr);
        }
        internal = new_extractor<true, oracle_>(mat, row, std::move(oracle), std::forward<Args_>(args)...);
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(i + shift, vbuffer, ibuffer);
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;
    Index_ shift;
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
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<int margin_, typename Value_, typename Index_>
class DelayedSubsetBlock : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param s Index of the start of the block. This should be a row index if `margin_ = 0` and a column index otherwise.
     * @param l Length of the block, in terms of the number of rows (if `margin_ = 0`) or columns (otherwise).
     */
    DelayedSubsetBlock(std::shared_ptr<const Matrix<Value_, Index_> > p, Index_ s, Index_ l) : mat(std::move(p)), block_start(s), block_length(l) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    Index_ block_start, block_length;

public:
    Index_ nrow() const {
        if constexpr(margin_==0) {
            return block_length;
        } else {
            return mat->nrow();
        }
    }

    Index_ ncol() const {
        if constexpr(margin_==0) {
            return mat->ncol();
        } else {
            return block_length;
        }
    }

    bool sparse() const {
        return mat->sparse();
    }

    double sparse_proportion() const {
        return mat->sparse_proportion();
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

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::sparse_column;

    using Matrix<Value_, Index_>::sparse_row;

    /************************
     ***** Myopic dense *****
     ************************/
private:
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, Args_&&... args) const {
        if (row != (margin_ == 0)) {
            return std::make_unique<DelayedSubsetBlock_internal::AlongDense<oracle_, Value_, Index_> >(mat.get(), block_start, block_length, row, std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::AcrossDense<oracle_, Value_, Index_> >(mat.get(), block_start, row, std::forward<Args_>(args)...);
        }
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

    /*************************
     ***** Myopic sparse *****
     *************************/
private:
    template<bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, Args_&&... args) const {
        if (row != (margin_ == 0)) {
            return std::make_unique<DelayedSubsetBlock_internal::AlongSparse<oracle_, Value_, Index_> >(mat.get(), block_start, block_length, row, std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::AcrossSparse<oracle_, Value_, Index_> >(mat.get(), block_start, row, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        return sparse_internal<false>(row, false, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices, const Options& opt) const {
        return sparse_internal<false>(row, false, std::move(indices), opt);
    }

    /**************************
     ***** Oracular dense *****
     **************************/
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

    /***************************
     ***** Oracular sparse *****
     ***************************/
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
 * @tparam margin_ Dimension along which the addition is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 *
 * @param p Pointer to the underlying (pre-subset) `Matrix`.
 * @param f Index of the start of the block. This should be a row index if `margin_ = 0` and a column index otherwise.
 * @param l Index of the one-past-the-end of the block.
 *
 * @return A pointer to a `DelayedSubsetBlock` instance.
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<const Matrix<Value_, Index_> > p, Index_ f, Index_ l) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedSubsetBlock<margin_, Value_, Index_>(std::move(p), f, l));
}

/**
 * @cond
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<Matrix<Value_, Index_> > p, Index_ f, Index_ l) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedSubsetBlock<margin_, Value_, Index_>(std::move(p), f, l));
}
/**
 * @endcond
 */

}

#endif
