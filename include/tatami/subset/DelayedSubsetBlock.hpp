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
void bump_indices(std::vector<Index_>& indices, Index_ subset_start) {
    if (subset_start) {
        for (auto& i : indices) {
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

template<typename Value_, typename Index_>
struct MyopicAlongDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicAlongDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, Index_ subset_length, const Options& opt) :
        internal(new_extractor<row_, false>(mat, subset_start, subset_length, opt)) {}

    template<bool row_>
    MyopicAlongDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, false>(mat, subset_start + block_start, block_length, opt)) {}

    template<bool row_>
    MyopicAlongDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, std::vector<Index_> indices, const Options& opt) {
        bump_indices(indices, subset_start); 
        internal = new_extractor<row_, false>(mat, std::move(indices), opt);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i, buffer);
    }

private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct OracularAlongDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_>
    OracularAlongDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, Index_ subset_length, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) :
        internal(new_extractor<row_, false>(mat, std::move(oracle), subset_start, subset_length, opt)) {}

    template<bool row_>
    OracularAlongDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, false>(mat, std::move(oracle), subset_start + block_start, block_length, opt)) {}

    template<bool row_>
    OracularAlongDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) {
        bump_indices(indices, subset_start); 
        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(indices), opt);
    }

public:
    const Value_* fetch(Value_* buffer) {
        return internal->fetch(buffer);
    }

private:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct MyopicAlongSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicAlongSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, Index_ subset_length, const Options& opt) :
        internal(new_extractor<row_, true>(mat, subset_start, subset_length, opt)), shift(subset_start) {}

    template<bool row_>
    MyopicAlongSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, true>(mat, subset_start + block_start, block_length, opt)), shift(subset_start) {}

    template<bool row_>
    MyopicAlongSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, std::vector<Index_> indices, const Options& opt) : 
        shift(subset_start) 
    {
        bump_indices(indices, subset_start); 
        internal = new_extractor<row_, true>(mat, std::move(indices), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto output = internal->fetch(i, vbuffer, ibuffer);
        debump_indices(output, ibuffer, shift);
        return output;
    }

private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    Index_ shift;
};

template<typename Value_, typename Index_>
struct OracularAlongSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_>
    OracularAlongSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, Index_ subset_length, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) :
        internal(new_extractor<row_, true>(mat, std::move(oracle), subset_start, subset_length, opt)), shift(subset_start) {}

    template<bool row_>
    OracularAlongSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, true>(mat, std::move(oracle), subset_start + block_start, block_length, opt)), shift(subset_start) {}

    template<bool row_>
    OracularAlongSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, [[maybe_unused]] Index_ subset_length, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) : 
        shift(subset_start) 
    {
        bump_indices(indices, subset_start); 
        internal = new_extractor<row_, true>(mat, std::move(oracle), std::move(indices), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        auto output = internal->fetch(vbuffer, ibuffer);
        debump_indices(output, ibuffer, shift);
        return output;
    }

private:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
    Index_ shift;
};

template<typename Value_, typename Index_>
struct MyopicAcrossDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    MyopicAcrossDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, Args_&& ... args) :
        internal(new_extractor<row_, false>(mat, std::forward<Args_>(args)...)), shift(subset_start) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i + shift, buffer);
    }

private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
    Index_ shift;
};

template<typename Value_, typename Index_>
struct MyopicAcrossSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    MyopicAcrossSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, Args_&& ... args) :
        internal(new_extractor<row_, true>(mat, std::forward<Args_>(args)...)), shift(subset_start) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(i + shift, vbuffer, ibuffer);
    }

private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    Index_ shift;
};

template<typename Index_>
struct SubsetOracle : public Oracle<Index_> {
    SubsetOracle(std::shared_ptr<Oracle<Index_> > input, Index_ shift) : input(std::move(input)), shift(shift) {}

    size_t total() const {
        return input->total();
    }

    Index_ get(size_t i) const {
        return input->get(i) + shift;
    }

private:
    std::shared_ptr<Oracle<Index_> > input;
    Index_ shift;
};

template<typename Value_, typename Index_>
struct OracularAcrossDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    OracularAcrossDense(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, std::shared_ptr<Oracle<Index_> > oracle, Args_&& ... args) :
        internal(new_extractor<row_, false>(mat, std::make_shared<SubsetOracle<Index_> > (std::move(oracle), subset_start), std::forward<Args_>(args)...)) {}

    const Value_* fetch(Value_* buffer) {
        return internal->fetch(buffer);
    }

private:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct OracularAcrossSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    OracularAcrossSparse(const Matrix<Value_, Index_>* mat, std::integral_constant<bool, row_>, Index_ subset_start, std::shared_ptr<Oracle<Index_> > oracle, Args_&& ... args) :
        internal(new_extractor<row_, true>(mat, std::make_shared<SubsetOracle<Index_> >(std::move(oracle), subset_start), std::forward<Args_>(args)...)) {}

    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        return internal->fetch(vbuffer, ibuffer);
    }

private:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
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
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > myopic_dense_internal(Args_&&... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ != (margin_ == 0)) {
            return std::make_unique<DelayedSubsetBlock_internal::MyopicAlongDense<Value_, Index_> >(mat.get(), flag, block_start, block_length, std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::MyopicAcrossDense<Value_, Index_> >(mat.get(), flag, block_start, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return myopic_dense_internal<true>(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return myopic_dense_internal<true>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return myopic_dense_internal<true>(std::move(indices), opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return myopic_dense_internal<false>(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return myopic_dense_internal<false>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return myopic_dense_internal<false>(std::move(indices), opt);
    }

    /*************************
     ***** Myopic sparse *****
     *************************/
private:
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > myopic_sparse_internal(Args_&&... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ != (margin_ == 0)) {
            return std::make_unique<DelayedSubsetBlock_internal::MyopicAlongSparse<Value_, Index_> >(mat.get(), flag, block_start, block_length, std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::MyopicAcrossSparse<Value_, Index_> >(mat.get(), flag, block_start, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return myopic_sparse_internal<true>(opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return myopic_sparse_internal<true>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return myopic_sparse_internal<true>(std::move(indices), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return myopic_sparse_internal<false>(opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return myopic_sparse_internal<false>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return myopic_sparse_internal<false>(std::move(indices), opt);
    }

    /**************************
     ***** Oracular dense *****
     **************************/
private:
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > oracular_dense_internal(Args_&&... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ != (margin_ == 0)) {
            return std::make_unique<DelayedSubsetBlock_internal::OracularAlongDense<Value_, Index_> >(mat.get(), flag, block_start, block_length, std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::OracularAcrossDense<Value_, Index_> >(mat.get(), flag, block_start, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return oracular_dense_internal<true>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return oracular_dense_internal<true>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return oracular_dense_internal<true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return oracular_dense_internal<false>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return oracular_dense_internal<false>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return oracular_dense_internal<false>(std::move(oracle), std::move(indices), opt);
    }

    /***************************
     ***** Oracular sparse *****
     ***************************/
private:
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > oracular_sparse_internal(Args_&&... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ != (margin_ == 0)) {
            return std::make_unique<DelayedSubsetBlock_internal::OracularAlongSparse<Value_, Index_> >(mat.get(), flag, block_start, block_length, std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetBlock_internal::OracularAcrossSparse<Value_, Index_> >(mat.get(), flag, block_start, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return oracular_sparse_internal<true>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return oracular_sparse_internal<true>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return oracular_sparse_internal<true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return oracular_sparse_internal<false>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return oracular_sparse_internal<false>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return oracular_sparse_internal<false>(std::move(oracle), std::move(indices), opt);
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
