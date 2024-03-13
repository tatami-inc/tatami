#ifndef TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP

#include "../base/Matrix.hpp"
#include "utils.hpp"

#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetSortedUnique.hpp
 *
 * @brief Delayed subsetting with sorted and unique row/column indices.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedSubsetSortedUnique_internal {

template<typename Index_, class IndexStorage_>
std::vector<Index_> create(const IndexStorage_& indices) {
    return std::vector<Index_>(indices.begin(), indices.end());
}

template<typename Index_, class IndexStorage_>
std::vector<Index_> slice(const IndexStorage_& indices, Index_ block_start, Index_ block_length) {
    auto pistart = indices.begin() + block_start;
    return std::vector<Index_>(pistart, pistart + block_length);
}

template<typename Index_, class IndexStorage_>
void reindex(const IndexStorage_& indices, std::vector<Index_>& idx) {
    for (auto& i : idx) {
        i = indices[i];
    }
}

template<typename Value_, typename Index_>
void remap(SparseRange<Value_, Index_>& out, const std::vector<Index_>& remapping, Index_* ibuffer) {
    if (out.index) {
        for (Index_ i = 0; i < out.number; ++i) {
            ibuffer[i] = remapping[out.index[i]];
        }
        out.index = ibuffer;
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, const Options& opt) :
        internal(new_extractor<row_, false>(mat, create<Index_>(indices), opt)) {}

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, false>(mat, slice(indices, block_start, block_length), opt)) {}

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::vector<Index_> idx, const Options& opt) {
        reindex(indices, idx);
        internal = new_extractor<row_, false>(mat, std::move(idx), opt);
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i, buffer);
    }

protected:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct MyopicParallelSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, std::integral_constant<bool, row_>, const Options& opt) : 
        internal(new_extractor<row_, true>(mat, create<Index_>(indices), opt)), remapping(remap) {}

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) : 
        internal(new_extractor<row_, true>(mat, slice(indices, block_start, block_length), opt), remapping(remap) {}

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, std::integral_constant<bool, row_>, std::vector<Index_> idx, const Options& opt) : remapping(remap) {
        reindex(indices, idx);
        internal = new_extractor<row_, true>(mat, std::move(idx), opt);
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto out = internal->fetch(i, vbuffer, ibuffer);
        remap(out, remapping, ibuffer);
        return out;
    }

protected:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    const std::vector<Index_>& remapping;
};

template<typename Value_, typename Index_>
struct OracularParallelDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) :
        internal(new_extractor<row_, false>(mat, std::move(oracle), create<Index_>(indices), opt)) {}

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, false>(mat, std::move(oracle), slice(indices, block_start, block_length), opt)) {}

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> idx, const Options& opt) {
        reindex(indices, idx);
        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(idx), opt);
    }

    const Value_* fetch(Value_* buffer) {
        return internal->fetch(buffer);
    }

protected:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
};

template<typename Value_, typename Index_>
struct OracularParallelSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, std::integral_constant<bool, row_>, const Options& opt) : 
        internal(new_extractor<row_, true>(mat, create<Index_>(indices), opt)), remapping(remap) {}

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<row_, true>(mat, slice(indices, block_start, block_length), opt)), remapping(remap) {}

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, std::integral_constant<bool, row_>, std::vector<Index_> idx, const Options& opt) : remapping(remap) {
        reindex(indices, idx);
        internal = new_extractor<row_, true>(mat, std::move(idx), opt);
    }

    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        auto out = internal->fetch(vbuffer, ibuffer);
        remap(out, remapping, ibuffer);
        return out;
    }

protected:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
    const std::vector<Index_>& remapping;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting of a matrix with sorted, unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of sorted and unique indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetSortedUnique : public Matrix<Value_, Index_> {
    static constexpr bool storage_has_data = has_data<Index_, IndexStorage_>::value;
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * This should be sorted and unique.
     * @param check Whether to check `idx` for sorted and unique values.
     */
    DelayedSubsetSortedUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : 
        mat(std::move(p)), indices(std::move(idx)) 
    {
        if (check) {
            for (Index_ i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] <= indices[i-1]) {
                    throw std::runtime_error("indices should be unique and sorted");
                }
            }
        }

        Index_ mapping_dim = margin_ == 0 ? mat->nrow() : mat->ncol();
        mapping_single.resize(mapping_dim);
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            mapping_single[indices[i]] = i;
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;
    std::vector<Index_> mapping_single;

public:
    Index_ nrow() const {
        if constexpr(margin_==0) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }

    Index_ ncol() const {
        if constexpr(margin_==0) {
            return mat->ncol();
        } else {
            return indices.size();
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

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool accrow_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate_myopic_dense(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return populate_myopic_dense<true>(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_dense<true>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return populate_myopic_dense<true>(std::move(indices), opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return populate_myopic_dense<false>(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_dense<false>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return populate_myopic_dense<false>(std::move(indices), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate_myopic_sparse(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularSparse<Value_, Index_, IndexStorage_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return populate_myopic_sparse<true>(opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_sparse<true>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return populate_myopic_sparse<true>(std::move(indices), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return populate_myopic_sparse<false>(opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_sparse<false>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return populate_myopic_sparse<false>(std::move(indices), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
private:
    template<bool accrow_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate_oracular_dense(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::OracularPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
    }

public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_dense<true>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_dense<true>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate_oracular_dense<true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_dense<false>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_dense<false>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate_oracular_dense<false>(std::move(oracle), std::move(indices), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
private:
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate_oracular_sparse(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::OracularPerpendicularSparse<Value_, Index_, IndexStorage_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_sparse<true>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_sparse<true>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate_oracular_sparse<true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_sparse<false>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_sparse<false>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate_oracular_sparse<false>(std::move(oracle), std::move(indices), opt);
    }
};

}

#endif
