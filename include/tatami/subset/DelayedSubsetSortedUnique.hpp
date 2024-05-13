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
VectorPtr<Index_> create(const IndexStorage_& indices) {
    return std::make_shared<std::vector<Index_> >(indices.begin(), indices.end());
}

template<typename Index_, class IndexStorage_>
VectorPtr<Index_> create(const IndexStorage_& indices, Index_ block_start, Index_ block_length) {
    auto pistart = indices.begin() + block_start;
    return std::make_shared<std::vector<Index_> >(pistart, pistart + block_length);
}

template<typename Index_, class IndexStorage_>
VectorPtr<Index_> create(const IndexStorage_& indices, const VectorPtr<Index_>& idx_ptr) {
    auto rawptr = std::make_shared<std::vector<Index_> >();
    VectorPtr<Index_> outptr(rawptr);
    auto& output = *rawptr;

    const auto& input = *idx_ptr;
    output.reserve(input.size());
    for (auto i : input) {
        output.push_back(indices[i]);
    }

    return outptr;
}

template<bool oracle_, typename Value_, typename Index_>
struct ParallelDense : public DenseExtractor<oracle_, Value_, Index_> {
    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) :
        internal(new_extractor<false, oracle_>(mat, row, std::move(oracle), create<Index_>(indices), opt)) {}

    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) :
        internal(new_extractor<false, oracle_>(mat, row, std::move(oracle), create<Index_>(indices, block_start, block_length), opt)) {}

    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> idx_ptr, const Options& opt) :
        internal(new_extractor<false, oracle_>(mat, row, std::move(oracle), create<Index_>(indices, idx_ptr), opt)) {}

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        return internal->fetch(i, buffer);
    }

protected:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;
};

template<bool oracle_, typename Value_, typename Index_>
struct ParallelSparse : public SparseExtractor<oracle_, Value_, Index_> {
    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) : 
        internal(new_extractor<true, oracle_>(mat, row, std::move(oracle), create<Index_>(indices), opt)), remapping(remap) {}

    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) : 
        internal(new_extractor<true, oracle_>(mat, row, std::move(oracle), create<Index_>(indices, block_start, block_length), opt)), remapping(remap) {}

    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> idx_ptr, const Options& opt) : 
        internal(new_extractor<true, oracle_>(mat, row, std::move(oracle), create<Index_>(indices, idx_ptr), opt)), remapping(remap) {}

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto out = internal->fetch(i, vbuffer, ibuffer);
        if (out.index) {
            for (Index_ i = 0; i < out.number; ++i) {
                ibuffer[i] = remapping[out.index[i]];
            }
            out.index = ibuffer;
        }
        return out;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;
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
 * This operation is "delayed" in that it is only evaluated when data is requested from the matrix.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 */
template<typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetSortedUnique : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `row = true`) or columns (otherwise).
     * This should be sorted and unique.
     * @param row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     * @param check Whether to check `idx` for sorted and unique values.
     */
    DelayedSubsetSortedUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool row, bool check = true) :
        mat(std::move(p)), indices(std::move(idx)), by_row(row)
    {
        if (check) {
            for (Index_ i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] <= indices[i-1]) {
                    throw std::runtime_error("indices should be unique and sorted");
                }
            }
        }

        Index_ mapping_dim = by_row ? mat->nrow() : mat->ncol();
        mapping_single.resize(mapping_dim);
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            mapping_single[indices[i]] = i;
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;
    bool by_row;
    std::vector<Index_> mapping_single;

public:
    Index_ nrow() const {
        if (by_row) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }

    Index_ ncol() const {
        if (by_row) {
            return mat->ncol();
        } else {
            return indices.size();
        }
    }

    bool is_sparse() const {
        return mat->is_sparse();
    }

    double is_sparse_proportion() const {
        return mat->is_sparse_proportion();
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
    template<typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > populate_myopic_dense(bool row, Args_&& ... args) const {
        if (row == by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, row, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::ParallelDense<false, Value_, Index_> >(mat.get(), indices, row, false, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options& opt) const {
        return populate_myopic_dense(row, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_dense(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return populate_myopic_dense(row, std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    template<typename ... Args_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > populate_myopic_sparse(bool row, Args_&& ... args) const {
        if (row == by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularSparse<Value_, Index_, IndexStorage_> >(mat.get(), indices, row, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::ParallelSparse<false, Value_, Index_> >(mat.get(), indices, mapping_single, row, false, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        return populate_myopic_sparse(row, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_sparse(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return populate_myopic_sparse(row, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
private:
    template<typename ... Args_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > populate_oracular_dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Args_&& ... args) const {
        if (row == by_row) {
            return std::make_unique<subset_utils::OracularPerpendicularDense<Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::ParallelDense<true, Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_dense(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_dense(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return populate_oracular_dense(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
private:
    template<typename ... Args_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > populate_oracular_sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Args_&& ... args) const {
        if (row == by_row) {
            return std::make_unique<subset_utils::OracularPerpendicularSparse<Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSortedUnique_internal::ParallelSparse<true, Value_, Index_> >(mat.get(), indices, mapping_single, row, std::move(oracle), std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_sparse(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_sparse(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return populate_oracular_sparse(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

}

#endif
