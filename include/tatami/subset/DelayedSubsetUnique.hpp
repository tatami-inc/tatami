#ifndef TATAMI_DELAYED_SUBSET_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_UNIQUE_HPP

#include "utils.hpp"
#include "../base/Matrix.hpp"
#include "../utils/copy.hpp"

#include <algorithm>
#include <numeric>
#include <memory>

/**
 * @file DelayedSubsetUnique.hpp
 *
 * @brief Delayed subsetting by unique row/column indices.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedSubsetUnique_internal {

template<typename Index_>
struct DenseParallelResults {
    std::vector<Index_> sorted;
    std::vector<Index_> permutation;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
DenseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(indices[to_index(i)], i);
    }
    std::sort(collected.begin(), collected.end());

    DenseParallelResults<Index_> output;
    output.sorted.reserve(len);
    output.permutation.reserve(len);
    for (const auto& pp : collected) {
        output.sorted.push_back(pp.first);
        output.permutation.push_back(pp.second);
    }

    return output;
}

template<bool oracle_, typename Value_, typename Index_>
struct ParallelDense : DenseExtractor<oracle_, Value_, Index_> {
    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, indices.size(), [&](Index_ i) -> Index_ { return i; });
        initialize(mat, std::move(processed), indices.size(), row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
        initialize(mat, std::move(processed), block_length, row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> subset_ptr, const Options& opt) {
        const auto& subset = *subset_ptr;
        auto processed = format_dense_parallel<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
        initialize(mat, std::move(processed), subset.size(), row, std::move(oracle), opt);
    }

private:
    void initialize(const Matrix<Value_, Index_>* mat, DenseParallelResults<Index_> processed, size_t extent, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        internal = new_extractor<false, oracle_>(mat, row, std::move(oracle), std::move(processed.sorted), opt);
        vholding.resize(extent);
        permutation = std::move(processed.permutation);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto src = internal->fetch(i, vholding.data());

        // 'input' and 'output' should not point to the same array. In theory, it
        // is possible to do an in-place permutation, but this requires another
        // array anyway to track the permutation status, so we'll just keep it simple.
        for (auto p : permutation) {
            buffer[p] = *src;
            ++src;
        }

        return buffer;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;
    std::vector<Value_> vholding;
    std::vector<Index_> permutation;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
std::vector<Index_> format_sparse_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
    std::vector<Index_> collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(indices[to_index(i)]);
    }
    std::sort(collected.begin(), collected.end());
    return collected;
}

template<bool oracle_, typename Value_, typename Index_>
struct ParallelSparse : public SparseExtractor<oracle_, Value_, Index_> {
    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) : remapping(remap) {
        auto processed = format_sparse_parallel<Index_>(indices, indices.size(), [&](Index_ i) -> Index_ { return i; });
        initialize(mat, std::move(processed), indices.size(), row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) : remapping(remap) {
        auto processed = format_sparse_parallel<Index_>(indices, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
        initialize(mat, std::move(processed), block_length, row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, const std::vector<Index_>& remap, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> subset_ptr, const Options& opt) : remapping(remap) {
        const auto& subset = *subset_ptr;
        auto processed = format_sparse_parallel<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
        initialize(mat, std::move(processed), subset.size(), row, std::move(oracle), opt);
    }

private:
    void initialize(const Matrix<Value_, Index_>* mat, std::vector<Index_> sorted, size_t extent, bool row, MaybeOracle<oracle_, Index_> oracle, Options opt) {
        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        needs_sort = opt.sparse_ordered_index;
 
        // The conditionals here mirror those in 'fetch',
        // to self-document the case where each of the temporaries are needed.
        if (!needs_sort) {
            if (needs_index) {
                ; // no 'iholding' required as a user-provided 'ibuffer' should be available.
            }

        } else if (needs_value) {
            opt.sparse_extract_index = true;
            sortspace.reserve(extent);
            if (needs_index) {
                ; // no 'iholding' required as a user-provided 'ibuffer' should be available.
            } else {
                iholding.resize(extent); // needs 'iholding' as user-provided 'ibuffer' may be NULL.
            }

        } else if (needs_index) {
            ; // no 'iholding' required as a user-provided 'ibuffer' should be available.
        }

        internal = new_extractor<true, oracle_>(mat, row, std::move(oracle), std::move(sorted), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto input = internal->fetch(i, vbuffer, (iholding.empty() ? ibuffer : iholding.data()));

        // Pointers in 'input' and the 'buffer' pointers may point to the same array,
        // as we're either just modifiying in place or we're copying to 'sortspace'.
        if (!needs_sort) {
            if (needs_index) {
                for (Index_ i = 0; i < input.number; ++i) {
                    ibuffer[i] = remapping[input.index[i]];
                }
                input.index = ibuffer;
            }

        } else if (needs_value) {
            // We assume that the indices have already been extracted for sorting
            // purposes, even if they weren't actually requested.
            sortspace.clear();
            for (Index_ i = 0; i < input.number; ++i) {
                sortspace.emplace_back(remapping[input.index[i]], input.value[i]);
            }
            std::sort(sortspace.begin(), sortspace.end());

            auto vcopy = vbuffer;
            for (const auto& ss : sortspace) {
                *vcopy = ss.second;
                ++vcopy;
            }
            input.value = vbuffer;

            if (needs_index) {
                auto icopy = ibuffer;
                for (const auto& ss : sortspace) {
                    *icopy = ss.first;
                    ++icopy;
                }
                input.index = ibuffer;
            } else {
                input.index = NULL;
            }

        } else if (needs_index) {
            for (Index_ i = 0; i < input.number; ++i) {
                ibuffer[i] = remapping[input.index[i]];
            }
            std::sort(ibuffer, ibuffer + input.number);
            input.index = ibuffer;
        }

        return input;
    }

private:
    const std::vector<Index_>& remapping;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;
    bool needs_value, needs_index, needs_sort;
    std::vector<std::pair<Index_, Value_> > sortspace;
    std::vector<Index_> iholding;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting of a matrix with unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of unique indices.
 * This operation is "delayed" in that it is only evaluated when rows or columns are requested from the matrix.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 */
template<typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetUnique : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `row = true`) or columns (otherwise).
     * This should be unique, but may be unsorted.
     * @param row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     * @param check Whether to check `idx` for unique values.
     */
    DelayedSubsetUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool row, bool check = true) : 
        mat(std::move(p)), indices(std::move(idx)), by_row(row)
    {
        Index_ fulldim = by_row ? mat->nrow() : mat->ncol();

        if (check) {
            std::vector<unsigned char> checks(fulldim);
            for (Index_ i = 0, end = indices.size(); i < end; ++i) {
                auto& found = checks[indices[i]];
                if (found) {
                    throw std::runtime_error("indices should be unique");
                } 
                found = 1;
            }
        }

        mapping_single.resize(fulldim);
        for (Index_  i = 0, end = indices.size(); i < end; ++i) {
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
    template<typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > populate_myopic_dense(bool row, Args_&& ... args) const {
        if (row == by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, row, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetUnique_internal::ParallelDense<false, Value_, Index_> >(mat.get(), indices, row, false, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelSparse<false, Value_, Index_> >(mat.get(), indices, mapping_single, row, false, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelDense<true, Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelSparse<true, Value_, Index_> >(mat.get(), indices, mapping_single, row, std::move(oracle), std::forward<Args_>(args)...);
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
