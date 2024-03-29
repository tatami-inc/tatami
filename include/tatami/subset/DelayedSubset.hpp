#ifndef TATAMI_DELAYED_SUBSET_HPP
#define TATAMI_DELAYED_SUBSET_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubset.hpp
 *
 * @brief Delayed subsetting by rows or columns.
 *
 * This is equivalent to the `DelayedSubset` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedSubset_internal {

template<typename Index_>
struct DenseParallelResults {
    std::vector<Index_> collapsed;
    std::vector<Index_> reindex;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
DenseParallelResults<Index_> format_dense_parallel_base(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(indices[to_index(i)], i);
    }
    std::sort(collected.begin(), collected.end());

    DenseParallelResults<Index_> output;
    if (collected.size()) {
        output.collapsed.reserve(len);
        output.reindex.resize(len);

        Index_ last = collected.front().first;
        output.collapsed.push_back(last);
        output.reindex[collected.front().second] = 0;

        Index_ counter = 0;
        for (Index_ i = 1; i < len; ++i) {
            const auto& pp = collected[i];
            if (pp.first != last) {
                last = pp.first;
                output.collapsed.push_back(last);
                ++counter;
            }
            output.reindex[pp.second] = counter;
        }
    }

    return output;
}

template<bool oracle_, typename Value_, typename Index_>
struct ParallelDense : public DenseExtractor<oracle_, Value_, Index_> {
    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        auto processed = format_dense_parallel_base<Index_>(indices, indices.size(), [&](Index_ i) -> Index_ { return i; });
        initialize(mat, std::move(processed), row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel_base<Index_>(indices, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
        initialize(mat, std::move(processed), row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> subset_ptr, const Options& opt) {
        const auto& subset = *subset_ptr;
        auto processed = format_dense_parallel_base<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
        initialize(mat, std::move(processed), row, std::move(oracle), opt);
    }

private:
    void initialize(const Matrix<Value_, Index_>* mat, DenseParallelResults<Index_> processed, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        vholding.resize(processed.collapsed.size());
        internal = new_extractor<false, oracle_>(mat, row, std::move(oracle), std::move(processed.collapsed), opt);
        reindex = std::move(processed.reindex);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto src = internal->fetch(i, vholding.data());

        // 'src' and 'buffer' should not point to the same array.
        auto copy = buffer;
        for (auto p : reindex) {
            *copy= src[p];
            ++copy;
        }

        return buffer;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > internal;
    std::vector<Value_> vholding;
    std::vector<Index_> reindex;
};

template<typename Index_>
struct SparseParallelReindex {
    // This is a bit complicated to explain.
    // Let 'x = pool_ptrs[i - offset]'.
    // Let 'y = pool_ptrs[i - offset + 1]'.
    // Let 'z' denote any integer in '[x, y)'.
    // In which case, 'indices[pool_indices[z]]' is equal to 'i'.
    // The general idea is that 'pool_indices[z]' can be used to fill the 'SparseRange::index' on output.
    std::vector<Index_> pool_ptrs;
    std::vector<Index_> pool_indices;
    Index_ offset;
};

template<typename Index_>
struct SparseParallelResults {
    std::vector<Index_> collapsed;
    SparseParallelReindex<Index_> reindex;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
SparseParallelResults<Index_> format_sparse_parallel_base(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        auto curdex = to_index(i);
        collected.emplace_back(indices[curdex], curdex);
    }
    std::sort(collected.begin(), collected.end());

    SparseParallelResults<Index_> output;

    if (collected.size()) {
        output.collapsed.reserve(len);
        output.reindex.pool_indices.reserve(len);
        Index_ first = collected.front().first;

        // 'pool_ptrs' is a vector that enables look-up according to the
        // indices of the underlying array. to avoid the need to allocate a
        // vector of length equal to the underlying array's dimension, we only
        // consider the extremes of 'indices'; we allocate 'pool_ptrs' to have
        // length equal to the range of 'indices' (plus 1, as we're storing
        // cumulative pointers). 'offset' defines the lower bound that must be
        // subtracted from the array indices to get an index into 'pool_ptrs'.
        output.reindex.offset = first;
        auto allocation = collected.back().first - output.reindex.offset + 1;
        output.reindex.pool_ptrs.resize(allocation + 1);

        Index_ counter = 0;
        output.reindex.pool_ptrs[counter] = 0; 
        ++counter;
        output.reindex.pool_indices.push_back(collected.front().second);
        output.reindex.pool_ptrs[counter] = 1;
        output.collapsed.push_back(first);
        auto last = first;

        for (Index_ i = 1; i < len; ++i) {
            const auto& pp = collected[i];
            auto current = pp.first;
            if (current == last) {
                output.reindex.pool_indices.push_back(pp.second);
                ++(output.reindex.pool_ptrs[counter]);
                continue;
            }

            Index_ pool_size = output.reindex.pool_indices.size();
            counter = current - output.reindex.offset;
            output.reindex.pool_ptrs[counter] = pool_size; // any overwrite is safe as the value is unchanged.
            ++counter;
            output.reindex.pool_indices.push_back(pp.second);
            output.reindex.pool_ptrs[counter] = pool_size + 1;
            output.collapsed.push_back(current);
            last = current;
        }
    }

    return output;
}

template<bool oracle_, typename Value_, typename Index_>
struct ParallelSparse : public SparseExtractor<oracle_, Value_, Index_> {
    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        auto processed = format_sparse_parallel_base<Index_>(indices, indices.size(), [](Index_ i) -> Index_ { return i; });
        initialize(mat, std::move(processed), indices.size(), row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_sparse_parallel_base<Index_>(indices, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
        initialize(mat, std::move(processed), block_length, row, std::move(oracle), opt);
    }

    template<class IndexStorage_>
    ParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> subset_ptr, const Options& opt) {
        const auto& subset = *subset_ptr;
        auto processed = format_sparse_parallel_base<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
        initialize(mat, std::move(processed), subset.size(), row, std::move(oracle), opt);
    }

private:
    void initialize(const Matrix<Value_, Index_>* mat, SparseParallelResults<Index_> processed, size_t extent, bool row, MaybeOracle<oracle_, Index_> oracle, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        needs_sort = opt.sparse_ordered_index;

        if (needs_sort && needs_value) {
            sortspace.reserve(extent);
        } 

        // We need to extract indices for sorting and expansion purposes, even
        // if they weren't actually requested.
        opt.sparse_extract_index = true;
        if (!needs_index) {
            iholding.resize(processed.collapsed.size());
        }

        internal = new_extractor<true, oracle_>(mat, row, std::move(oracle), std::move(processed.collapsed), opt);
        reindex = std::move(processed.reindex);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto vinit = (needs_value ? vbuffer + shift : NULL);
        auto iinit = (needs_index ? ibuffer + shift : iholding.data());
        auto input = internal->fetch(i, vinit, iinit);

        if (!needs_sort) {
            // Pointers in 'input' and the two 'buffer' pointers may optionally point
            // to overlapping arrays as long as each 'buffer' pointer precedes its
            // corresponding pointer in 'input'.  The idea is that the expansion of
            // values into, e.g., 'vbuffer' will cause it to catch up to 'input.value'
            // without clobbering any values in the latter. This assumes that
            // 'input.value' has been shifted enough to make space for expansion; the
            // required shift depends on the number of duplicates.
            Index_ count = 0;
            auto vcopy = vbuffer;
            auto icopy = ibuffer;

            auto vsrc = input.value;
            bool replace_value = needs_value && vsrc != vcopy;

            for (Index_ i = 0; i < input.number; ++i) {
                auto lookup = input.index[i] - reindex.offset;
                auto start = reindex.pool_ptrs[lookup];
                auto num = reindex.pool_ptrs[lookup + 1] - start;
                count += num;

                if (replace_value) {
                    auto val = *vsrc; // make a copy just in case 'vcopy' and 'input.value' overlap.
                    std::fill_n(vcopy, num, val);
                    vcopy += num;
                    ++vsrc;
                    replace_value = (vcopy != vsrc); // if we've caught up, there no need to do this replacement.
                }

                if (needs_index) {
                    // Again, 'icopy' will eventually catch up to 'input.index' if
                    // they point to overlapping arrays. But we still need to
                    // replace values once we've managed to catch up, so we can't
                    // short-circuit like we did with 'replace_value'.
                    std::copy_n(reindex.pool_indices.begin() + start, num, icopy);
                    icopy += num;
                }
            }

            input.number = count;
            if (needs_value) {
                input.value = vbuffer;
            }
            if (needs_index) {
                input.index = ibuffer;
            } else {
                input.index = NULL;
            }

        } else if (needs_value) {
            // This does not require any careful consideration of the overlaps
            // between 'input' and 'buffers', as we're copying things into
            // 'sortspace' anyway before copying them back into 'buffer'.
            sortspace.clear();
            for (Index_ i = 0; i < input.number; ++i) {
                auto val = input.value[i];
                auto lookup = input.index[i] - reindex.offset;
                auto start = reindex.pool_ptrs[lookup];
                auto end = reindex.pool_ptrs[lookup + 1];
                for (Index_ j = start; j < end; ++j) {
                    sortspace.emplace_back(reindex.pool_indices[j], val);
                }
            }
            std::sort(sortspace.begin(), sortspace.end());
            input.number = sortspace.size();

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

        } else {
            // Again, 'input.index' and 'ibuffer' may point to overlapping arrays,
            // as long as the latter precedes the former; expansion into the latter
            // will allow it to catch up to the former without clobbering, assuming 
            // that the latter was shifted back to provide enough space. 
            Index_ count = 0;
            auto icopy = ibuffer;

            for (Index_ i = 0; i < input.number; ++i) {
                auto lookup = input.index[i] - reindex.offset;
                auto start = reindex.pool_ptrs[lookup];
                auto num = reindex.pool_ptrs[lookup + 1] - start;
                count += num;

                if (needs_index) {
                    std::copy_n(reindex.pool_indices.begin() + start, num, icopy);
                    icopy += num;
                }
            }

            input.number = count;
            if (needs_index) {
                std::sort(ibuffer, ibuffer + count);
                input.index = ibuffer;
            } else {
                input.index = NULL;
            }
        }

        return input;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > internal;
    bool needs_value, needs_index, needs_sort;
    SparseParallelReindex<Index_> reindex;
    std::vector<std::pair<Index_, Value_> > sortspace;
    std::vector<Index_> iholding;
    size_t shift;
};

}
/**
 * @endcond
 */

/**
 * @brief delayed subsetting of a matrix with general indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of arbitrary indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubset : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * These may be duplicated and/or unsorted.
     */
    DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx) : mat(std::move(p)), indices(std::move(idx)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

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
    template<typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > populate_myopic_dense(bool row, Args_&& ... args) const {
        if (row == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, row, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelDense<false, Value_, Index_> >(mat.get(), indices, row, false, std::forward<Args_>(args)...);
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
        if (row == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularSparse<Value_, Index_, IndexStorage_> >(mat.get(), indices, row, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelSparse<false, Value_, Index_> >(mat.get(), indices, row, false, std::forward<Args_>(args)...);
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
        if (row == (margin_ == 0)) {
            return std::make_unique<subset_utils::OracularPerpendicularDense<Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelDense<true, Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...);
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
        if (row == (margin_ == 0)) {
            return std::make_unique<subset_utils::OracularPerpendicularSparse<Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelSparse<true, Value_, Index_> >(mat.get(), indices, row, std::move(oracle), std::forward<Args_>(args)...);
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
