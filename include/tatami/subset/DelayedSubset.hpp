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
DenseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
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

template<typename Index_, class IndexStorage_>
DenseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices) {
    return format_dense_parallel<Index_>(indices, indices.size(), [&](Index_ i) -> Index_ { return i; });
}

template<typename Index_, class IndexStorage_>
DenseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices, Index_ start, Index_ length) {
    return format_dense_parallel<Index_>(indices, length, [&](Index_ i) -> Index_ { return i + start; });
}

template<typename Index_, class IndexStorage_>
DenseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices, const std::vector<Index_>& subset) {
    return format_dense_parallel<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
}

template<typename Value_, typename Index_>
void reorder_dense_parallel(const Value_* input, Value_* output, const std::vector<Index_>& reindex) {
    // 'input' and 'output' should not point to the same array.
    for (auto p : reindex) {
        *output = input[p];
        ++output;
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelDense : MyopicDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices);
        initialize<row_>(processed, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, subset);
        initialize<row_>(processed, opt);
    }

private:
    template<bool row_>
    void initialize(DenseParallelResults<Index_>& processed, const Options& opt) {
        vholding.resize(processed.collapsed.size());
        internal = new_extractor<row_, false>(mat, std::move(processed.collapsed), opt);
        reindex = std::move(processed.reindex);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto src = internal->fetch(i, vholding.data());
        expand_dense_parallel(src, buffer, reindex);
        return buffer;
    }

private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
    std::vector<Value_> vholding;
    std::vector<Index_> reindex;
};

/**********************
 *** Oracular dense ***
 **********************/

template<typename Value_, typename Index_>
struct OracularParallelDense : OracularDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, 0, indices.size());
        initialize<row_>(processed, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, subset);
        initialize<row_>(processed, std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(DenseParallelResults<Index_>& processed, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) { 
        vholding.resize(processed.collapsed.size());
        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(processed.collapsed), opt);
        reindex = std::move(processed.reindex);
    }

public:
    const Value_* fetch(Value_* buffer) {
        auto src = internal->fetch(vholding.data());
        expand_dense_parallel(src, buffer, reindex);
        return buffer;
    }

private:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
    std::vector<Value_> vholding;
    std::vector<Index_> reindex;
};

/*********************
 *** Myopic sparse ***
 *********************/

template<typename Index_>
struct SparseParallelReindex {
    // This is a bit complicated to explain.
    // Let 'x = pool_ptrs[i - offset]'.
    // Let 'y = pool_ptrs[i - offset + 1]'.
    // Let 'z' denote any integer in '[x, y)'.
    // In which case, 'indices[pool_index[z]]' is equal to 'i'.
    std::vector<Index_> pool_ptrs;
    std::vector<Index_> pool_index;
    Index_ offset;
};

template<typename Index_>
struct SparseParallelResults {
    std::vector<Index_> collapsed;
    SparseParallelReindex reindex;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
SparseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(indices[to_index(i)], i);
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

        for (index_ i = 1 i < len; ++i) {
            const auto& pp = collected[i];
            auto current = pp.first;
            if (current == last) {
                output.reindex.pool_indices.push_back(pp.second);
                ++(output.reindex.pool_ptrs[counter])
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

template<typename Index_, class IndexStorage_>
SparseParallelResults<Index_> format_sparse_parallel(const IndexStorage_& indices) {
    return format_sparse_parallel<Index_>(indices, indices.size(), [&](Index_ i) -> Index_ { return i; });
}

template<typename Index_, class IndexStorage_>
SparseParallelResults<Index_> format_sparse_parallel(const IndexStorage_& indices, Index_ start, Index_ length) {
    return format_sparse_parallel<Index_>(indices, length, [&](Index_ i) -> Index_ { return i + start; });
}

template<typename Index_, class IndexStorage_>
SparseParallelResults<Index_> format_sparse_parallel(const IndexStorage_& indices, const std::vector<Index_>& subset) {
    return format_sparse_parallel<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
}

template<typename Value_, typename Index_>
SparseRange<Value_, Index_> reorder_sparse_parallel(
    SparseRange<Value_, Index_>& input,
    Value_* vbuffer, 
    Index_* ibuffer, 
    bool needs_value,
    bool needs_index,
    const SparseParallelReindex<Index_>& reindex)
{
    // Pointers in 'input' and the two 'buffer' pointers may optionally point
    // to overlapping arrays as long as each 'buffer' pointer precede its
    // corresponding pointer in 'input'.  The idea is that the expansion of
    // values into 'buffer' will cause it to "catch up" to 'input' without
    // clobbering any values in the latter. This assumes that 'input' has been
    // shifted enough to make space for expansion; the required shift depends
    // on the number of duplicates.
    if (!needs_sort) {
        Index_ count = 0;
        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        bool replace_value = needs_value;

        for (Index_ i = 0; i < input.number; ++i) {
            auto lookup = input.index[i] - reindex.offset;
            auto start = reindex.pool_ptrs[lookup];
            auto num = reindex.pool_ptrs[lookup + 1] - start;
            count += num;

            if (replace_value) {
                auto ivptr = input.value + i;
                auto val = *ivptr; // copy it out just in case 'vcopy' and 'input.value' overlap.
                std::fill_n(vcopy, vcopy + num, val);
                vcopy += num;
                replace_value = (vcopy != ivptr); // if we've caught up, there no need to do this replacement.
            }

            if (needs_index) {
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
        // We assume that the indices have already been extracted for sorting
        // purposes, even if they weren't actually requested.
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
        // Again, we assume that the indices have already been extracted for
        // sorting purposes, even if they weren't actually requested.
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
}

template<typename Value_, typename Index_>
struct MyopicParallelSparse : MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices);
        initialize<row_>(processed, indices.size(), opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(
        const Matrix<Value_, Index_>* mat, 
        const IndexStorage_& indices, 
        const std::vector<Index_>& remap, 
        std::integral_constant<bool, row_>, 
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        remapping(remap) 
    {
        auto processed = format_sparse_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, block_length, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(
        const Matrix<Value_, Index_>* mat, 
        const IndexStorage_& indices, 
        const std::vector<Index_>& remap, 
        std::integral_constant<bool, row_>, 
        std::vector<Index_> subset, 
        const Options& opt) :
        remapping(remap)
    {
        auto processed = format_sparse_parallel<Index_>(indices, subset);
        initialize<row_>(processed, subset.size(), opt);
    }

private:
    template<bool row_>
    void initialize(SparseParallelResults<Index_>& processed, size_t extent, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        needs_sort = opt.sparse_ordered_index;

        if (needs_sort && needs_value) {
            sortspace.reserve(extent);
        } 
        if (!needs_index) {
            iholding.resize(processed.collapsed.size());
        }

        internal = new_extractor<row_, false>(mat, std::move(processed.collapsed), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto vinit = (needs_value ? vbuffer + shift : NULL);
        auto iinit = (needs_index ? ibuffer + shift : iholding.data());
        auto src = internal->fetch(i, vinit, iinit);
        reorder_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, needs_sort, sortspace); 
        return src;
    }

private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    bool needs_value, needs_index, needs_sort;
    std::vector<std::pair<Index_, Value_> > sortspace;
    std::vector<Index_> iholding;
    size_t shift;
};

/***********************
 *** Oracular sparse ***
 ***********************/

template<typename Value_, typename Index_>
struct OracularParallelSparse : OracularSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelSparse(
        const Matrix<Value_, Index_>* mat, 
        const IndexStorage_& indices, 
        std::integral_constant<bool, row_>, 
        std::shared_ptr<Oracle<Index_> > oracle, 
        const Options& opt) 
    {
        auto processed = format_sparse_parallel<Index_>(indices);
        initialize<row_>(processed, indices.size(), std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(
        const Matrix<Value_, Index_>* mat, 
        const IndexStorage_& indices, 
        const std::vector<Index_>& remap, 
        std::integral_constant<bool, row_>, 
        std::shared_ptr<Oracle<Index_> > oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        remapping(remap) 
    {
        auto processed = format_sparse_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, block_length, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(
        const Matrix<Value_, Index_>* mat, 
        const IndexStorage_& indices, 
        const std::vector<Index_>& remap, 
        std::integral_constant<bool, row_>, 
        std::shared_ptr<Oracle<Index_> > oracle, 
        std::vector<Index_> subset, 
        const Options& opt) :
        remapping(remap)
    {
        auto processed = format_sparse_parallel<Index_>(indices, subset);
        initialize<row_>(processed, subset.size(), std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(SparseParallelResults<Index_>& processed, size_t extent, std::shared_ptr<Oracle<Index_ > > oracle, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        needs_sort = opt.sparse_ordered_index;

        if (needs_sort && needs_value) {
            sortspace.reserve(extent);
        } 
        if (!needs_index) {
            iholding.resize(processed.collapsed.size());
        }

        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(processed.collapsed), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto vinit = (needs_value ? vbuffer + shift : NULL);
        auto iinit = (needs_index ? ibuffer + shift : iholding.data());
        auto src = internal->fetch(i, vinit, iinit);
        reorder_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, needs_sort, sortspace); 
        return src;
    }

private:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
    bool needs_value, needs_index, needs_sort;
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
private:
    static void finish_assembly(
        const std::vector<std::pair<Index_, Index_> >& collected,
        const IndexStorage_& indices, 
        std::vector<Index_>& reverse_mapping,
        std::vector<Index_>& unique_and_sorted,
        Index_ mapping_dim,
        std::vector<std::pair<Index_, Index_> >& mapping_duplicates,
        std::vector<Index_>& mapping_duplicates_pool
    ) {
        unique_and_sorted.reserve(indices.size());
        reverse_mapping.resize(indices.size());

        mapping_duplicates.resize(mapping_dim);
        mapping_duplicates_pool.reserve(indices.size());

        for (Index_ i = 0, end = collected.size(); i < end; ++i) {
            const auto& current = collected[i];
            auto& range = mapping_duplicates[current.first];
            if (unique_and_sorted.empty() || current.first != unique_and_sorted.back()) {
                unique_and_sorted.push_back(current.first);
                range.first = mapping_duplicates_pool.size();
            }

            mapping_duplicates_pool.push_back(current.second);
            reverse_mapping[current.second] = unique_and_sorted.size() - 1;
            ++range.second;
        }
    }

public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * These may be duplicated and/or unsorted.
     */
    DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx) : mat(std::move(p)), indices(std::move(idx)) {
        std::vector<std::pair<Index_, Index_> > collected;
        collected.reserve(indices.size());
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            collected.emplace_back(indices[i], i);
        }
        std::sort(collected.begin(), collected.end());

        finish_assembly(
            collected,
            indices, 
            reverse_mapping,
            unique_and_sorted,
            margin_ == 0 ? mat->nrow() : mat->ncol(),
            mapping_duplicates,
            mapping_duplicates_pool
        );

        return;
    }

    /**
     * @cond
     */
    DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > p, const std::vector<std::pair<Index_, Index_> >& collected, IndexStorage_ idx) : 
        mat(std::move(p)), indices(std::move(idx)) 
    {
        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            unique_and_sorted,
            margin_ == 0 ? mat->nrow() : mat->ncol(),
            mapping_duplicates,
            mapping_duplicates_pool
        );
    }
    /**
     * @endcond
     */

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

    std::vector<Index_> reverse_mapping;
    std::vector<Index_> unique_and_sorted;
    std::vector<std::pair<Index_, Index_> > mapping_duplicates; // holds (position, size) in the pool.
    std::vector<Index_> mapping_duplicates_pool; 

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
            return std::make_unique<DelayedSubset_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubset_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubset_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubset_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...);
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
