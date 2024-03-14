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

/********************
 *** Myopic dense ***
 ********************/

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
void reorder_dense_parallel(const Value_* input, Value_* output, const std::vector<Index_>& permutation) {
    // 'input' and 'output' should not point to the same array. In theory, it
    // is possible to do an in-place permutation, but this requires another
    // array anyway to track the permutation status, so we'll just keep it simple.
    for (auto p : permutation) {
        output[p] = *input;
        ++input;
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelDense : MyopicDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices);
        initialize<row_>(mat, std::move(processed), indices.size(), opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(mat, std::move(processed), block_length, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, subset);
        initialize<row_>(mat, std::move(processed), subset.size(), opt);
    }

private:
    template<bool row_>
    void initialize(const Matrix<Value_, Index_>* mat, DenseParallelResults<Index_> processed, size_t extent, const Options& opt) {
        internal = new_extractor<row_, false>(mat, std::move(processed.sorted), opt);
        vholding.resize(extent);
        permutation = std::move(processed.permutation);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto src = internal->fetch(i, vholding.data());
        reorder_dense_parallel(src, buffer, permutation);
        return buffer;
    }

private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
    std::vector<Value_> vholding;
    std::vector<Index_> permutation;
};

/**********************
 *** Oracular dense ***
 **********************/

template<typename Value_, typename Index_>
struct OracularParallelDense : OracularDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices);
        initialize<row_>(mat, std::move(processed), indices.size(), std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(mat, std::move(processed), block_length, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, subset);
        initialize<row_>(mat, std::move(processed), subset.size(), std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(const Matrix<Value_, Index_>* mat, DenseParallelResults<Index_> processed, size_t extent, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) { 
        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(processed.sorted), opt);
        vholding.resize(extent);
        permutation = std::move(processed.permutation);
    }

public:
    const Value_* fetch(Value_* buffer) {
        auto src = internal->fetch(vholding.data());
        reorder_dense_parallel(src, buffer, permutation);
        return buffer;
    }

private:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
    std::vector<Value_> vholding;
    std::vector<Index_> permutation;
};

/*********************
 *** Myopic sparse ***
 *********************/

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

template<typename Index_, class IndexStorage_>
std::vector<Index_> format_sparse_parallel(const IndexStorage_& indices) {
    return format_sparse_parallel<Index_>(indices, indices.size(), [&](Index_ i) -> Index_ { return i; });
}

template<typename Index_, class IndexStorage_>
std::vector<Index_> format_sparse_parallel(const IndexStorage_& indices, Index_ start, Index_ length) {
    return format_sparse_parallel<Index_>(indices, length, [&](Index_ i) -> Index_ { return i + start; });
}

template<typename Index_, class IndexStorage_>
std::vector<Index_> format_sparse_parallel(const IndexStorage_& indices, const std::vector<Index_>& subset) {
    return format_sparse_parallel<Index_>(indices, subset.size(), [&](Index_ i) -> Index_ { return subset[i]; });
}

template<typename Value_, typename Index_>
void reorder_sparse_parallel(
    SparseRange<Value_, Index_>& input, 
    Value_* vbuffer, 
    Index_* ibuffer,
    bool needs_value,
    bool needs_index,
    bool needs_sort,
    std::vector<std::pair<Index_, Value_> >& sortspace,
    const std::vector<Index_>& remapping) 
{
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
}

template<typename Value_, typename Index_>
void allocate_sparse_parallel(
    bool needs_value,
    bool needs_index,
    bool needs_sort,
    std::vector<std::pair<Index_, Value_> >& sortspace,
    std::vector<Index_>& iholding,
    Options& opt,
    size_t extent)
{
    // The conditionals here mirror those in 'reorder_sparse_parallel',
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
}

template<typename Value_, typename Index_>
struct MyopicParallelSparse : MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(
        const Matrix<Value_, Index_>* mat, 
        const IndexStorage_& indices, 
        const std::vector<Index_>& remap, 
        std::integral_constant<bool, row_>, 
        const Options& opt) : 
        remapping(remap) 
    {
        auto processed = format_sparse_parallel<Index_>(indices);
        initialize<row_>(mat, std::move(processed), indices.size(), opt);
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
        initialize<row_>(mat, std::move(processed), block_length, opt);
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
        initialize<row_>(mat, std::move(processed), subset.size(), opt);
    }

private:
    template<bool row_>
    void initialize(const Matrix<Value_, Index_>* mat, std::vector<Index_> sorted, size_t extent, Options opt) {
        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        needs_sort = opt.sparse_ordered_index;
        allocate_sparse_parallel(needs_value, needs_index, needs_sort, sortspace, iholding, opt, extent);
        internal = new_extractor<row_, true>(mat, std::move(sorted), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto src = internal->fetch(i, vbuffer, (iholding.empty() ? ibuffer : iholding.data()));
        reorder_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, needs_sort, sortspace, remapping); 
        return src;
    }

private:
    const std::vector<Index_>& remapping;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    bool needs_value, needs_index, needs_sort;
    std::vector<std::pair<Index_, Value_> > sortspace;
    std::vector<Index_> iholding;
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
        const std::vector<Index_>& remap, 
        std::integral_constant<bool, row_>, 
        std::shared_ptr<Oracle<Index_> > oracle,
        const Options& opt) : 
        remapping(remap) 
    {
        auto processed = format_sparse_parallel<Index_>(indices);
        initialize<row_>(mat, std::move(processed), indices.size(), std::move(oracle), opt);
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
        initialize<row_>(mat, std::move(processed), block_length, std::move(oracle), opt);
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
        initialize<row_>(mat, std::move(processed), subset.size(), std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(const Matrix<Value_, Index_>* mat, std::vector<Index_> sorted, size_t extent, std::shared_ptr<Oracle<Index_> > oracle, Options opt) {
        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        needs_sort = opt.sparse_ordered_index;
        allocate_sparse_parallel(needs_value, needs_index, needs_sort, sortspace, iholding, opt, extent);
        internal = new_extractor<row_, true>(mat, std::move(oracle), std::move(sorted), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        auto src = internal->fetch(vbuffer, (iholding.empty() ? ibuffer : iholding.data()));
        reorder_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, needs_sort, sortspace, remapping); 
        return src;
    }

private:
    const std::vector<Index_>& remapping;
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
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
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetUnique : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * This should be unique, but may be unsorted.
     * @param check Whether to check `idx` for unique values.
     */
    DelayedSubsetUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        Index_ fulldim = margin_ == 0 ? mat->nrow() : mat->ncol();

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
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > populate_myopic_dense(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetUnique_internal::MyopicParallelDense<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
        }
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
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > populate_myopic_sparse(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularSparse<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetUnique_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, mapping_single, flag, std::forward<Args_>(args)...);
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
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > populate_oracular_dense(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::OracularPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetUnique_internal::OracularParallelDense<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
        }
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
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > populate_oracular_sparse(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::OracularPerpendicularSparse<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetUnique_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, mapping_single, flag, std::forward<Args_>(args)...);
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
