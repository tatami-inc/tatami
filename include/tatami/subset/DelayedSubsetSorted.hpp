#ifndef TATAMI_DELAYED_SUBSET_SORTED_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_HPP

#include "utils.hpp"
#include "../base/Matrix.hpp"

#include <algorithm>
#include <numeric>
#include <memory>

/**
 * @file DelayedSubsetSorted.hpp
 *
 * @brief Delayed subsetting with sorted row/column indices.
 *
 * This is equivalent to the `DelayedSubset` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedSubsetSorted_internal {

/********************
 *** Myopic dense ***
 ********************/

template<typename Index_>
struct DenseParallelResults {
    std::vector<Index_> collapsed;
    std::vector<Index_> expansion;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
DenseParallelResults<Index_> format_dense_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ to_index) {
    DenseParallelResults<Index_> output;
    output.expansion.reserve(len);
    output.collapsed.reserve(len);

    if (len) {
        Index_ last = indices[to_index(0)];
        output.expansion.push_back(1);
        output.collapsed.push_back(last);

        for (Index_ i = 1; i < len; ++i) {
            auto current = indices[to_index(i)];
            if (current == last) {
                ++(output.expansion.back());
            } else {
                last = current;
                output.expansion.push_back(1);
                output.collapsed.push_back(last);
            }
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
void expand_dense_parallel(const Value_* input, Value_* output, const std::vector<Index_>& expansion) {
    // 'input' and 'output' may point to overlapping arrays as long as 'output'
    // precedes 'input'. The idea is that the expansion of values into 'output'
    // will cause it to "catch up" to 'input' without clobbering any values in
    // the latter. This assumes that 'input' has been shifted enough to make
    // space for expansion; the required shift depends on the number of
    // duplicates in 'expansion'.
    for (auto e : expansion) {
        auto val = *input;
        std::fill_n(output, e, val);
        ++input;
        output += e;

        // Once we've caught up, everything else must be a non-duplicate,
        // otherwise we'd be clobbering as-yet-unread values from the input.
        // So we might as well just quit at this point.
        if (input == output) {
            return;
        }
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelDense : MyopicDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices);
        initialize<row_>(processed, indices.size(), opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, block_length, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, subset);
        initialize<row_>(processed, subset.size(), opt);
    }

private:
    template<bool row_>
    void initialize(DenseParallelResults<Index_>& processed, size_t extent, const Options& opt) {
        shift = extent - processed.collapsed.size();
        internal = new_extractor<row_, false>(mat, std::move(processed.collapsed), opt);
        expansion = std::move(processed.expansion);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        if (shift == 0) {
            return internal->fetch(i, buffer);
        } 
        // Shifting so that there's enough space for expansion.
        auto src = internal->fetch(i, buffer + shift);
        expand_dense_parallel(src, buffer, expansion);
        return buffer;
    }

private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > internal;
    std::vector<Index_> expansion;
    size_t shift;
};

/**********************
 *** Oracular dense ***
 **********************/

template<typename Value_, typename Index_>
struct OracularParallelDense : OracularDenseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices);
        initialize<row_>(processed, indices.size(), std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, block_length, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelDense(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(indices, subset);
        initialize<row_>(processed, subset.size(), std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(DenseParallelResults<Index_>& processed, size_t extent, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) { 
        shift = extent - processed.collapsed.size();
        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(processed.collapsed), opt);
        expansion = std::move(processed.expansion);
    }

public:
    const Value_* fetch(Value_* buffer) {
        if (shift == 0) {
            return internal->fetch(buffer);
        }
        // Shifting so that there's enough space for expansion.
        auto src = internal->fetch(buffer + shift);
        expand_dense_parallel(src, buffer, expansion);
        return buffer;
    }

private:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > internal;
    std::vector<Index_> expansion;
    size_t shift;
};

/*********************
 *** Myopic sparse ***
 *********************/

template<typename Index_>
struct SparseParallelResults {
    std::vector<Index_> collapsed;
    std::vector<Index_> expansion_start, expansion_length;
    Index_ expansion_offset = 0;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
SparseParallelResults<Index_> format_sparse_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ toindex) {
    SparseParallelResults<Index_> output;

    if (len) {
        auto last = indices[toindex(0)];

        auto allocation = indices[toindex(len - 1)] - last + 1;
        output.expansion_start.resize(allocation);
        output.expansion_length.resize(allocation);
        output.collapsed.reserve(len);

        output.collapsed.push_back(last);
        output.expansion_start[0] = 0;
        output.expansion_length[0] = 1;

        for (Index_ i = 1; i < len; ++i) {
            auto current = indices[toindex(i)];
            if (current == last) {
                ++(output.expansion_length[current - start]);
            } else {
                last = current;
                output.expansion_start[current - start] = i;
                output.expansion_length[current - start] = 1;
                output.collapsed.push_back(last);
            }
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
SparseRange<Value_, Index_> expand_sparse_parallel(
    const SparseRange<Value_, Index_>& input, 
    Value_* vbuffer, 
    Index_* ibuffer, 
    bool needs_value,
    bool needs_index,
    const std::vector<Index_>& expansion_start, 
    const std::vector<Index_>& expansion_length,
    Index_ expansion_offset) 
{
    auto vcopy = vbuffer;
    auto icopy = ibuffer;
    Index_ count = 0;
    bool replace_value = needs_value;

    // Pointers in 'input' and the two 'buffer' pointers may point to
    // overlapping arrays as long as both 'buffer' pointers precede 'input'.
    // The idea is that the expansion of values into 'buffer' will cause it to
    // "catch up" to 'input' without clobbering any values in the latter. This
    // assumes that 'input' has been shifted enough to make space for
    // expansion; the required shift depends on the number of duplicates.
    for (Index i = 0; i < input.number; ++i) {
        auto expansion_index = input.index[i] - expansion_offset;
        auto nexpand = expansion_length[expansion_index];
        count += nexpand;

        if (replace_value) {
            auto v = input.value[i];
            std::fill_n(vcopy, nexpand, v);
            vcopy += nexpand;
            replace_value = (vcopy != input.value); // if we've caught up, there no need to do this replacement.
        }

        if (needs_index) {
            auto sexpand = expansion_start[expansion_index];
            std::iota(icopy, icopy + nexpand, sexpand);
            icopy += nexpand;
        }
    }

    return SparseRange<Value_, Index_>(
        count, 
        (needs_value ? vbuffer : NULL),
        (needs_index ? ibuffer : NULL)
    );
}

template<typename Value_, typename Index_>
struct MyopicParallelSparse : MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices);
        initialize<row_>(processed, indices.size(), opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, block_length, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, subset);
        initialize<row_>(processed, subset.size(), opt);
    }

private:
    template<bool row_>
    void initialize(SparseParallelResults<Index_>& processed, size_t extent, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        opt.sparse_extract_index = true; // must extract the indices for proper expansion.
        if (!needs_index) {
            iholding.reserve(processed.collapsed.size()); // need a holding space for indices if 'ibuffer' is not supplied.
        }

        internal = new_extractor<row_, false>(mat, std::move(processed.collapsed), opt);
        expansion_start = std::move(processed.expansion_start);
        expansion_length = std::move(processed.expansion_length);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        if (shift == 0) {
            return internal->fetch(i, vbuffer, ibuffer);
        } 
        // Shifting so that there's enough space for expansion, but only doing
        // so if we actually are guaranteed non-NULL pointers.
        auto src = internal->fetch(i, (needs_value ? vbuffer + shift : NULL), (needs_index ? ibuffer + shift : iholding.data()));
        return expand_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, expansion_start, expansion_length);
    }

private:
    bool needs_value, needs_index;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    std::vector<Index_> iholding;
    std::vector<Index_> expansion, expansion_length;
    size_t shift;
};

/***********************
 *** Oracular sparse ***
 ***********************/

template<typename Value_, typename Index_>
struct OracularParallelSparse : OracularSparseExtractor<Value_, Index_> {
    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices);
        initialize<row_>(processed, indices.size(), std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(processed, block_length, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, subset);
        initialize<row_>(processed, subset.size(), std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(SparseParallelResults<Index_>& processed, size_t extent, std::shared_ptr<Oracle<Index_> > oracle, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        opt.sparse_extract_index = true; // must extract the indices for proper expansion.
        if (!needs_index) {
            iholding.reserve(processed.collapsed.size()); // need a holding space for indices if 'ibuffer' is not supplied.
        }

        internal = new_extractor<row_, false>(mat, std::move(oracle), std::move(processed.collapsed), opt);
        expansion_start = std::move(processed.expansion_start);
        expansion_length = std::move(processed.expansion_length);
    }

public:
    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        if (shift == 0) {
            return internal->fetch(vbuffer, ibuffer);
        } 
        // Shifting so that there's enough space for expansion, but only doing
        // so if we actually are guaranteed non-NULL pointers.
        auto src = internal->fetch((needs_value ? vbuffer + shift : NULL), (needs_index ? ibuffer + shift : iholding.data()));
        return expand_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, expansion_start, expansion_length);
    }

private:
    bool needs_value, needs_index;
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
    std::vector<Index_> iholding;
    std::vector<Index_> expansion, expansion_length;
    size_t shift;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting of a matrix with sorted indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of sorted indices.
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
class DelayedSubsetSorted : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * This should be sorted, but may be duplicated.
     * @param check Whether to check `idx` for sorted values.
     */
    DelayedSubsetSorted(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        if (check) {
            for (Index_ i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] < indices[i-1]) {
                    throw std::runtime_error("indices should be sorted");
                }
            }
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

    Index_ get_mapping_dim() const {
        if constexpr(margin_ == 0) {
            return mat->nrow();
        } else {
            return mat->ncol();
        }
    }

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
            return std::make_unique<DelayedSubsetSorted_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetSorted_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetSorted_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetSorted_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, remapping, flag, std::forward<Args_>(args)...);
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
