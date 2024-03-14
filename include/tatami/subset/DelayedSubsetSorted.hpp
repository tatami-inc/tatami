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
    // 'input' and 'output' may optionally point to overlapping arrays as long
    // as 'output' precedes 'input'. The idea is that the expansion of values
    // into 'output' will cause it to "catch up" to 'input' without clobbering
    // any values in the latter. This assumes that 'input' has been shifted
    // enough to make space for expansion; the required shift depends on the
    // number of duplicates in 'expansion'.
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
struct SparseParallelExpansion {
    // Holds the start position in 'indices' for each value, i.e., 
    // 'indices[start[i - offset]]' is the first entry equal to 'i'.
    std::vector<Index_> start;

    // Holds the length of the stretch of identical values in 'indices', i.e., 
    // 'indices[start[i - offset] + lengths[i - offset] - 1]' is the last entry equal to 'i'.
    std::vector<Index_> length;

    Index_ offset = 0;
};

template<typename Index_>
struct SparseParallelResults {
    std::vector<Index_> collapsed;
    SparseParallelExpansion<Index_> expansion;
};

template<typename Index_, class IndexStorage_, class ToIndex_>
SparseParallelResults<Index_> format_sparse_parallel(const IndexStorage_& indices, Index_ len, ToIndex_ toindex) {
    SparseParallelResults<Index_> output;

    if (len) {
        output.collapsed.reserve(len);
        Index_ first = indices[toindex(0)];

        // 'start' and 'length' are vectors that enable look-up according to
        // the indices of the underlying array. To avoid the need to allocate a
        // vector of length equal to the underlying array's dimension, we only
        // consider the extremes of 'indices'; we allocate the two vectors to
        // have length equal to the range of 'indices'. The 'offset' defines
        // the lower bound that must be subtracted from the array indices to
        // get an index into 'start' or 'length'.
        output.expansion.offset = first;
        auto allocation = indices[toindex(len - 1)] - output.expansion.offset + 1;
        output.expansion.start.resize(allocation);
        output.expansion.length.resize(allocation);

        Index_ lookup = 0;
        output.expansion.start[0] = 0;
        output.expansion.length[0] = 1;
        output.collapsed.push_back(first);
        auto last = first;

        for (Index_ i = 1; i < len; ++i) {
            auto current = indices[toindex(i)];
            if (current == last) {
                ++(output.expansion.length[lookup]);
                continue;
            } 

            lookup = current - output.expansion.offset;
            output.expansion.start[lookup] = i;
            output.expansion.length[lookup] = 1;
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
SparseRange<Value_, Index_> expand_sparse_parallel(
    const SparseRange<Value_, Index_>& input, 
    Value_* vbuffer, 
    Index_* ibuffer, 
    bool needs_value,
    bool needs_index,
    const SparseParallelExpansion<Index_>& expansion)
{
    auto vcopy = vbuffer;
    auto icopy = ibuffer;
    Index_ count = 0;

    auto vsrc = input.value;
    bool replace_value = needs_value && vsrc != vcopy;

    // Pointers in 'input' and the two 'buffer' pointers may optionally point
    // to overlapping arrays as long as each 'buffer' pointer precede its
    // corresponding pointer in 'input'.  The idea is that the expansion of
    // values into 'buffer' will cause it to "catch up" to 'input' without
    // clobbering any values in the latter. This assumes that 'input' has been
    // shifted enough to make space for expansion; the required shift depends
    // on the number of duplicates.
    for (Index_ i = 0; i < input.number; ++i) {
        auto eindex = input.index[i] - expansion.offset;
        auto nexpand = expansion.length[eindex];
        count += nexpand;

        if (replace_value) {
            auto v = *vsrc; // make a copy just in case 'vcopy' and 'input.value' overlap.
            std::fill_n(vcopy, nexpand, v);
            vcopy += nexpand;
            ++vsrc;
            replace_value = (vcopy != vsrc); // if we've caught up, there no need to do this replacement.
        }

        if (needs_index) {
            auto sexpand = expansion.start[eindex];
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
        initialize<row_>(mat, std::move(processed), indices.size(), opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(mat, std::move(processed), block_length, opt);
    }

    template<bool row_, class IndexStorage_>
    MyopicParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, subset);
        initialize<row_>(mat, std::move(processed), subset.size(), opt);
    }

private:
    template<bool row_>
    void initialize(const Matrix<Value_, Index_>* mat, SparseParallelResults<Index_> processed, size_t extent, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        opt.sparse_extract_index = true; // must extract the indices for proper expansion.
        if (!needs_index) {
            iholding.reserve(processed.collapsed.size()); // need a holding space for indices if 'ibuffer' is not supplied.
        }

        internal = new_extractor<row_, true>(mat, std::move(processed.collapsed), opt);
        expansion = std::move(processed.expansion);
    }

public:
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        if (shift == 0) {
            return internal->fetch(i, vbuffer, ibuffer);
        } 

        // Shifting so that there's enough space for expansion, but only doing
        // so if these pointers are guaranteed to be non-NULL.
        auto vinit = (needs_value ? vbuffer + shift : NULL);
        auto iinit = (needs_index ? ibuffer + shift : iholding.data());
        auto src = internal->fetch(i, vinit, iinit);
        return expand_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, expansion);
    }

private:
    bool needs_value, needs_index;
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > internal;
    std::vector<Index_> iholding;
    SparseParallelExpansion<Index_> expansion;
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
        initialize<row_>(mat, std::move(processed), indices.size(), std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, block_start, block_length);
        initialize<row_>(mat, std::move(processed), block_length, std::move(oracle), opt);
    }

    template<bool row_, class IndexStorage_>
    OracularParallelSparse(const Matrix<Value_, Index_>* mat, const IndexStorage_& indices, std::integral_constant<bool, row_>, std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> subset, const Options& opt) {
        auto processed = format_sparse_parallel<Index_>(indices, subset);
        initialize<row_>(mat, std::move(processed), subset.size(), std::move(oracle), opt);
    }

private:
    template<bool row_>
    void initialize(const Matrix<Value_, Index_>* mat, SparseParallelResults<Index_> processed, size_t extent, std::shared_ptr<Oracle<Index_> > oracle, Options opt) {
        shift = extent - processed.collapsed.size();

        needs_value = opt.sparse_extract_value;
        needs_index = opt.sparse_extract_index;
        opt.sparse_extract_index = true; // must extract the indices for proper expansion.
        if (!needs_index) {
            iholding.reserve(processed.collapsed.size()); // need a holding space for indices if 'ibuffer' is not supplied.
        }

        internal = new_extractor<row_, true>(mat, std::move(oracle), std::move(processed.collapsed), opt);
        expansion = std::move(processed.expansion);
    }

public:
    SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) {
        if (shift == 0) {
            return internal->fetch(vbuffer, ibuffer);
        } 

        // Shifting so that there's enough space for expansion, but only doing
        // so if these pointers are guaranteed to be non-NULL.
        auto vinit = (needs_value ? vbuffer + shift : NULL);
        auto iinit = (needs_index ? ibuffer + shift : iholding.data());
        auto src = internal->fetch(vinit, iinit);
        return expand_sparse_parallel(src, vbuffer, ibuffer, needs_value, needs_index, expansion);
    }

private:
    bool needs_value, needs_index;
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > internal;
    std::vector<Index_> iholding;
    SparseParallelExpansion<Index_> expansion;
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
    template<bool accrow_, typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > populate_myopic_dense(Args_&& ... args) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, IndexStorage_> >(mat.get(), indices, flag, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSorted_internal::MyopicParallelDense<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetSorted_internal::MyopicParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetSorted_internal::OracularParallelDense<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
            return std::make_unique<DelayedSubsetSorted_internal::OracularParallelSparse<Value_, Index_> >(mat.get(), indices, flag, std::forward<Args_>(args)...);
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
