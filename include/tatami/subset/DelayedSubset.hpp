#ifndef TATAMI_DELAYED_SUBSET_HPP
#define TATAMI_DELAYED_SUBSET_HPP

#include "utils.hpp"
#include "../utils/Index_to_container.hpp"

#include <algorithm>
#include <memory>

#include "sanisizer/sanisizer.hpp"

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

template<typename Index_, class SubsetStorage_, class ToIndex_>
DenseParallelResults<Index_> format_dense_parallel_base(const SubsetStorage_& subset, const Index_ len, const ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(subset[to_index(i)], i);
    }
    std::sort(collected.begin(), collected.end());

    DenseParallelResults<Index_> output;
    if (collected.size()) {
        output.collapsed.reserve(len);
        resize_container_to_Index_size(output.reindex, len);

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
class ParallelDense final : public DenseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelDense(
        const Matrix<Value_, Index_>& matrix,
        const SubsetStorage_& subset,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) {
        auto processed = format_dense_parallel_base<Index_>(subset, subset.size(), [&](const Index_ i) -> Index_ { return i; });
        initialize(matrix, std::move(processed), row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelDense(
        const Matrix<Value_, Index_>& matrix,
        const SubsetStorage_& subset,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) {
        auto processed = format_dense_parallel_base<Index_>(subset, block_length, [&](const Index_ i) -> Index_ { return i + block_start; });
        initialize(matrix, std::move(processed), row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelDense(
        const Matrix<Value_, Index_>& matrix,
        const SubsetStorage_& subset,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) {
        const auto& indices = *indices_ptr;
        auto processed = format_dense_parallel_base<Index_>(subset, indices.size(), [&](const Index_ i) -> Index_ { return indices[i]; });
        initialize(matrix, std::move(processed), row, std::move(oracle), opt);
    }

private:
    void initialize(
        const Matrix<Value_, Index_>& matrix,
        DenseParallelResults<Index_> processed,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) {
        resize_container_to_Index_size(my_holding_vbuffer, processed.collapsed.size()); // processed.collapsed.size() should fit in an Index_, so this cast is safe.
        my_ext = new_extractor<false, oracle_>(matrix, row, std::move(oracle), std::move(processed.collapsed), opt);
        my_reindex.swap(processed.reindex);
    }

public:
    const Value_* fetch(const Index_ i, Value_* const buffer) {
        const auto src = my_ext->fetch(i, my_holding_vbuffer.data());

        // 'src' and 'buffer' should not point to the same array.
        auto copy = buffer;
        for (auto p : my_reindex) {
            *copy= src[p];
            ++copy;
        }

        return buffer;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;
    std::vector<Value_> my_holding_vbuffer;
    std::vector<Index_> my_reindex;
};

template<typename Index_>
struct SparseParallelReindex {
    // This is a bit complicated to explain.
    // Let 'x = pool_ptrs[i - offset]'.
    // Let 'y = pool_ptrs[i - offset + 1]'.
    // Let 'z' denote any integer in '[x, y)'.
    // In which case, 'indices[pool_indices[z]]' is equal to 'i'.
    // The general idea is that 'pool_indices[z]' can be used to fill the 'SparseRange::index' on output.
    std::vector<Index_> pool_ptrs; // this can be Index_ as the length of 'pool_indices' is no greater than the output dimension extent.
    std::vector<Index_> pool_indices;
    Index_ offset;
};

template<typename Index_>
struct SparseParallelResults {
    std::vector<Index_> collapsed;
    SparseParallelReindex<Index_> reindex;
};

template<typename Index_, class SubsetStorage_, class ToIndex_>
SparseParallelResults<Index_> format_sparse_parallel_base(const SubsetStorage_& indices, const Index_ len, const ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        const auto curdex = to_index(i);
        collected.emplace_back(indices[curdex], curdex);
    }
    std::sort(collected.begin(), collected.end());

    SparseParallelResults<Index_> output;

    if (collected.size()) {
        output.collapsed.reserve(len);
        output.reindex.pool_indices.reserve(len);
        const Index_ first = collected.front().first;

        // 'pool_ptrs' is a vector that enables look-up according to the indices of the underlying array.
        // To avoid the need to allocate a vector of length equal to the underlying array's dimension, we only consider the extremes of 'indices'.
        // We allocate 'pool_ptrs' to have length equal to the range of 'indices'... plus 1, as we're storing cumulative pointers.
        // 'offset' defines the lower bound that must be subtracted from the array indices to get an index into 'pool_ptrs'.
        output.reindex.offset = first;
        const Index_ allocation = collected.back().first - output.reindex.offset + 1;
        output.reindex.pool_ptrs.resize(sanisizer::sum<I<decltype(output.reindex.pool_ptrs.size())> >(attest_for_Index(allocation), 1));

        Index_ counter = 0; // this can never be larger than 'len', so using Index_ will not lead to overflows.
        output.reindex.pool_ptrs[counter] = 0; 
        ++counter;
        output.reindex.pool_indices.push_back(collected.front().second);
        output.reindex.pool_ptrs[counter] = 1;
        output.collapsed.push_back(first);
        auto last = first;

        for (Index_ i = 1; i < len; ++i) {
            const auto& pp = collected[i];
            const auto current = pp.first;
            if (current == last) {
                output.reindex.pool_indices.push_back(pp.second);
                ++(output.reindex.pool_ptrs[counter]); // contents of pool_ptrs will never be greater than len, so this won't overflow.
                continue;
            }

            const Index_ pool_size = output.reindex.pool_indices.size();
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
class ParallelSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelSparse(
        const Matrix<Value_, Index_>& mat,
        const SubsetStorage_& subset,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) {
        auto processed = format_sparse_parallel_base<Index_>(subset, subset.size(), [](const Index_ i) -> Index_ { return i; });
        initialize(mat, std::move(processed), subset.size(), row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelSparse(
        const Matrix<Value_, Index_>& mat,
        const SubsetStorage_& subset,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) {
        auto processed = format_sparse_parallel_base<Index_>(subset, block_length, [&](const Index_ i) -> Index_ { return i + block_start; });
        initialize(mat, std::move(processed), block_length, row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelSparse(
        const Matrix<Value_, Index_>& mat,
        const SubsetStorage_& subset,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) {
        const auto& indices = *indices_ptr;
        auto processed = format_sparse_parallel_base<Index_>(subset, indices.size(), [&](const Index_ i) -> Index_ { return indices[i]; });
        initialize(mat, std::move(processed), indices.size(), row, std::move(oracle), opt);
    }

private:
    void initialize(
        const Matrix<Value_, Index_>& mat,
        SparseParallelResults<Index_> processed,
        const Index_ extent,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Options opt
    ) {
        const Index_ num_collapsed = processed.collapsed.size(); // number of unique subset indices should be no greater than the extent.
        my_shift = extent - num_collapsed;

        my_needs_value = opt.sparse_extract_value;
        my_needs_index = opt.sparse_extract_index;
        my_needs_sort = opt.sparse_ordered_index;

        if (my_needs_sort && my_needs_value) {
            my_sortspace.reserve(extent);
        } 

        // We need to extract indices for sorting and expansion purposes, even if they weren't actually requested.
        opt.sparse_extract_index = true;
        if (!my_needs_index) {
            resize_container_to_Index_size(my_holding_ibuffer, num_collapsed);
        }

        my_ext = new_extractor<true, oracle_>(mat, row, std::move(oracle), std::move(processed.collapsed), opt);
        my_reindex = std::move(processed.reindex);
    }

public:
    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const vbuffer, Index_* const ibuffer) {
        const auto vinit = (my_needs_value ? vbuffer + my_shift : NULL);
        const auto iinit = (my_needs_index ? ibuffer + my_shift : my_holding_ibuffer.data());
        auto input = my_ext->fetch(i, vinit, iinit);

        if (!my_needs_sort) {
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
            bool replace_value = my_needs_value && vsrc != vcopy;

            for (Index_ i = 0; i < input.number; ++i) {
                const auto lookup = input.index[i] - my_reindex.offset;
                const auto start = my_reindex.pool_ptrs[lookup];
                const auto num = my_reindex.pool_ptrs[lookup + 1] - start;
                count += num;

                if (replace_value) {
                    auto val = *vsrc; // make a copy just in case 'vcopy' and 'input.value' overlap.
                    std::fill_n(vcopy, num, val);
                    vcopy += num;
                    ++vsrc;
                    replace_value = (vcopy != vsrc); // if we've caught up, there no need to do this replacement.
                }

                if (my_needs_index) {
                    // Again, 'icopy' will eventually catch up to 'input.index' if
                    // they point to overlapping arrays. But we still need to
                    // replace values once we've managed to catch up, so we can't
                    // short-circuit like we did with 'replace_value'.
                    std::copy_n(my_reindex.pool_indices.begin() + start, num, icopy);
                    icopy += num;
                }
            }

            input.number = count;
            if (my_needs_value) {
                input.value = vbuffer;
            }
            if (my_needs_index) {
                input.index = ibuffer;
            } else {
                input.index = NULL;
            }

        } else if (my_needs_value) {
            // This does not require any careful consideration of the overlaps
            // between 'input' and 'buffers', as we're copying things into
            // 'my_sortspace' anyway before copying them back into 'buffer'.
            my_sortspace.clear();
            for (Index_ i = 0; i < input.number; ++i) {
                const auto val = input.value[i];
                const auto lookup = input.index[i] - my_reindex.offset;
                const auto start = my_reindex.pool_ptrs[lookup];
                const auto end = my_reindex.pool_ptrs[lookup + 1];
                for (Index_ j = start; j < end; ++j) {
                    my_sortspace.emplace_back(my_reindex.pool_indices[j], val);
                }
            }
            std::sort(my_sortspace.begin(), my_sortspace.end());
            input.number = my_sortspace.size();

            auto vcopy = vbuffer;
            for (const auto& ss : my_sortspace) {
                *vcopy = ss.second;
                ++vcopy;
            }
            input.value = vbuffer;

            if (my_needs_index) {
                auto icopy = ibuffer;
                for (const auto& ss : my_sortspace) {
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
                const auto lookup = input.index[i] - my_reindex.offset;
                const auto start = my_reindex.pool_ptrs[lookup];
                const auto num = my_reindex.pool_ptrs[lookup + 1] - start;
                count += num;

                if (my_needs_index) {
                    std::copy_n(my_reindex.pool_indices.begin() + start, num, icopy);
                    icopy += num;
                }
            }

            input.number = count;
            if (my_needs_index) {
                std::sort(ibuffer, ibuffer + count);
                input.index = ibuffer;
            } else {
                input.index = NULL;
            }
        }

        return input;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;
    bool my_needs_value, my_needs_index, my_needs_sort;
    SparseParallelReindex<Index_> my_reindex;
    std::vector<std::pair<Index_, Value_> > my_sortspace;
    std::vector<Index_> my_holding_ibuffer;
    Index_ my_shift;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting of a matrix with general indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of arbitrary indices.
 * This operation is "delayed" in that it is only evaluated when data is extracted from the matrix.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam SubsetStorage_ Vector containing the subset indices.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 */
template<typename Value_, typename Index_, class SubsetStorage_>
class DelayedSubset final : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying (pre-subset) matrix.
     * @param subset Vector of 0-based indices to use for subsetting on the rows (if `by_row = true`) or columns (otherwise).
     * These may be duplicated and/or unsorted.
     * @param by_row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     */
    DelayedSubset(
        std::shared_ptr<const Matrix<Value_, Index_> > matrix,
        SubsetStorage_ subset,
        const bool by_row
    ) : 
        my_matrix(std::move(matrix)),
        my_subset(std::move(subset)),
        my_by_row(by_row)
    {
        // Check that we can still report the dimension extents of the subsetted matrix.
        sanisizer::can_cast<Index_>(my_subset.size());
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > my_matrix;
    SubsetStorage_ my_subset;
    bool my_by_row;

public:
    Index_ nrow() const {
        if (my_by_row) {
            return my_subset.size();
        } else {
            return my_matrix->nrow();
        }
    }

    Index_ ncol() const {
        if (my_by_row) {
            return my_matrix->ncol();
        } else {
            return my_subset.size();
        }
    }

    bool is_sparse() const {
        return my_matrix->is_sparse();
    }

    double is_sparse_proportion() const {
        return my_matrix->is_sparse_proportion();
    }

    bool prefer_rows() const {
        return my_matrix->prefer_rows();
    }

    double prefer_rows_proportion() const {
        return my_matrix->prefer_rows_proportion();
    }

    bool uses_oracle(const bool row) const {
        return my_matrix->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense;

    using Matrix<Value_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<typename ... Args_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > populate_myopic_dense(
        const bool row,
        Args_&& ... args
    ) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, SubsetStorage_> >(
                *my_matrix,
                my_subset,
                row,
                std::forward<Args_>(args)...
            ); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelDense<false, Value_, Index_> >(
                *my_matrix,
                my_subset,
                row,
                false,
                std::forward<Args_>(args)...
            );
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        const Options& opt
    ) const {
        return populate_myopic_dense(row, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return populate_myopic_dense(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        VectorPtr<Index_> my_subset_ptr,
        const Options& opt
    ) const {
        return populate_myopic_dense(row, std::move(my_subset_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    template<typename ... Args_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > populate_myopic_sparse(
        const bool row,
        Args_&& ... args
    ) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularSparse<Value_, Index_, SubsetStorage_> >(
                *my_matrix,
                my_subset,
                row,
                std::forward<Args_>(args)...
            ); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelSparse<false, Value_, Index_> >(
                *my_matrix,
                my_subset,
                row,
                false,
                std::forward<Args_>(args)...
            );
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        const Options& opt
    ) const {
        return populate_myopic_sparse(row, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return populate_myopic_sparse(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        VectorPtr<Index_> my_subset_ptr,
        const Options& opt
    ) const {
        return populate_myopic_sparse(row, std::move(my_subset_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
private:
    template<typename ... Args_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > populate_oracular_dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        Args_&& ... args
    ) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::OracularPerpendicularDense<Value_, Index_> >(
                *my_matrix,
                my_subset,
                row,
                std::move(oracle),
                std::forward<Args_>(args)...
            ); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelDense<true, Value_, Index_> >(
                *my_matrix,
                my_subset,
                row,
                std::move(oracle),
                std::forward<Args_>(args)...
            );
        }
    }

public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        return populate_oracular_dense(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return populate_oracular_dense(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> my_subset_ptr,
        const Options& opt
    ) const {
        return populate_oracular_dense(row, std::move(oracle), std::move(my_subset_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
private:
    template<typename ... Args_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > populate_oracular_sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        Args_&& ... args
    ) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::OracularPerpendicularSparse<Value_, Index_> >(
                *my_matrix,
                my_subset,
                row,
                std::move(oracle),
                std::forward<Args_>(args)...
            ); 
        } else {
            return std::make_unique<DelayedSubset_internal::ParallelSparse<true, Value_, Index_> >(
                *my_matrix,
                my_subset,
                row,
                std::move(oracle),
                std::forward<Args_>(args)...
            );
        }
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        return populate_oracular_sparse(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        return populate_oracular_sparse(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> my_subset_ptr,
        const Options& opt
    ) const {
        return populate_oracular_sparse(row, std::move(oracle), std::move(my_subset_ptr), opt);
    }
};

}

#endif
