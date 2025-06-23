#ifndef TATAMI_DELAYED_SUBSET_SORTED_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_HPP

#include "utils.hpp"
#include "../base/Matrix.hpp"
#include "../utils/Index_to_container.hpp"

#include <algorithm>
#include <numeric>
#include <memory>

#include "sanisizer/sanisizer.hpp"

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

template<typename Index_>
struct DenseParallelResults {
    std::vector<Index_> collapsed;
    std::vector<Index_> expansion;
};

template<typename Index_, class SubsetStorage_, class ToIndex_>
DenseParallelResults<Index_> format_dense_parallel(const SubsetStorage_& indices, Index_ len, ToIndex_ to_index) {
    DenseParallelResults<Index_> output;
    output.expansion.reserve(len);
    output.collapsed.reserve(len);

    if (len) {
        auto last = indices[to_index(0)];
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

template<bool oracle_, typename Value_, typename Index_>
class ParallelDense final : public DenseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelDense(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(subset, subset.size(), [&](Index_ i) -> Index_ { return i; });
        initialize(matrix, std::move(processed), subset.size(), row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelDense(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) {
        auto processed = format_dense_parallel<Index_>(subset, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
        initialize(matrix, std::move(processed), block_length, row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelDense(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) {
        const auto& indices = *indices_ptr;
        auto processed = format_dense_parallel<Index_>(subset, indices.size(), [&](Index_ i) -> Index_ { return indices[i]; });
        initialize(matrix, std::move(processed), indices.size(), row, std::move(oracle), opt);
    }

private:
    void initialize(const Matrix<Value_, Index_>* mat, DenseParallelResults<Index_> processed, Index_ extent, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) {
        my_shift = extent - processed.collapsed.size();
        my_ext = new_extractor<false, oracle_>(mat, row, std::move(oracle), std::move(processed.collapsed), opt);
        my_expansion = std::move(processed.expansion);
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto input = my_ext->fetch(i, buffer + my_shift);

        // 'input' and 'buffer' may optionally point to overlapping arrays as long
        // as 'buffer' precedes 'input'. The idea is that the expansion of values
        // into 'buffer' will cause it to "catch up" to 'input' without clobbering
        // any values in the latter. This assumes that 'input' has been shifted
        // enough to make space for expansion; the required shift depends on the
        // number of duplicates in 'expansion'.

        auto copy = buffer;
        for (auto e : my_expansion) {
            // Once we've caught up, everything else must be a non-duplicate,
            // otherwise we'd be clobbering as-yet-unread values from the input.
            // So we might as well just quit at this point.
            if (input == copy) {
                break;
            }

            auto val = *input;
            std::fill_n(copy, e, val);
            ++input;
            copy += e;
        }

        return buffer;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;
    std::vector<Index_> my_expansion;
    Index_ my_shift;
};

template<typename Index_>
struct SparseParallelExpansion {
    // This is a bit complicated to explain.
    // Let 'x = start[i - offset]'.
    // Let 'y = lengths[i - offset]'.
    // Let 'z' denote any integer in '[x, x + y)'.
    // Let 'f' be the selection-specific function such that 'f(a)' is the a-th element of the selection
    // (i.e., 'a' for full selection, 'a + start' for block selection and 'subset[a]' for indexed selection).
    // In which case, 'indices[f(z)]' is equal to 'i'.
    // The general idea is that 'f(z)' can be used to fill the 'SparseRange::index' on output.
    std::vector<Index_> start;
    std::vector<Index_> length;

    Index_ offset = 0;
};

template<typename Index_>
struct SparseParallelResults {
    std::vector<Index_> collapsed;
    SparseParallelExpansion<Index_> expansion;
};

template<typename Index_, class SubsetStorage_, class ToIndex_>
SparseParallelResults<Index_> format_sparse_parallel(const SubsetStorage_& indices, Index_ len, ToIndex_ to_index) {
    SparseParallelResults<Index_> output;

    if (len) {
        output.collapsed.reserve(len);
        auto first = indices[to_index(0)];

        // 'start' and 'length' are vectors that enable look-up according to the indices of the underlying array.
        // To avoid the need to allocate a vector of length equal to the underlying array's dimension, we only consider the extremes of 'indices'.
        // We allocate both start/length vectors to have length equal to the range of 'indices'.
        // The 'offset' defines the lower bound that must be subtracted from the array indices to get an index into 'start' or 'length'.
        output.expansion.offset = first;
        Index_ allocation = indices[to_index(len - 1)] - output.expansion.offset + 1; // fits in an Index_ as it should be no greater than the input matrix's dimension extent.
        resize_container_to_Index_size(output.expansion.start, allocation);
        output.expansion.length.resize(output.expansion.start.size());

        Index_ lookup = 0;
        output.expansion.start[0] = 0;
        output.expansion.length[0] = 1;
        output.collapsed.push_back(first);
        auto last = first;

        for (Index_ i = 1; i < len; ++i) {
            auto current = indices[to_index(i)];
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

template<bool oracle_, typename Value_, typename Index_>
class ParallelSparseCore {
public:
    template<class SubsetStorage_, class ToIndex_>
    ParallelSparseCore(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, Index_ extent, bool row, MaybeOracle<oracle_, Index_> oracle, Options opt, ToIndex_ to_index) {
        auto processed = format_sparse_parallel<Index_>(subset, extent, std::forward<ToIndex_>(to_index));
        my_shift = extent - processed.collapsed.size();

        my_needs_value = opt.sparse_extract_value;
        my_needs_index = opt.sparse_extract_index;
        opt.sparse_extract_index = true; // must extract the indices for proper my_expansion.
        if (!my_needs_index) {
            // Need a holding space for indices if 'ibuffer' is not supplied.
            resize_container_to_Index_size(my_holding_ibuffer, processed.collapsed.size()); // processed.collapsed.size() should fit in an Index_, hence the cast is safe.
        }

        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::move(processed.collapsed), opt);
        my_expansion = std::move(processed.expansion);
    }

    template<class ToIndex_>
    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer, ToIndex_ to_index) {
        // Shifting so that there's enough space for my_expansion, but only doing
        // so if these pointers are guaranteed to be non-NULL.
        auto vinit = (my_needs_value ? value_buffer + my_shift : NULL);
        auto iinit = (my_needs_index ? index_buffer + my_shift : my_holding_ibuffer.data());
        auto input = my_ext->fetch(i, vinit, iinit);

        auto vcopy = value_buffer;
        auto icopy = index_buffer;
        Index_ count = 0;

        auto vsrc = input.value;
        bool replace_value = my_needs_value && vsrc != vcopy;

        // Pointers in 'input' and the two 'buffer' pointers may optionally point
        // to overlapping arrays as long as each 'buffer' pointer precede its
        // corresponding pointer in 'input'.  The idea is that the expansion of
        // values into 'buffer' will cause it to "catch up" to 'input' without
        // clobbering any values in the latter. This assumes that 'input' has been
        // shift enough to make space for expansion; the required shift depends
        // on the number of duplicates.
        for (Index_ i = 0; i < input.number; ++i) {
            auto eindex = input.index[i] - my_expansion.offset;
            auto nexpand = my_expansion.length[eindex];
            count += nexpand;

            if (replace_value) {
                auto v = *vsrc; // make a copy just in case 'vcopy' and 'input.value' overlap.
                std::fill_n(vcopy, nexpand, v);
                vcopy += nexpand;
                ++vsrc;
                replace_value = (vcopy != vsrc); // if we've caught up, there no need to do this replacement.
            }

            if (my_needs_index) {
                auto sexpand = my_expansion.start[eindex];
                for (Index_ e = 0; e < nexpand; ++e, ++icopy) {
                    *icopy = to_index(sexpand + e);
                }
            }
        }

        return SparseRange<Value_, Index_>(
            count, 
            (my_needs_value ? value_buffer : NULL),
            (my_needs_index ? index_buffer : NULL)
        );
    }

private:
    bool my_needs_value, my_needs_index;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;
    std::vector<Index_> my_holding_ibuffer;
    SparseParallelExpansion<Index_> my_expansion;
    Index_ my_shift;
};

template<bool oracle_, typename Value_, typename Index_>
class ParallelFullSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelFullSparse(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, MaybeOracle<oracle_, Index_> oracle, const Options& opt) :
        my_core(matrix, subset, subset.size(), row, std::move(oracle), opt, [](Index_ i) -> Index_ { return i; }) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        return my_core.fetch(i, value_buffer, index_buffer, [](Index_ i) -> Index_ { return i; });
    }

private:
    ParallelSparseCore<oracle_, Value_, Index_> my_core;
};

template<bool oracle_, typename Value_, typename Index_>
class ParallelBlockSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelBlockSparse(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, MaybeOracle<oracle_, Index_> oracle, Index_ block_start, Index_ block_length, const Options& opt) : 
        my_core(matrix, subset, block_length, row, std::move(oracle), opt, [&](Index_ i) -> Index_ { return i + block_start; }),
        my_block_start(block_start)
    {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        return my_core.fetch(i, value_buffer, index_buffer, [&](Index_ i) -> Index_ { return i + my_block_start; });
    }

private:
    ParallelSparseCore<oracle_, Value_, Index_> my_core;
    Index_ my_block_start;
};

template<bool oracle_, typename Value_, typename Index_>
class ParallelIndexSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelIndexSparse(const Matrix<Value_, Index_>* matrix, const SubsetStorage_& subset, bool row, MaybeOracle<oracle_, Index_> oracle, VectorPtr<Index_> indices_ptr, const Options& opt) : 
        my_core(matrix, subset, indices_ptr->size(), row, std::move(oracle), opt, [&](Index_ i) -> Index_ { return indices_ptr->operator[](i); }),
        my_indices_ptr(std::move(indices_ptr))
    {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        const auto& indices = *my_indices_ptr;
        return my_core.fetch(i, value_buffer, index_buffer, [&](Index_ i) -> Index_ { return indices[i]; });
    }

private:
    ParallelSparseCore<oracle_, Value_, Index_> my_core;
    VectorPtr<Index_> my_indices_ptr;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting of a matrix with sorted indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of sorted indices.
 * This operation is "delayed" in that it is only evaluated when data is extracted from the matrix.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam SubsetStorage_ Vector containing the subset indices.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 */
template<typename Value_, typename Index_, class SubsetStorage_>
class DelayedSubsetSorted final : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying (pre-subset) matrix.
     * @param subset Vector of 0-based indices to use for subsetting on the rows (if `by_row = true`) or columns (otherwise).
     * This should be sorted, but may be duplicated.
     * @param by_row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     * @param check Whether to check `subset` for sorted values.
     */
    DelayedSubsetSorted(std::shared_ptr<const Matrix<Value_, Index_> > matrix, SubsetStorage_ subset, bool by_row, bool check = true) : 
        my_matrix(std::move(matrix)), my_subset(std::move(subset)), my_by_row(by_row) 
    {
        sanisizer::can_cast<Index_>(my_subset.size());
        if (check) {
            for (Index_ i = 1, end = my_subset.size(); i < end; ++i) {
                if (my_subset[i] < my_subset[i-1]) {
                    throw std::runtime_error("my_subset should be sorted");
                }
            }
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > my_matrix;
    SubsetStorage_ my_subset;
    bool my_by_row;

    Index_ get_mapping_dim() const {
        if (my_by_row) {
            return my_matrix->nrow();
        } else {
            return my_matrix->ncol();
        }
    }

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

    bool uses_oracle(bool row) const {
        return my_matrix->uses_oracle(row);
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
        if (row == my_by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularDense<Value_, Index_, SubsetStorage_> >(my_matrix.get(), my_subset, row, std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSorted_internal::ParallelDense<false, Value_, Index_> >(my_matrix.get(), my_subset, row, false, std::forward<Args_>(args)...);
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
    template<DimensionSelectionType selection_, bool oracle_, typename ... Args_>
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > populate_sparse(bool row, MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) const {
        if constexpr(selection_ == DimensionSelectionType::FULL) {
            return std::make_unique<DelayedSubsetSorted_internal::ParallelFullSparse<oracle_, Value_, Index_> >(my_matrix.get(), my_subset, row, std::move(oracle), std::forward<Args_>(args)...);
        } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
            return std::make_unique<DelayedSubsetSorted_internal::ParallelBlockSparse<oracle_, Value_, Index_> >(my_matrix.get(), my_subset, row, std::move(oracle), std::forward<Args_>(args)...);
        } else {
            return std::make_unique<DelayedSubsetSorted_internal::ParallelIndexSparse<oracle_, Value_, Index_> >(my_matrix.get(), my_subset, row, std::move(oracle), std::forward<Args_>(args)...);
        }
    }

    template<DimensionSelectionType selection_, typename ... Args_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > populate_myopic_sparse(bool row, Args_&& ... args) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::MyopicPerpendicularSparse<Value_, Index_, SubsetStorage_> >(my_matrix.get(), my_subset, row, std::forward<Args_>(args)...); 
        } else {
            return populate_sparse<selection_, false>(row, false, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        return populate_myopic_sparse<DimensionSelectionType::FULL>(row, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_myopic_sparse<DimensionSelectionType::BLOCK>(row, block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return populate_myopic_sparse<DimensionSelectionType::INDEX>(row, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
private:
    template<typename ... Args_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > populate_oracular_dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Args_&& ... args) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::OracularPerpendicularDense<Value_, Index_> >(my_matrix.get(), my_subset, row, std::move(oracle), std::forward<Args_>(args)...); 
        } else {
            return std::make_unique<DelayedSubsetSorted_internal::ParallelDense<true, Value_, Index_> >(my_matrix.get(), my_subset, row, std::move(oracle), std::forward<Args_>(args)...);
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
    template<DimensionSelectionType selection_, typename ... Args_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > populate_oracular_sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Args_&& ... args) const {
        if (row == my_by_row) {
            return std::make_unique<subset_utils::OracularPerpendicularSparse<Value_, Index_> >(my_matrix.get(), my_subset, row, std::move(oracle), std::forward<Args_>(args)...); 
        } else {
            return populate_sparse<selection_, true>(row, std::move(oracle), std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return populate_oracular_sparse<DimensionSelectionType::FULL>(row, std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate_oracular_sparse<DimensionSelectionType::BLOCK>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return populate_oracular_sparse<DimensionSelectionType::INDEX>(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

}

#endif
