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

template<typename Index_, class SubsetStorage_, class ToIndex_>
DenseParallelResults<Index_> format_dense_parallel(const SubsetStorage_& subset, const Index_ len, const ToIndex_ to_index) {
    std::vector<std::pair<Index_, Index_> > collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(subset[to_index(i)], i);
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
        auto processed = format_dense_parallel<Index_>(subset, subset.size(), [&](Index_ i) -> Index_ { return i; });
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
        auto processed = format_dense_parallel<Index_>(subset, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
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
        auto processed = format_dense_parallel<Index_>(subset, indices.size(), [&](Index_ i) -> Index_ { return indices[i]; });
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
        resize_container_to_Index_size(my_holding_vbuffer, processed.sorted.size()); // processed.sorted.size() should fit in an Index_, hence the cast is safe.
        my_ext = new_extractor<false, oracle_>(matrix, row, std::move(oracle), std::move(processed.sorted), opt);
        my_permutation = std::move(processed.permutation);
    }

public:
    const Value_* fetch(const Index_ i, Value_* const buffer) {
        auto src = my_ext->fetch(i, my_holding_vbuffer.data());

        // 'input' and 'output' should not point to the same array. In theory, it
        // is possible to do an in-place permutation, but this requires another
        // array anyway to track the permutation status, so we'll just keep it simple.
        for (const auto p : my_permutation) {
            buffer[p] = *src;
            ++src;
        }

        return buffer;
    }

private:
    std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > my_ext;
    std::vector<Value_> my_holding_vbuffer;
    std::vector<Index_> my_permutation;
};

template<typename Index_, class SubsetStorage_, class ToIndex_>
std::vector<Index_> format_sparse_parallel(const SubsetStorage_& subset, const Index_ len, const ToIndex_ to_index) {
    std::vector<Index_> collected;
    collected.reserve(len);
    for (Index_ i = 0; i < len; ++i) {
        collected.emplace_back(subset[to_index(i)]);
    }
    std::sort(collected.begin(), collected.end());
    return collected;
}

template<bool oracle_, typename Value_, typename Index_>
class ParallelSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    template<class SubsetStorage_>
    ParallelSparse(
        const Matrix<Value_, Index_>& matrix,
        const SubsetStorage_& subset,
        const std::vector<Index_>& remap,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) : 
        my_remapping(remap) 
    {
        auto processed = format_sparse_parallel<Index_>(subset, subset.size(), [](Index_ i) -> Index_ { return i; });
        initialize(matrix, std::move(processed), row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelSparse(
        const Matrix<Value_, Index_>& matrix,
        const SubsetStorage_& subset,
        const std::vector<Index_>& remap,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) : 
        my_remapping(remap) 
    {
        auto processed = format_sparse_parallel<Index_>(subset, block_length, [&](Index_ i) -> Index_ { return i + block_start; });
        initialize(matrix, std::move(processed), row, std::move(oracle), opt);
    }

    template<class SubsetStorage_>
    ParallelSparse(
        const Matrix<Value_, Index_>& matrix,
        const SubsetStorage_& subset,
        const std::vector<Index_>& remap,
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) : 
        my_remapping(remap)
    {
        const auto& indices = *indices_ptr;
        auto processed = format_sparse_parallel<Index_>(subset, indices.size(), [&](Index_ i) -> Index_ { return indices[i]; });
        initialize(matrix, std::move(processed), row, std::move(oracle), opt);
    }

private:
    void initialize(const Matrix<Value_, Index_>& matrix, std::vector<Index_> sorted, const bool row, MaybeOracle<oracle_, Index_> oracle, Options opt) {
        my_needs_value = opt.sparse_extract_value;
        my_needs_index = opt.sparse_extract_index;
        my_needs_sort = opt.sparse_ordered_index;

        // The conditionals here mirror those in 'fetch',
        // to self-document the case where each of the temporaries are needed.
        if (!my_needs_sort) {
            if (my_needs_index) {
                ; // no 'my_holding_ibuffer' required as a user-provided 'index_buffer' should be available.
            }

        } else if (my_needs_value) {
            opt.sparse_extract_index = true;
            my_sortspace.reserve(sorted.size());
            if (my_needs_index) {
                // no 'my_holding_ibuffer' required as a user-provided 'index_buffer' should be available.
            } else {
                // Needs 'my_holding_ibuffer' as user-provided 'index_buffer' may be NULL.
                resize_container_to_Index_size(my_holding_ibuffer, sorted.size()); // sorted.size() should fit in an Index_, hence the cast is safe.
            }

        } else if (my_needs_index) {
            ; // no 'my_holding_ibuffer' required as a user-provided 'index_buffer' should be available.
        }

        my_ext = new_extractor<true, oracle_>(matrix, row, std::move(oracle), std::move(sorted), opt);
    }

public:
    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        auto input = my_ext->fetch(i, value_buffer, (my_holding_ibuffer.empty() ? index_buffer : my_holding_ibuffer.data()));

        // Pointers in 'input' and the 'buffer' pointers may point to the same array,
        // as we're either just modifiying in place or we're copying to 'my_sortspace'.
        if (!my_needs_sort) {
            if (my_needs_index) {
                for (Index_ i = 0; i < input.number; ++i) {
                    index_buffer[i] = my_remapping[input.index[i]];
                }
                input.index = index_buffer;
            }

        } else if (my_needs_value) {
            // We assume that the indices have already been extracted for sorting
            // purposes, even if they weren't actually requested.
            my_sortspace.clear();
            for (Index_ i = 0; i < input.number; ++i) {
                my_sortspace.emplace_back(my_remapping[input.index[i]], input.value[i]);
            }
            std::sort(my_sortspace.begin(), my_sortspace.end());

            auto vcopy = value_buffer;
            for (const auto& ss : my_sortspace) {
                *vcopy = ss.second;
                ++vcopy;
            }
            input.value = value_buffer;

            if (my_needs_index) {
                auto icopy = index_buffer;
                for (const auto& ss : my_sortspace) {
                    *icopy = ss.first;
                    ++icopy;
                }
                input.index = index_buffer;
            } else {
                input.index = NULL;
            }

        } else if (my_needs_index) {
            for (Index_ i = 0; i < input.number; ++i) {
                index_buffer[i] = my_remapping[input.index[i]];
            }
            std::sort(index_buffer, index_buffer + input.number);
            input.index = index_buffer;
        }

        return input;
    }

private:
    const std::vector<Index_>& my_remapping;
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;
    bool my_needs_value, my_needs_index, my_needs_sort;
    std::vector<std::pair<Index_, Value_> > my_sortspace;
    std::vector<Index_> my_holding_ibuffer;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed subsetting of a matrix with unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of unique indices.
 * This operation is "delayed" in that it is only evaluated when data is extracted from the matrix.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type of index value.
 * @tparam SubsetStorage_ Vector containing the subset indices.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 */
template<typename Value_, typename Index_, class SubsetStorage_>
class DelayedSubsetUnique final : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying (pre-subset) matrix.
     * @param subset Vector of 0-based indices to use for subsetting on the rows (if `by_row = true`) or columns (otherwise).
     * This should be unique, but may be unsorted.
     * @param by_row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     * @param check Whether to check `subset` for unique values.
     */
    DelayedSubsetUnique(    
        std::shared_ptr<const Matrix<Value_, Index_> > matrix,
        SubsetStorage_ subset,
        const bool by_row,
        const bool check = true
    ) : 
        my_matrix(std::move(matrix)),
        my_subset(std::move(subset)),
        my_by_row(by_row)
    {
        const Index_ fulldim = my_by_row ? my_matrix->nrow() : my_matrix->ncol();
        const auto nsub = my_subset.size();

        if (check) {
            auto checks = create_container_of_Index_size<std::vector<unsigned char> >(fulldim);
            for (I<decltype(nsub)> i = 0; i < nsub; ++i) {
                auto& found = checks[my_subset[i]];
                if (found) {
                    throw std::runtime_error("my_subset should be unique");
                } 
                found = 1;
            }
        }

        resize_container_to_Index_size(my_mapping_single, fulldim);
        for (I<decltype(nsub)> i = 0; i < nsub; ++i) {
            my_mapping_single[my_subset[i]] = i;
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > my_matrix;
    SubsetStorage_ my_subset;
    bool my_by_row;
    std::vector<Index_> my_mapping_single;

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

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::sparse_column;

    using Matrix<Value_, Index_>::sparse_row;

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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelDense<false, Value_, Index_> >(
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
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return populate_myopic_dense(row, std::move(indices_ptr), opt);
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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelSparse<false, Value_, Index_> >(
                *my_matrix,
                my_subset,
                my_mapping_single,
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
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return populate_myopic_sparse(row, std::move(indices_ptr), opt);
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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelDense<true, Value_, Index_> >(
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
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return populate_oracular_dense(row, std::move(oracle), std::move(indices_ptr), opt);
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
            return std::make_unique<DelayedSubsetUnique_internal::ParallelSparse<true, Value_, Index_> >(
                *my_matrix,
                my_subset,
                my_mapping_single,
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
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        return populate_oracular_sparse(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

}

#endif
