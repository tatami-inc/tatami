#ifndef TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP

#include <algorithm>
#include <memory>

#include "../base/Matrix.hpp"
#include "../utils/copy.hpp"
#include "utils.hpp"

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

template<typename Index_, class SubsetStorage_>
VectorPtr<Index_> create(const SubsetStorage_& subset) {
    return std::make_shared<std::vector<Index_> >(subset.begin(), subset.end());
}

template<typename Index_, class SubsetStorage_>
VectorPtr<Index_> create(const SubsetStorage_& subset, const Index_ block_start, const Index_ block_length) {
    auto pistart = subset.begin() + block_start;
    return std::make_shared<std::vector<Index_> >(pistart, pistart + block_length);
}

template<typename Index_, class SubsetStorage_>
VectorPtr<Index_> create(const SubsetStorage_& subset, const VectorPtr<Index_>& indices_ptr) {
    auto rawptr = new std::vector<Index_>;
    VectorPtr<Index_> outptr(rawptr);
    auto& output = *rawptr;

    const auto& input = *indices_ptr;
    output.reserve(input.size());
    for (auto i : input) {
        output.push_back(subset[i]);
    }

    return outptr;
}

template<bool oracle_, typename Value_, typename Index_, class SubsetStorage_>
std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > create_parallel_dense(
    const Matrix<Value_, Index_>& matrix,
    const SubsetStorage_& subset, 
    const bool row, 
    MaybeOracle<oracle_, Index_> oracle, 
    const Options& opt
) {
    return new_extractor<false, oracle_>(matrix, row, std::move(oracle), create<Index_>(subset), opt);
}

template<bool oracle_, typename Value_, typename Index_, class SubsetStorage_>
std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > create_parallel_dense(
    const Matrix<Value_, Index_>& matrix,
    const SubsetStorage_& subset,
    const bool row,
    MaybeOracle<oracle_, Index_> oracle,
    const Index_ block_start,
    const Index_ block_length,
    const Options& opt)
{
    return new_extractor<false, oracle_>(matrix, row, std::move(oracle), create<Index_>(subset, block_start, block_length), opt);
}

template<bool oracle_, typename Value_, typename Index_, class SubsetStorage_>
std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > create_parallel_dense(
    const Matrix<Value_, Index_>& matrix,
    const SubsetStorage_& subset,
    const bool row,
    MaybeOracle<oracle_, Index_> oracle,
    VectorPtr<Index_> indices_ptr,
    const Options& opt) 
{
    return new_extractor<false, oracle_>(matrix, row, std::move(oracle), create<Index_>(subset, indices_ptr), opt);
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
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), create<Index_>(subset), opt)),
        my_remapping(remap)
    {}

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
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), create<Index_>(subset, block_start, block_length), opt)),
        my_remapping(remap)
    {}

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
        my_ext(new_extractor<true, oracle_>(matrix, row, std::move(oracle), create<Index_>(subset, indices_ptr), opt)),
        my_remapping(remap)
    {}

public:
    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        auto out = my_ext->fetch(i, value_buffer, index_buffer);
        if (out.index) {
            for (Index_ i = 0; i < out.number; ++i) {
                index_buffer[i] = my_remapping[out.index[i]];
            }
            out.index = index_buffer;
        }
        return out;
    }

private:
    std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > my_ext;
    const std::vector<Index_>& my_remapping;
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
 * @tparam SubsetStorage_ Vector containing the subset indices.
 */
template<typename Value_, typename Index_, class SubsetStorage_>
class DelayedSubsetSortedUnique final : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrix Pointer to the underlying (pre-subset) matrix.
     * @param subset Vector of 0-based indices to use for subsetting on the rows (if `by_row = true`) or columns (otherwise).
     * This should be sorted and unique.
     * @param by_row Whether to apply the subset to the rows.
     * If false, the subset is applied to the columns.
     * @param check Whether to check `subset` for sorted and unique values.
     */
    DelayedSubsetSortedUnique(
        std::shared_ptr<const Matrix<Value_, Index_> > matrix,
        SubsetStorage_ subset,
        const bool by_row,
        const bool check = true
    ) :
        my_matrix(std::move(matrix)),
        my_subset(std::move(subset)),
        my_by_row(by_row)
    {
        const auto nsub = my_subset.size();
        if (check) {
            for (I<decltype(nsub)> i = 1; i < nsub; ++i) {
                if (my_subset[i] <= my_subset[i-1]) {
                    throw std::runtime_error("subset should be unique and sorted");
                }
            }
        }

        const Index_ mapping_dim = my_by_row ? my_matrix->nrow() : my_matrix->ncol();
        resize_container_to_Index_size(my_mapping_single, mapping_dim);
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
            return DelayedSubsetSortedUnique_internal::create_parallel_dense<false>(
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
            return std::make_unique<DelayedSubsetSortedUnique_internal::ParallelSparse<false, Value_, Index_> >(
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
            return DelayedSubsetSortedUnique_internal::create_parallel_dense<true>(
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
            return std::make_unique<DelayedSubsetSortedUnique_internal::ParallelSparse<true, Value_, Index_> >(
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
