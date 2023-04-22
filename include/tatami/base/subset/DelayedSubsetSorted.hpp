#ifndef TATAMI_DELAYED_SUBSET_SORTED_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_HPP

#include "utils.hpp"
#include <algorithm>
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
 */
template<int MARGIN, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetSorted : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `MARGIN = 0`) or columns (if `MARGIN = 1`).
     * This should be sorted, but may be duplicated.
     * @param check Whether to check `idx` for sorted values.
     */
    DelayedSubsetSorted(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        if (check) {
            for (size_t i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] < indices[i-1]) {
                    throw std::runtime_error("indices should be sorted");
                }
            }
        }

        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
        unique_and_sorted.reserve(indices.size());
        reverse_mapping.reserve(indices.size());
        duplicate_starts.resize(mapping_dim);
        duplicate_lengths.resize(mapping_dim);

        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            Index_ curdex = indices[i];
            if (unique_and_sorted.empty() || curdex != unique_and_sorted.back()) {
                unique_and_sorted.push_back(curdex);
                duplicates_start[curdex] = mapping_duplicates_pool.size();
            }
            reverse_mapping.push_back(unique_and_sorted.size() - 1);
            ++(duplicates_length[curdex]);
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

    std::vector<Index_> unique_and_sorted;
    std::vector<size_t> reverse_mapping;
    std::vector<size_t> duplicate_starts; // holds the start position and the number of duplicates.
    std::vector<size_t> duplicate_lengths; // holds the start position and the number of duplicates.

public:
    Index_ nrow() const {
        if constexpr(MARGIN==0) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }

    Index_ ncol() const {
        if constexpr(MARGIN==0) {
            return mat->ncol();
        } else {
            return indices.size();
        }
    }

    bool sparse() const {
        return mat->sparse();
    }

    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    std::pair<double, double> dimension_preference() const {
        return mat->dimension_preference();
    }

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::sparse_column;

    using Matrix<Value_, Index_>::sparse_row;

private:
    template<class Extractor_>
    static SparseRange<Value_, Index_> remap_sparse_duplicates(
        Extractor_* ex,
        Index_ i, 
        Value_* vbuffer, 
        Index_* ibuffer, 
        const std::vector<std::pair<Index_, size_t> >& dups,
        ) 
    {
        // Allocation status of work->vbuffer depends on the extraction mode used to construct work->internal.
        Value_* vin = work->vbuffer.data();

        // work->ibuffer should always be allocated, as we need this to get the expanded counts.
        Index_* iin = work->ibuffer.data();
        auto raw = extract_sparse<WORKROW>(mat, i, vin, iin, work->internal.get());

        if (!raw.value) {
            vbuffer = NULL;
        }
        if (!(work->report_index)) {
            ibuffer = NULL;
        }

        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        size_t counter = 0;

        for (size_t i = 0; i < raw.number; ++i) {
            const auto& pool_pos = dups[raw.index[i]];
            counter += pool_pos.second;

            if (vcopy) {
                std::fill(vcopy, vcopy + pool_pos.second, raw.value[i]);
                vcopy += pool_pos.second;
            }

            if (icopy) {
                auto istart = pool.begin() + pool_pos.first;
                std::copy(istart, istart + pool_pos.second, icopy);
                icopy += pool_pos.second;
            }
        }

        return SparseRange<Value_, Index_>(counter, vbuffer, ibuffer);
    }

    /**************************************************
     ************ Full parallel extraction ************
     **************************************************/
private:
    template<bool sparse_>
    struct FullParallelWorkspaceBase : public Extractor<DimensionSelectionType::FULL, sparse_, Value_, Index_> {
        FullParallelWorkspaceBase(const DelayedSubsetSorted* p, const Options<Index_>& opt) : parent(p) {
            this->full_length = parent->indices.size();
            internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), parent->unique_and_sorted.data(), parent->unique_and_sorted.size(), opt);
        }
    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        const DelayedSubsetSorted* parent;
    };

    struct DenseFullParallelWorkspaceBase : public FullParallelWorkspaceBase<false> {
        FullParallelWorkspaceBase(const DelayedSubsetSorted* p, const Options<Index_>& opt) : parent(p) {

        }
    };

    /***************************************************
     ************ Block parallel extraction ************
     ***************************************************/
private:
    template<bool sparse_>
    struct BlockParallelWorkspaceBase : public Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> {
        BlockWorkspaceBase(const DelayedSubsetSorted* p, const Options<Index_>& opt, Index_ bs, Index_ bl) : parent(p) {
            this->block_start = bs;
            this->block_length = bl;

            // Finding the edges of the unique subset.
            size_t from = 0, to = 0;
            if (bl) {
                auto left = parent->indices[bs];
                size_t end = unique_and_sorted.size();
                for (; from < end; ++from) {
                    if (unique_and_sorted[from] > left) {
                        --from;
                        break;
                    }
                }

                auto blast = bs + bl - 1;
                auto right = parent->indices[blast];
                for (to = from; to < end; ++to) {
                    if (unique_and_sorted[to] > right) {
                        break;
                    }
                }

                if constexpr(sparse_) {
                    left_dups_lost = parent->duplicate_lengths[left] - 1;
                    for (Index_ b = 1; b < bl; ++b) {
                        if (parent->indices[bs + b] == left) {
                            --left_dups_lost;
                        } else {
                            break;
                        }
                    }

                    right_dups_lost = parent->duplicate_lengths[right] - 1;
                    for (Index_ b = 1; b < bl; b) {
                        if (parent->indices[blast - b] == right) {
                            --right_dups_lost;
                        } else {
                            break;
                        }
                    }
                }
            }

            internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), parent->indices.data() + from, to - from, opt);
        }
    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        const DelayedSubsetSorted* parent;
        size_t left_dups_lost = 0, right_dups_lost = 0;
    };

    template<bool sparse_>
    struct IndexParallelWorkspaceBase : public Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> {
        ParallelWorkspaceBase(const DelayedSubsetSorted* parent, const Options<Index_>& opt, const Index_* is, size_t il) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                indices.reserve(il);
                for (size_t i = 0; i < il; ++i) {
                    indices.push_back(parent->indices[is[i]]);
                }
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), indices.data(), il, opt);

                this->index_length = il;
                std::copy(is, is + il, indices.begin());
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        typename std::conditional<ConditionalSupplement<sparse_, Index_>, const DelayedSubsetSorted*> supplementary;
        typename std::conditional<DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
    };


    template<bool WORKROW, bool SPARSE>
    struct AlongWorkspace : public Workspace<WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        AlongWorkspace(size_t bufsize, const WorkspaceOptions& opt) : ConditionalBase<SPARSE>(bufsize, opt) {}

        std::shared_ptr<IndexWorkspace<Index_, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<Workspace<WORKROW, SPARSE> > create_new_workspace(const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), opt);
        } else {
            auto ptr = new AlongWorkspace<WORKROW, SPARSE>(unique_and_sorted.size(), opt);
            std::shared_ptr<Workspace<WORKROW, SPARSE> > output(ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(unique_and_sorted, opt); // don't move the unique_and_sorted; it's a deliberate copy.
            return output;
        }
    }

    template<class InputWorkspace>
    const T* get_dense(size_t i, T* buffer, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), indices[i], buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<WORKROW, false>*>(work);
            return subset_utils::remap_dense<WORKROW>(mat.get(), i, buffer, wptr, reverse_mapping);
        }
    }

    template<class InputWorkspace>
    SparseRange<Value_, Index_> get_sparse(size_t i, T* out_values, Index_* out_indices, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), indices[i], out_values, out_indices, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<WORKROW, true>*>(work);
            return subset_utils::remap_sparse_duplicates<WORKROW>(mat.get(), i, out_values, out_indices, wptr, mapping_duplicates, mapping_duplicates_pool);
        }
    }

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(opt);
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        return get_dense(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return get_dense(c, buffer, work);
    }

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(opt);
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(opt);
    }

    SparseRange<Value_, Index_> row(size_t r, T* out_values, Index_* out_indices, SparseRowWorkspace* work) const {
        return get_sparse(r, out_values, out_indices, work);
    }

    SparseRange<Value_, Index_> column(size_t c, T* out_values, Index_* out_indices, SparseColumnWorkspace* work) const {
        return get_sparse(c, out_values, out_indices, work);
    }

    /*****************************************
     *********** Block extraction ************
     *****************************************/
private:
    struct DenseBase2 : public DenseBase, public subset_utils::DenseSupplement {
        DenseBase2(size_t bufsize, const WorkspaceOptions& opt, subset_utils::DenseSupplement host) : DenseBase(bufsize, opt), subset_utils::DenseSupplement(std::move(host)) {}
    };

    struct SparseBase2 : public SparseBase, public subset_utils::SparseSupplement<Index_> {
        SparseBase2(size_t bufsize, const WorkspaceOptions& opt, subset_utils::SparseSupplement<Index_> host) : SparseBase(bufsize, opt), subset_utils::SparseSupplement<Index_>(std::move(host)) {}
    };

    template<bool SPARSE>
    using ConditionalBase2 = typename std::conditional<SPARSE, SparseBase2, DenseBase2>::type;

    template<template<bool, bool> class TargetWorkspace, class InputWorkspace>
    const T* get_dense2(size_t i, T* buffer, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), indices[i], buffer, work);
        } else {
            auto wptr = static_cast<TargetWorkspace<WORKROW, false>*>(work);
            return subset_utils::remap_dense<WORKROW>(mat.get(), i, buffer, wptr, wptr->reverse_mapping);
        }
    }

    template<template<bool, bool> class TargetWorkspace, class InputWorkspace>
    SparseRange<Value_, Index_> get_sparse2(size_t i, T* out_values, Index_* out_indices, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), indices[i], out_values, out_indices, work);
        } else {
            auto wptr = static_cast<TargetWorkspace<WORKROW, true>*>(work);
            return subset_utils::remap_sparse_duplicates<WORKROW>(mat.get(), i, out_values, out_indices, wptr, wptr->mapping_duplicates, wptr->mapping_duplicates_pool);
        }
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongBlockWorkspace(size_t start, size_t length, size_t bufsize, const WorkspaceOptions& opt, subset_utils::ConditionalSupplement<Index_, SPARSE> host) : 
            BlockWorkspace<WORKROW, SPARSE>(start, length), ConditionalBase2<SPARSE>(bufsize, opt, std::move(host)) {}

        std::shared_ptr<IndexWorkspace<Index_, WORKROW, SPARSE> > internal;
    };

    template<bool SPARSE, class Function>
    std::vector<Index_> transplant_indices(subset_utils::ConditionalSupplement<Index_, SPARSE>& temp, size_t length, Function to_index) const {
        if constexpr(SPARSE) {
            size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
            temp.mapping_duplicates.resize(mapping_dim);
            temp.mapping_duplicates_pool.reserve(length);
        } else {
            temp.reverse_mapping.reserve(length);
        }

        std::vector<Index_> local;
        local.reserve(length);

        for (size_t x = 0; x < length; ++x) {
            auto i = to_index(x);
            auto s = indices[i];
            bool diff = local.empty() || s != local.back();
            if (diff) {
                local.push_back(s);
            }

            if constexpr(SPARSE) {
                auto& range = temp.mapping_duplicates[s];
                if (diff) {
                    range.first = temp.mapping_duplicates_pool.size();
                }
                ++range.second;
                temp.mapping_duplicates_pool.push_back(i);
            } else {
                temp.reverse_mapping.push_back(local.size() - 1);
            }
        }

        return local;
    }

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), start, length, opt);

        } else {
            // Need to create a temporary object so that we can get the number of unique indices
            // within this block, for use in the constructor for the _actual_ Workspace. 
            subset_utils::ConditionalSupplement<Index_, SPARSE> temp;
            auto local = transplant_indices<SPARSE>(temp, length, [&](size_t i) -> size_t { return i + start; });

            auto ptr = new AlongBlockWorkspace<WORKROW, SPARSE>(start, length, local.size(), opt, std::move(temp));
            std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > output(ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(std::move(local), opt); 

            return output;
        }
    }

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(start, length, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(start, length, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        return get_dense2<AlongBlockWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return get_dense2<AlongBlockWorkspace>(c, buffer, work);
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(start, length, opt);
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(start, length, opt);
    }

    SparseRange<Value_, Index_> row(size_t r, T* out_values, Index_* out_indices, SparseRowBlockWorkspace* work) const {
        return get_sparse2<AlongBlockWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<Value_, Index_> column(size_t c, T* out_values, Index_* out_indices, SparseColumnBlockWorkspace* work) const {
        return get_sparse2<AlongBlockWorkspace>(c, out_values, out_indices, work);
    }

    /*****************************************
     *********** Index extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct AlongIndexWorkspace : public IndexWorkspace<Index_, WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongIndexWorkspace(std::vector<Index_> i, size_t bufsize, const WorkspaceOptions& opt, subset_utils::ConditionalSupplement<Index_, SPARSE> host) : 
            IndexWorkspace<Index_, WORKROW, SPARSE>(i.size()), ConditionalBase2<SPARSE>(bufsize, opt, std::move(host)), indices_(std::move(i)) {}

        std::vector<Index_> indices_;
        const std::vector<Index_>& indices() const { return indices_; }

        std::shared_ptr<IndexWorkspace<Index_, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<IndexWorkspace<Index_, WORKROW, SPARSE> > create_new_workspace(std::vector<Index_> subset_, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset_), opt);

        } else {
            // Need to create a temporary object so that we can get the number of unique indices
            // within this block, for use in the constructor for the _actual_ Workspace. 
            subset_utils::ConditionalSupplement<Index_, SPARSE> temp;
            size_t length = subset_.size();
            auto local = transplant_indices<SPARSE>(temp, subset_.size(), [&](size_t i) -> Index_ { return subset_[i]; });

            auto ptr = new AlongIndexWorkspace<WORKROW, SPARSE>(std::move(subset_), local.size(), opt, std::move(temp));
            std::shared_ptr<IndexWorkspace<Index_, WORKROW, SPARSE> > output(ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(std::move(local), opt); 
            return output;
        }
    }

public:    
    std::shared_ptr<DenseRowIndexWorkspace<Index_> > dense_row_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(std::move(subset), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<Index_> > dense_column_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(std::move(subset), opt);
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<Index_>* work) const {
        return get_dense2<AlongIndexWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<Index_>* work) const {
        return get_dense2<AlongIndexWorkspace>(c, buffer, work);
    }

    std::shared_ptr<SparseRowIndexWorkspace<Index_> > sparse_row_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(std::move(subset), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<Index_> > sparse_column_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(std::move(subset), opt);
    }

    SparseRange<Value_, Index_> row(size_t r, T* out_values, Index_* out_indices, SparseRowIndexWorkspace<Index_>* work) const {
        return get_sparse2<AlongIndexWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<Value_, Index_> column(size_t c, T* out_values, Index_* out_indices, SparseColumnIndexWorkspace<Index_>* work) const {
        return get_sparse2<AlongIndexWorkspace>(c, out_values, out_indices, work);
    }
};

}

#endif
