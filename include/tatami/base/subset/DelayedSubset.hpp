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
 * @brief Delayed subsetting of a matrix with general indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of arbitrary indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam T Type of matrix value.
 * @tparam V Vector containing the subset indices.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX, class V>
class DelayedSubset : public Matrix<T, IDX> {
private:
    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<V>()[0])>::type>::type V_type;

    static void finish_assembly(
        const std::vector<std::pair<V_type, size_t> >& collected,
        const V& indices, 
        std::vector<size_t>& reverse_mapping,
        std::vector<IDX>& unique_and_sorted,
        size_t mapping_dim,
        std::vector<std::pair<size_t, size_t> >& mapping_duplicates,
        std::vector<IDX>& mapping_duplicates_pool
    ) {
        unique_and_sorted.reserve(indices.size());
        reverse_mapping.resize(indices.size());

        mapping_duplicates.resize(mapping_dim);
        mapping_duplicates_pool.reserve(indices.size());

        for (size_t i = 0, end = collected.size(); i < end; ++i) {
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
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `MARGIN = 0`) or columns (if `MARGIN = 1`).
     * These may be duplicated and/or unsorted.
     */
    DelayedSubset(std::shared_ptr<const Matrix<T, IDX> > p, V idx) : mat(std::move(p)), indices(std::move(idx)) {
        std::vector<std::pair<V_type, size_t> > collected;
        collected.reserve(indices.size());
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            collected.emplace_back(indices[i], i);
        }
        std::sort(collected.begin(), collected.end());

        finish_assembly(
            collected,
            indices, 
            reverse_mapping,
            unique_and_sorted,
            MARGIN == 0 ? mat->nrow() : mat->ncol(),
            mapping_duplicates,
            mapping_duplicates_pool
        );

        return;
    }

    /**
     * @cond
     */
    DelayedSubset(std::shared_ptr<const Matrix<T, IDX> > p, const std::vector<std::pair<V_type, size_t> >& collected, V idx) : mat(std::move(p)), indices(std::move(idx)) {
        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            unique_and_sorted,
            MARGIN == 0 ? mat->nrow() : mat->ncol(),
            mapping_duplicates,
            mapping_duplicates_pool
        );
    }
    /**
     * @endcond
     */

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    V indices;

    std::vector<size_t> reverse_mapping;
    std::vector<IDX> unique_and_sorted;
    std::vector<std::pair<size_t, size_t> > mapping_duplicates; // holds (position, size) in the pool.
    std::vector<IDX> mapping_duplicates_pool; 

public:
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }

    size_t ncol() const {
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

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

private:
    struct DenseBase {
        DenseBase(size_t bufsize, const WorkspaceOptions& opt) : buffer(bufsize) {}
        std::vector<T> buffer;
    };

    struct SparseBase {
        /*
         * The behavior of the buffers depends on whether we want sorted output.
         *
         * If we want sorted output, then the logic roughly follows that of DelayedSubsetUnique::SparseBase.
         * We can directly extract into the user-supplied ibuffer/vbuffer and then copy everything into sortspace for resorting. 
         * The only exception is if the doesn't want the indices, in which case we need to do some allocation.
         *
         * If we don't want sorted output, then the logic follows that of DelayedSubsetSorted::SparseBase.
         * Here, we can't extract into the user-supplied ibuffer/vbuffer, because we need some
         * place to hold the to-be-processed values while we expand the duplicates.
         * Of course, we can omit the vbuffer allocation if the user doesn't want the values.
         */ 

        SparseBase(size_t bufsize, const WorkspaceOptions& opt) : report_index(sparse_extract_index(opt.sparse_extract_mode)), needs_sort(opt.sparse_ordered_index) {
            if (needs_sort) {
                if (!sparse_extract_index(opt.sparse_extract_mode)) {
                    ibuffer.resize(bufsize);
                }
                sortspace.reserve(bufsize);
            } else {
                if (sparse_extract_value(opt.sparse_extract_mode)) {
                    vbuffer.resize(bufsize);
                }
                ibuffer.resize(bufsize);
            }
        }

        bool report_index;
        bool needs_sort;
        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;
        std::vector<std::pair<IDX, T> > sortspace;
    };

    template<bool SPARSE>
    using ConditionalBase = typename std::conditional<SPARSE, SparseBase, DenseBase>::type;

    template<bool WORKROW, bool SPARSE>
    auto define_internal_workspace(std::vector<IDX> subset, WorkspaceOptions opt) const {
        if (opt.sparse_extract_mode == SparseExtractMode::VALUE) {
            // Making sure we extract the indices to do the expansion.
            opt.sparse_extract_mode = SparseExtractMode::BOTH;
        } else if (opt.sparse_extract_mode == SparseExtractMode::NONE) {
            // Still need indices to get the duplicate counts.
            opt.sparse_extract_mode = SparseExtractMode::INDEX;
        }

        // Turning off the sorting to enable possible optimizations in the underlying matrix.
        // We don't need sorted output as we'll be resorting ourselves later.
        opt.sparse_ordered_index = false;
        return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), opt);
    }

    template<class InputWorkspace>
    SparseRange<T, IDX> reorganize_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work, const std::vector<std::pair<size_t, size_t> >& dups, const std::vector<IDX>& pool) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if (!work->needs_sort) {
            // Just calling the expansion method used in DelayedSubsetSorted,
            // and not bothering to copy everything into sortspace and back
            // when we don't need to sort it.
            return subset_utils::remap_sparse_duplicates<WORKROW>(mat.get(), i, vbuffer, ibuffer, work, dups, pool);
        }

        IDX* iin = (work->ibuffer.empty() ? ibuffer : work->ibuffer.data());
        auto raw = extract_sparse<WORKROW>(mat.get(), i, vbuffer, iin, work->internal.get());

        auto& sortspace = work->sortspace;
        sortspace.clear();

        for (size_t i = 0; i < raw.number; ++i) {
            const auto& pool_pos = dups[raw.index[i]];
            size_t pool_end = pool_pos.first + pool_pos.second;

            if (raw.value) {
                for (size_t j = pool_pos.first; j < pool_end; ++j) {
                    sortspace.emplace_back(pool[j], raw.value[i]);
                }
            } else {
                for (size_t j = pool_pos.first; j < pool_end; ++j) {
                    sortspace.emplace_back(pool[j], 0);
                }
            }
        }

        std::sort(sortspace.begin(), sortspace.end());

        if (work->report_index) {
            auto icopy = ibuffer;
            for (const auto& x : sortspace) {
                *icopy= x.first;
                ++icopy;
            }
        } else {
            ibuffer = NULL;
        }

        if (raw.value) {
            auto vcopy = vbuffer;
            for (const auto& x : sortspace) {
                *vcopy = x.second;
                ++vcopy;
            }
        } else {
            vbuffer = NULL;
        }

        return SparseRange<T, IDX>(sortspace.size(), vbuffer, ibuffer);
    }

    /*****************************************
     ************ Full extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct AlongWorkspace : public Workspace<WORKROW, SPARSE>, public ConditionalBase<SPARSE> {
        AlongWorkspace(size_t bufsize, const WorkspaceOptions& opt) : ConditionalBase<SPARSE>(bufsize, opt) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<Workspace<WORKROW, SPARSE> > create_new_workspace(const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), opt);
        } else {
            auto ptr = new AlongWorkspace<WORKROW, SPARSE>(unique_and_sorted.size(), opt); 
            std::shared_ptr<Workspace<WORKROW, SPARSE> > output(ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(unique_and_sorted, opt); // copy is deliberate.
            return output;
        }
    }

    template<bool WORKROW>
    const T* get_dense(size_t i, T* buffer, DenseWorkspace<WORKROW>* work) const {
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), indices[i], buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<WORKROW, false>*>(work);
            return subset_utils::remap_dense<WORKROW>(mat.get(), i, buffer, wptr, reverse_mapping);
        }
    }

    template<bool WORKROW>
    SparseRange<T, IDX> get_sparse(size_t i, T* out_values, IDX* out_indices, SparseWorkspace<WORKROW>* work) const {
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), indices[i], out_values, out_indices, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<WORKROW, true>*>(work);
            return reorganize_sparse(i, out_values, out_indices, wptr, mapping_duplicates, mapping_duplicates_pool);
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

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowWorkspace* work) const {
        return get_sparse(r, out_values, out_indices, work);
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnWorkspace* work) const {
        return get_sparse(c, out_values, out_indices, work);
    }

    /*****************************************
     *********** Block extraction ************
     *****************************************/
private:
    struct DenseBase2 : public DenseBase, public subset_utils::DenseSupplement {
        DenseBase2(size_t bufsize, const WorkspaceOptions& opt, subset_utils::DenseSupplement host) : DenseBase(bufsize, opt), subset_utils::DenseSupplement(std::move(host)) {}
    };

    struct SparseBase2 : public SparseBase, public subset_utils::SparseSupplement<IDX> {
        SparseBase2(size_t bufsize, const WorkspaceOptions& opt, subset_utils::SparseSupplement<IDX> host) : SparseBase(bufsize, opt), subset_utils::SparseSupplement<IDX>(std::move(host)) {}
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
    SparseRange<T, IDX> get_sparse2(size_t i, T* out_values, IDX* out_indices, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), indices[i], out_values, out_indices, work);
        } else {
            auto wptr = static_cast<TargetWorkspace<WORKROW, true>*>(work);
            return reorganize_sparse(i, out_values, out_indices, wptr, wptr->mapping_duplicates, wptr->mapping_duplicates_pool);
        }
    }

    template<class InputWorkspace, class Function>
    std::vector<IDX> transplant_indices(std::vector<std::pair<V_type, size_t> >& collected, InputWorkspace* ptr, Function to_index) const {
        std::sort(collected.begin(), collected.end());

        constexpr bool SPARSE = InputWorkspace::sparse;
        if constexpr(SPARSE) {
            size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
            ptr->mapping_duplicates.resize(mapping_dim);
            ptr->mapping_duplicates_pool.reserve(collected.size());
        } else {
            ptr->reverse_mapping.resize(collected.size());
        }

        std::vector<IDX> local;
        local.reserve(collected.size());

        for (size_t i = 0, end = collected.size(); i < end; ++i) {
            const auto& current = collected[i];
            bool diff = (local.empty() || current.first != local.back());
            if (diff) {
                local.push_back(current.first);
            }

            if constexpr(SPARSE) {
                auto& range = ptr->mapping_duplicates[current.first];
                if (diff) {
                    range.first = ptr->mapping_duplicates_pool.size();
                }
                ++range.second;
                ptr->mapping_duplicates_pool.push_back(to_index(current.second));
            } else {
                ptr->reverse_mapping[current.second] = local.size() - 1;
            }
        }

        return local;
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongBlockWorkspace(size_t s, size_t l, size_t bufsize, const WorkspaceOptions& opt, subset_utils::ConditionalSupplement<IDX, SPARSE> host) : 
            BlockWorkspace<WORKROW, SPARSE>(s, l), ConditionalBase2<SPARSE>(bufsize, opt, std::move(host)) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), start, length, opt);

        } else {
            size_t end = start + length;
            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(length);
            for (size_t i = start; i < end; ++i) {
                collected.emplace_back(indices[i], i - start);
            }

            // Creating the temporary because we need to know the number of
            // unique indices in 'local' before constructing our actual workspace.
            subset_utils::ConditionalSupplement<IDX, SPARSE> temp;
            auto local = transplant_indices(collected, &temp, [&](size_t i) -> size_t { return i + start; });

            auto ptr = new AlongBlockWorkspace<WORKROW, SPARSE>(start, length, local.size(), opt, std::move(temp));
            std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > output(ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(std::move(local), opt);

            return output;
        }
    }

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(s, l, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(s, l, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        return get_dense2<AlongBlockWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return get_dense2<AlongBlockWorkspace>(c, buffer, work);
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(s, l, opt);
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t s, size_t l, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(s, l, opt);
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowBlockWorkspace* work) const {
        return get_sparse2<AlongBlockWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnBlockWorkspace* work) const {
        return get_sparse2<AlongBlockWorkspace>(c, out_values, out_indices, work);
    }

    /*****************************************
     *********** Index extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongIndexWorkspace(std::vector<IDX> i, size_t bufsize, const WorkspaceOptions& opt, subset_utils::ConditionalSupplement<IDX, SPARSE> host) : 
            IndexWorkspace<IDX, WORKROW, SPARSE>(i.size()), ConditionalBase2<SPARSE>(bufsize, opt, std::move(host)), indices_(std::move(i)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > create_new_workspace(std::vector<IDX> subset_, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset_), opt);

        } else {
            size_t length = subset_.size();
            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(length);
            for (size_t i = 0; i < length; ++i) {
                collected.emplace_back(indices[subset_[i]], i);
            }

            // Creating the temporary because we need to know the number of
            // unique indices in 'local' before constructing our actual workspace.
            subset_utils::ConditionalSupplement<IDX, SPARSE> temp;
            auto local = transplant_indices(collected, &temp, [&](size_t i) -> size_t { return subset_[i]; });

            auto ptr = new AlongIndexWorkspace<WORKROW, SPARSE>(std::move(subset_), local.size(), opt, std::move(temp));
            std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > output(ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(std::move(local), opt);

            return output;
        }
    }


public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, false>(std::move(subset), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, false>(std::move(subset), opt);
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        return get_dense2<AlongIndexWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        return get_dense2<AlongIndexWorkspace>(c, buffer, work);
    }

    std::shared_ptr<SparseRowIndexWorkspace<IDX>> sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(std::move(subset), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX>> sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(std::move(subset), opt);
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowIndexWorkspace<IDX>* work) const {
        return get_sparse2<AlongIndexWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnIndexWorkspace<IDX>* work) const {
        return get_sparse2<AlongIndexWorkspace>(c, out_values, out_indices, work);
    }
};

}

#endif
