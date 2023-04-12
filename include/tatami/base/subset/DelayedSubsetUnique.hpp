#ifndef TATAMI_DELAYED_SUBSET_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_UNIQUE_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetUnique.hpp
 *
 * @brief Delayed subsetting by unique row/column indices.
 */

namespace tatami {

/**
 * @brief Delayed subsetting of a matrix with unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of unique indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam T Type of matrix value.
 * @tparam V Vector containing the subset indices.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX, class V>
class DelayedSubsetUnique : public Matrix<T, IDX> {
private:
    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<V>()[0])>::type>::type V_type;

    static void finish_assembly(
        const std::vector<std::pair<V_type, size_t> >& collected,
        const V& indices, 
        std::vector<size_t>& reverse_mapping,
        std::vector<IDX>& unique_and_sorted,
        size_t mapping_dim,
        std::vector<IDX>& mapping_single
    ) {
        unique_and_sorted.reserve(indices.size());
        reverse_mapping.resize(indices.size());

        for (size_t i = 0, end = collected.size(); i < end; ++i) {
            const auto& current = collected[i];
            unique_and_sorted.push_back(current.first);
            reverse_mapping[current.second] = unique_and_sorted.size() - 1;
        }

        mapping_single.resize(mapping_dim);
        for (IDX  i = 0, end = indices.size(); i < end; ++i) {
            mapping_single[indices[i]] = i;
        }
    }

public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `MARGIN = 0`) or columns (if `MARGIN = 1`).
     * This should be unique, but may be unsorted.
     * @param check Whether to check `idx` for unique values.
     */
    DelayedSubsetUnique(std::shared_ptr<const Matrix<T, IDX> > p, V idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        std::vector<std::pair<V_type, size_t> > collected;
        collected.reserve(indices.size());
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            collected.emplace_back(indices[i], i);
        }
        std::sort(collected.begin(), collected.end());

        if (check) {
            for (size_t i = 0, end = collected.size(); i < end; ++i) {
                if (i) {
                    if (collected[i].first == collected[i-1].first) {
                        throw std::runtime_error("indices should be unique");
                        break;
                    }
                }
            }
        }

        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            unique_and_sorted,
            MARGIN == 0 ? mat->nrow() : mat->ncol(),
            mapping_single
        );
    }

    /**
     * @cond
     */
    DelayedSubsetUnique(std::shared_ptr<const Matrix<T, IDX> > p, const std::vector<std::pair<V_type, size_t> >& collected, V idx) : mat(std::move(p)), indices(std::move(idx)) {
        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            unique_and_sorted,
            MARGIN == 0 ? mat->nrow() : mat->ncol(),
            mapping_single
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
    std::vector<IDX> mapping_single;

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
         * There's no need to have a vbuffer here, because we can directly extract into the user-supplied vbuffer if we want the values. 
         * The inner extraction is guaranteed to be no greater than the vbuffer length, as it's a 1:1 mapping between unique indices in the row/column call and the unique internal subset indices.
         * Temporaries are not a concern because we end up copying everything into sortspace for re-sorting anyway.
         *
         * Similarly, an ibuffer is _usually_ unnecessary as we can write directly into the user-supplied ibuffer.
         * There's only a need to allocate our own ibuffer if (i) the user wants the values but not the indices, and (ii) they want the output to be sorted. 
         * In this case, the user-supplied ibuffer (i) may not be valid but (ii) we still need the indices, so we need to create our ibuffer.
         */ 

        SparseBase(size_t bufsize, const WorkspaceOptions& opt) : 
            report_index(sparse_extract_index(opt.mode)),
            needs_sort(opt.sorted),
            ibuffer(opt.mode == SparseExtractMode::VALUE && needs_sort ? bufsize : 0) 
        {
            if (needs_sort) {
                sortspace.reserve(bufsize);
            }
        }

        bool report_index;
        bool needs_sort;
        std::vector<IDX> ibuffer;
        std::vector<std::pair<IDX, T> > sortspace;
    };

    template<bool SPARSE>
    using ConditionalBase = typename std::conditional<SPARSE, SparseBase, DenseBase>::type;

    template<bool WORKROW, bool SPARSE>
    auto define_internal_workspace(std::vector<IDX> subset, WorkspaceOptions opt) const {
        if (opt.mode == SparseExtractMode::VALUE && opt.sorted) {
            // Making sure we extract the indices if we want the sorted values.
            opt.mode = SparseExtractMode::BOTH;
        }

        // Turning off the sorting to enable possible optimizations in the underlying matrix.
        // We don't need sorted output as we'll be resorting ourselves later.
        opt.sorted = false;
        return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), opt);
    }

    template<bool WORKROW, class InputWorkspace>
    SparseRange<T, IDX> reorganize_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work) const {
        if (!(work->needs_sort)) {
            // When we don't need sorting, validity of ibuffer and vbuffer/ should be consistent 
            // with the extraction mode used to construct work->internal.
            auto raw = extract_sparse<WORKROW>(mat.get(), i, vbuffer, ibuffer, work->internal.get());
            if (raw.index) {
                auto icopy = ibuffer;
                for (size_t i = 0; i < raw.number; ++i, ++icopy) {
                    *icopy = mapping_single[raw.index[i]];
                }
            }
            return raw;
        }

        // If the workspace's ibuffer is empty, we're either extracting the indices
        // directly into the user's ibuffer, or we don't need the indices at all.
        // Either way, it doesn't hurt to use the user's ibuffer.
        IDX* iin = (work->ibuffer.empty() ? ibuffer : work->ibuffer.data()); 

        auto raw = extract_sparse<WORKROW>(mat.get(), i, vbuffer, iin, work->internal.get());
        if (!raw.value && !raw.index) {
            return raw;
        }

        auto& sortspace = work->sortspace;
        sortspace.clear();

        // raw.index is guaranteed to be valid as this point, as the internal workspace
        // can never be set to VALUE-only if it's requesting sorted output.
        if (raw.value) {
            for (size_t i = 0; i < raw.number; ++i) {
                sortspace.emplace_back(mapping_single[raw.index[i]], raw.value[i]);
            }
        } else {
            for (size_t i = 0; i < raw.number; ++i) {
                sortspace.emplace_back(mapping_single[raw.index[i]], 0);
            }
        }

        std::sort(sortspace.begin(), sortspace.end());

        if (raw.value) {
            auto vcopy = vbuffer;
            for (const auto& x : sortspace) {
                *vcopy = x.second;
                ++vcopy;
            }
        } else {
            vbuffer = NULL;
        }

        if (work->report_index) {
            auto icopy = ibuffer;
            for (const auto& x : sortspace) {
                *icopy = x.first;
                ++icopy;
            }
        } else {
            ibuffer = NULL;
        }

        return SparseRange<T, IDX>(raw.number, vbuffer, ibuffer);
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
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(unique_and_sorted, opt); // deliberate copy here.
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

    template<template<bool, bool> class TargetWorkspace, class InputWorkspace>
    SparseRange<T, IDX> get_sparse(size_t i, T* out_values, IDX* out_indices, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), indices[i], out_values, out_indices, work);
        } else {
            auto wptr = static_cast<TargetWorkspace<WORKROW, true>*>(work);
            return reorganize_sparse<WORKROW>(i, out_values, out_indices, wptr);
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
        return get_sparse<AlongWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnWorkspace* work) const {
        return get_sparse<AlongWorkspace>(c, out_values, out_indices, work);
    }

    /*****************************************
     *********** Block extraction ************
     *****************************************/
private:
    template<class InputWorkspace>
    std::vector<IDX> transplant_indices(std::vector<std::pair<V_type, size_t> >& collected, InputWorkspace* work) const {
        std::sort(collected.begin(), collected.end());

        std::vector<IDX> local;
        local.reserve(collected.size());
        if constexpr(!InputWorkspace::sparse) {
            work->reverse_mapping.resize(collected.size());
        }

        for (size_t i = 0, end = collected.size(); i < end; ++i) {
            local.push_back(collected[i].first);
            if constexpr(!InputWorkspace::sparse) {
                work->reverse_mapping[collected[i].second] = local.size() - 1;
            }
        }

        return local;
    }

    struct DenseBase2 : public DenseBase {
        DenseBase2(size_t bufsize, const WorkspaceOptions& opt) : DenseBase(bufsize, opt) {}
        std::vector<size_t> reverse_mapping;
    };

    template<bool SPARSE>
    using ConditionalBase2 = typename std::conditional<SPARSE, SparseBase, DenseBase2>::type;

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongBlockWorkspace(size_t s, size_t l, const WorkspaceOptions& opt) : BlockWorkspace<WORKROW, SPARSE>(s, l), ConditionalBase2<SPARSE>(l, opt) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), start, length, opt);
        } else {
            auto ptr = new AlongBlockWorkspace<WORKROW, SPARSE>(start, length, opt);
            std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > output(ptr);

            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(length);
            size_t end = start + length;
            for (size_t i = start; i < end; ++i) {
                collected.emplace_back(indices[i], i - start);
            }

            auto local = transplant_indices(collected, ptr);
            ptr->internal = define_internal_workspace<WORKROW, SPARSE>(std::move(local), opt); 
            return output;
        }
    }

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

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowBlockWorkspace* work) const {
        return get_sparse<AlongBlockWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnBlockWorkspace* work) const {
        return get_sparse<AlongBlockWorkspace>(c, out_values, out_indices, work);
    }

    /*****************************************
     *********** Index extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongIndexWorkspace(std::vector<IDX> i, const WorkspaceOptions& opt) : 
            IndexWorkspace<IDX, WORKROW, SPARSE>(i.size()), ConditionalBase2<SPARSE>(i.size(), opt), indices_(std::move(i)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > create_new_workspace(std::vector<IDX> subset_, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset_), opt);
        } else {
            auto ptr = new AlongIndexWorkspace<WORKROW, SPARSE>(std::move(subset_), opt);
            std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > output(ptr);

            const auto& subset = ptr->indices_;
            size_t length = subset.size();

            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(length);
            for (size_t i = 0; i < length; ++i) {
                collected.emplace_back(indices[subset[i]], i);
            }

            auto local = transplant_indices(collected, ptr);
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

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(std::move(subset), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, true>(std::move(subset), opt);
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowIndexWorkspace<IDX>* work) const {
        return get_sparse<AlongIndexWorkspace>(r, out_values, out_indices, work);
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnIndexWorkspace<IDX>* work) const {
        return get_sparse<AlongIndexWorkspace>(c, out_values, out_indices, work);
    }
};

}

#endif
