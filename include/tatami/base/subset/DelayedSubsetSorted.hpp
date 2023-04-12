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
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam T Type of matrix value.
 * @tparam V Vector containing the subset indices.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX, class V>
class DelayedSubsetSorted : public Matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `MARGIN = 0`) or columns (if `MARGIN = 1`).
     * This should be sorted, but may be duplicated.
     * @param check Whether to check `idx` for sorted values.
     */
    DelayedSubsetSorted(std::shared_ptr<const Matrix<T, IDX> > p, V idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
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
        mapping_duplicates.resize(mapping_dim);
        mapping_duplicates_pool.reserve(indices.size());

        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            IDX curdex = indices[i];
            auto& range = mapping_duplicates[curdex];
            if (unique_and_sorted.empty() || curdex != unique_and_sorted.back()) {
                unique_and_sorted.push_back(curdex);
                range.first = mapping_duplicates_pool.size();
            }
            mapping_duplicates_pool.push_back(i);
            reverse_mapping.push_back(unique_and_sorted.size() - 1);
            ++range.second;
        }
    }

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    V indices;

    std::vector<IDX> unique_and_sorted;
    std::vector<size_t> reverse_mapping;
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
         * In general, we can't extract indices/values directly into the user-supplied ibuffer/vbuffer, 
         * because we need some place to hold the to-be-processed values while we expand the duplicates. 
         * Thus, we need to allocate the workspace's ibuffer/vbuffer to the length of the unique indices. 
         *
         * Obviously, if the values are not required, we can skip the allocation of vbuffer.
         * However, ibuffer is still required, even if the indices are not required, to figure out the expansion of each value.
         * In fact, even if both are not required, ibuffer is still required because we still need the number of times each value is expanded to get the count.
         */ 

        SparseBase(size_t bufsize, const WorkspaceOptions& opt) : 
            report_index(sparse_extract_index(opt.mode)),
            vbuffer(sparse_extract_value(opt.mode) ? bufsize : 0),
            ibuffer(bufsize) 
        {}

        bool report_index;
        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;
    };

    template<bool SPARSE>
    using ConditionalBase = typename std::conditional<SPARSE, SparseBase, DenseBase>::type;

    template<bool WORKROW, bool SPARSE>
    auto define_internal_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        if (opt.mode == SparseExtractMode::VALUE) {
            // Making sure we extract the indices,
            // as we need this to do the reverse mapping of duplicates.
            auto copy = opt;
            copy.mode = SparseExtractMode::BOTH;
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), copy);
        } else if (opt.mode == SparseExtractMode::NONE) {
            // Still need the indices to get the duplicate counts.
            auto copy = opt;
            copy.mode = SparseExtractMode::INDEX;
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), copy);
        } else {
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), opt);
        }
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
    SparseRange<T, IDX> get_sparse(size_t i, T* out_values, IDX* out_indices, InputWorkspace* work) const {
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
            return subset_utils::remap_sparse_duplicates<WORKROW>(mat.get(), i, out_values, out_indices, wptr, wptr->mapping_duplicates, wptr->mapping_duplicates_pool);
        }
    }

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongBlockWorkspace(size_t start, size_t length, size_t bufsize, const WorkspaceOptions& opt, subset_utils::ConditionalSupplement<IDX, SPARSE> host) : 
            BlockWorkspace<WORKROW, SPARSE>(start, length), ConditionalBase2<SPARSE>(bufsize, opt, std::move(host)) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool SPARSE, class Function>
    std::vector<IDX> transplant_indices(subset_utils::ConditionalSupplement<IDX, SPARSE>& temp, size_t length, Function to_index) const {
        if constexpr(SPARSE) {
            size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
            temp.mapping_duplicates.resize(mapping_dim);
            temp.mapping_duplicates_pool.reserve(length);
        } else {
            temp.reverse_mapping.reserve(length);
        }

        std::vector<IDX> local;
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
            subset_utils::ConditionalSupplement<IDX, SPARSE> temp;
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
            // Need to create a temporary object so that we can get the number of unique indices
            // within this block, for use in the constructor for the _actual_ Workspace. 
            subset_utils::ConditionalSupplement<IDX, SPARSE> temp;
            size_t length = subset_.size();
            auto local = transplant_indices<SPARSE>(temp, subset_.size(), [&](size_t i) -> IDX { return subset_[i]; });

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

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, true>(std::move(subset), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
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
