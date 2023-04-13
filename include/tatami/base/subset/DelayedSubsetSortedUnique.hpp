#ifndef TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetSortedUnique.hpp
 *
 * @brief Delayed subsetting with sorted and unique row/column indices.
 */

namespace tatami {

/**
 * @brief Delayed subsetting of a matrix with sorted, unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of sorted and unique indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam T Type of matrix value.
 * @tparam V Vector containing the subset indices.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX, class V>
class DelayedSubsetSortedUnique : public Matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `MARGIN = 0`) or columns (if `MARGIN = 1`).
     * This should be sorted and unique.
     * @param check Whether to check `idx` for sorted and unique values.
     */
    DelayedSubsetSortedUnique(std::shared_ptr<const Matrix<T, IDX> > p, V idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        if (check) {
            for (size_t i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] <= indices[i-1]) {
                    throw std::runtime_error("indices should be unique and sorted");
                }
            }
        }

        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
        mapping_single.resize(mapping_dim);
        for (IDX i = 0, end = indices.size(); i < end; ++i) {
            mapping_single[indices[i]] = i;
        }
    }

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    V indices;
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

    /*****************************************
     ************ Full extraction ************
     *****************************************/
private:
    template<bool WORKROW, bool SPARSE>
    struct AlongWorkspace : public Workspace<WORKROW, SPARSE> {
        AlongWorkspace(std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > i) : internal(std::move(i)) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<Workspace<WORKROW, SPARSE> > create_new_workspace(const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), opt);
        } else {
            // Deliberate copy of the vector here.
            auto ptr = new AlongWorkspace<WORKROW, SPARSE>(new_workspace<WORKROW, SPARSE>(mat.get(), std::vector<IDX>(indices.begin(), indices.end()), opt));
            return std::shared_ptr<Workspace<WORKROW, SPARSE> >(ptr);
        }
    } 

    template<template<bool, bool> class TargetWorkspace, class InputWorkspace>
    const T* get_dense(size_t i, T* buffer, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_dense<WORKROW>(mat.get(), indices[i], buffer, work);
        } else {
            auto wptr = static_cast<TargetWorkspace<WORKROW, false>*>(work);
            return extract_dense<WORKROW>(mat.get(), i, buffer, wptr->internal.get());
        }
    }

    template<template<bool, bool> class TargetWorkspace, class InputWorkspace>
    SparseRange<T, IDX> get_sparse(size_t i, T* out_values, IDX* out_indices, InputWorkspace* work) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if constexpr((MARGIN==0) == WORKROW) {
            return extract_sparse<WORKROW>(mat.get(), indices[i], out_values, out_indices, work);
        } else {
            auto wptr = static_cast<TargetWorkspace<WORKROW, true>*>(work);
            auto raw = extract_sparse<WORKROW>(mat.get(), i, out_values, out_indices, wptr->internal.get());

            // Only updating the indices if we actually extracted some;
            // otherwise setting the pointer to NULL for the return.
            if (raw.index) {
                auto icopy = out_indices;
                for (size_t i = 0; i < raw.number; ++i, ++icopy) {
                    *icopy = mapping_single[raw.index[i]];
                }
            } else {
                out_indices = NULL;
            }

            return SparseRange<T, IDX>(raw.number, raw.value, out_indices);
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
        return get_dense<AlongWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        return get_dense<AlongWorkspace>(c, buffer, work);
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
    template<bool WORKROW, bool SPARSE>
    struct AlongBlockWorkspace : public BlockWorkspace<WORKROW, SPARSE> {
        AlongBlockWorkspace(size_t start, size_t length, const WorkspaceOptions& opt) :
            BlockWorkspace<WORKROW, SPARSE>(start, length) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), start, length, opt);

        } else {
            auto ptr = new AlongBlockWorkspace<WORKROW, SPARSE>(start, length, opt);
            std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > output(ptr);

            auto left = indices.begin() + start, right = left + length;
            ptr->internal = new_workspace<WORKROW, SPARSE>(mat.get(), std::vector<IDX>(left, right), opt);

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
        return get_dense<AlongBlockWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        return get_dense<AlongBlockWorkspace>(c, buffer, work);
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
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE> {
        AlongIndexWorkspace(std::vector<IDX> subset, const WorkspaceOptions& opt) : 
            IndexWorkspace<IDX, WORKROW, SPARSE>(subset.size()), indices_(std::move(subset)) {}

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

            auto local = ptr->indices_;
            for (auto& x : local) {
                x = indices[x];
            }

            ptr->internal = new_workspace<WORKROW, SPARSE>(mat.get(), std::move(local), opt);
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
        return get_dense<AlongIndexWorkspace>(r, buffer, work);
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        return get_dense<AlongIndexWorkspace>(c, buffer, work);
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
