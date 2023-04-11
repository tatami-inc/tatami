#ifndef TATAMI_DELAYED_SUBSET
#define TATAMI_DELAYED_SUBSET

#include "Matrix.hpp"
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

/**
 * @cond
 */
namespace subset_utils {

template<bool WORKROW, typename T, typename IDX, class InputWorkspace>
const T* remap_dense(const Matrix<T, IDX>* mat, size_t i, T* buffer, InputWorkspace* work, const std::vector<size_t>& rmapping) {
    const T* dump = extract_dense<WORKROW>(mat, i, work->buffer.data(), work->internal.get());
    auto temp = buffer;
    for (auto i : rmapping) {
        *temp = dump[i];
        ++temp;
    } 
    return buffer;
}

template<bool WORKROW, typename T, typename IDX, class InputWorkspace>
SparseRange<T, IDX> remap_sparse_duplicates(
    const Matrix<T, IDX>* mat, 
    size_t i, 
    T* vbuffer, 
    IDX* ibuffer, 
    InputWorkspace* work, 
    const std::vector<std::pair<size_t, size_t> >& dups, 
    const std::vector<IDX>& pool)
{
    T* vin = work->vbuffer.data();
    IDX* iin = work->ibuffer.data();
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

    return SparseRange<T, IDX>(counter, vbuffer, ibuffer);
}

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
         * In general, we can't extract into the user-supplied ibuffer/vbuffer, because we need some
         * place to hold the to-be-processed values while we expand the duplicates.
         *
         * This can be subjected to some optimizations depending on the extraction mode:
         *
         * - NONE: both vbuffer and ibuffer are empty. The internal workspace's mode is mirrored to NONE, 
         *   so it won't attempt to use vbuffer/ibuffer, and remap_sparse_duplicates will return early.
         * - INDEX: vbuffer is empty and ibuffer is allocated. The internal workspace's mode is mirrored to INDEX, 
         *   so it won't attempt to use vbuffer. remap_sparse_duplicates will also skip the null'd value.
         *
         * For VALUE, both vbuffer and ibuffer are allocated, as we need to go use the indices to figure out
         * how much to expand each value. remap_sparse_duplicates will ignore the null'd index, though.
         *
         * For BOTH, both vbuffer and ibuffer are allocated, obviously.
         */ 

        SparseBase(size_t bufsize, const WorkspaceOptions& opt) : 
            report_index(sparse_extract_index(opt.mode)),
            vbuffer(sparse_extract_value(opt.mode) ? bufsize : 0),
            ibuffer(opt.mode != SparseExtractMode::NONE ? bufsize : 0) 
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
            // Making sure we extract the indices if we want the values,
            // as we need this to do the reverse mapping of duplicates.
            auto copy = opt;
            copy.mode = SparseExtractMode::BOTH;
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), copy);
        } else {
            return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), opt);
        }
    }

public:
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

private:
    struct DenseSupplement {
        std::vector<size_t> reverse_mapping;
        static constexpr bool sparse = false;
    };

    struct SparseSupplement {
        std::vector<std::pair<size_t, size_t> > mapping_duplicates; 
        std::vector<IDX> mapping_duplicates_pool; 
        static constexpr bool sparse = true;
    };

    template<bool SPARSE>
    using ConditionalSupplement = typename std::conditional<SPARSE, SparseSupplement, DenseSupplement>::type;

    struct DenseBase2 : public DenseBase, public DenseSupplement {
        DenseBase2(size_t bufsize, const WorkspaceOptions& opt, DenseSupplement host) : DenseBase(bufsize, opt), DenseSupplement(std::move(host)) {}
    };

    struct SparseBase2 : public SparseBase, public SparseSupplement {
        SparseBase2(size_t bufsize, const WorkspaceOptions& opt, SparseSupplement host) : SparseBase(bufsize, opt), SparseSupplement(std::move(host)) {}
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
        AlongBlockWorkspace(size_t start, size_t length, size_t bufsize, const WorkspaceOptions& opt, ConditionalSupplement<SPARSE> host) : 
            BlockWorkspace<WORKROW, SPARSE>(start, length), ConditionalBase2<SPARSE>(bufsize, opt, std::move(host)) {}

        std::shared_ptr<IndexWorkspace<IDX, WORKROW, SPARSE> > internal;
    };

    template<bool WORKROW, bool SPARSE>
    std::shared_ptr<BlockWorkspace<WORKROW, SPARSE> > create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            return new_workspace<WORKROW, SPARSE>(mat.get(), start, length, opt);

        } else {
            // Need to create a temporary object so that we can get the number of unique indices
            // within this block, for use in the constructor for the _actual_ Workspace. 
            ConditionalSupplement<SPARSE> temp;

            if constexpr(SPARSE) {
                size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
                temp.mapping_duplicates.resize(mapping_dim);
                temp.mapping_duplicates_pool.reserve(length);
            } else {
                temp.reverse_mapping.reserve(length);
            }

            std::vector<IDX> local;
            local.reserve(length);

            size_t end = start + length;
            for (size_t i = start; i < end; ++i) {
                bool diff = local.empty() || indices[i] != local.back();
                if (diff) {
                    local.push_back(indices[i]);
                }

                if constexpr(SPARSE) {
                    auto& range = temp.mapping_duplicates[indices[i]];
                    if (diff) {
                        range.first = temp.mapping_duplicates_pool.size();
                    }
                    ++range.second;
                    temp.mapping_duplicates_pool.push_back(i);
                } else {
                    temp.reverse_mapping.push_back(local.size() - 1);
                }
            }

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

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongIndexWorkspace(std::vector<IDX> i, size_t bufsize, const WorkspaceOptions& opt, ConditionalSupplement<SPARSE> host) : 
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
            ConditionalSupplement<SPARSE> temp;

            size_t length = subset_.size();
            if constexpr(SPARSE) {
                size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
                temp.mapping_duplicates.resize(mapping_dim);
                temp.mapping_duplicates_pool.reserve(indices.size());
            } else {
                temp.reverse_mapping.reserve(length);
            }

            std::vector<IDX> local;
            local.reserve(length);

            for (size_t i = 0; i < length; ++i) {
                auto s = subset_[i];
                bool diff = (local.empty() || indices[s] != local.back());
                if (diff) {
                    local.push_back(indices[s]);
                }

                if constexpr(SPARSE) {
                    auto& range = temp.mapping_duplicates[indices[s]];
                    if (diff) {
                        range.first = temp.mapping_duplicates_pool.size();
                    }
                    ++range.second;
                    temp.mapping_duplicates_pool.push_back(s);
                } else {
                    temp.reverse_mapping.push_back(local.size() - 1);
                }
            }

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
         * There's no need to have a vbuffer here, because we can directly extract into
         * the user-supplied vbuffer if we want the values. The inner extraction is guaranteed
         * to be a subset of the vbuffer length, as we uniquify the indices in the internal workspace.
         * Overwriting is not a concern because we end up copying everything into sortspace for re-sorting anyway.
         *
         * Similarly, an ibuffer is _usually_ unnecessary as we can write directly into the user-supplied ibuffer.
         * There's only a need to allocate our own ibuffer if (i) the user wants the values but not 
         * the indices, and (ii) they want the output to be sorted. In this case, the 
         * user-supplied ibuffer may not be valid, so we need to create our own.
         */ 

        SparseBase(size_t bufsize, const WorkspaceOptions& opt) : 
            report_index(sparse_extract_index(opt.mode)),
            needs_sort(opt.sorted),
            ibuffer(opt.mode == SparseExtractMode::VALUE && needs_sort ? bufsize : 0) 
        {}

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
        sortspace.reserve(raw.number);

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
         * If we want sorted output, then the logic follows that of DelayedSubsetUnique::SparseBase.
         * We can directly extract into the user-supplied ibuffer/vbuffer and then copy everything into
         * sortspace for resorting. The only exception is if the user wants the values but not the
         * indices, in which case we need to do some allocation.
         *
         * If we don't want sorted output, then the logic follows that of DelayedSubsetSorted::SparseBase.
         * Here, we can't extract into the user-supplied ibuffer/vbuffer, because we need some
         * place to hold the to-be-processed values while we expand the duplicates.
         * This can be subjected to some optimizations depending on the extraction mode and sorting:
         *
         * - NONE: both vbuffer and ibuffer are empty. The internal workspace's mode is mirrored to NONE, 
         *   so it won't attempt to use vbuffer/ibuffer, and reorganize_sparse will return early.
         * - INDEX: vbuffer is empty and ibuffer is allocated. The internal workspace's mode is mirrored to INDEX, 
         *   so it won't attempt to use vbuffer. remap_sparse_duplicates will also skip the null'd value.
         *
         * For VALUE, both vbuffer and ibuffer are allocated, as we need to go use the indices to figure out
         * how much to expand each value. remap_sparse_duplicates will ignore the null'd index, though.
         *
         * For BOTH, both vbuffer and ibuffer are allocated, obviously.
         */ 

        SparseBase(size_t bufsize, const WorkspaceOptions& opt) : report_index(sparse_extract_index(opt.mode)), needs_sort(opt.sorted) {
            if (needs_sort) {
                if (opt.mode == SparseExtractMode::VALUE) {
                    ibuffer.resize(bufsize);
                }
            } else {
                if (sparse_extract_value(opt.mode)) {
                    vbuffer.resize(bufsize);
                }
                if (opt.mode != SparseExtractMode::NONE) {
                    ibuffer.resize(bufsize);
                }
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
        if (opt.mode == SparseExtractMode::VALUE) {
            // Making sure we extract the indices to do the deduplication.
            opt.mode = SparseExtractMode::BOTH;
        }

        // Turning off the sorting to enable possible optimizations in the underlying matrix.
        // We don't need sorted output as we'll be resorting ourselves later.
        opt.sorted = false;
        return new_workspace<WORKROW, SPARSE>(mat.get(), std::move(subset), opt);
    }

    template<class InputWorkspace>
    SparseRange<T, IDX> reorganize_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work, const std::vector<std::pair<size_t, size_t> >& dups, const std::vector<IDX>& pool) const {
        constexpr bool WORKROW = InputWorkspace::row;
        if (!work->needs_sort) {
            return subset_utils::remap_sparse_duplicates<WORKROW>(mat.get(), i, vbuffer, ibuffer, work, dups, pool);
        }

        IDX* iin = (work->ibuffer.empty() ? ibuffer : work->ibuffer.data());
        auto raw = extract_sparse<WORKROW>(mat.get(), i, vbuffer, iin, work->internal.get());

        auto& sortspace = work->sortspace;
        sortspace.clear();
        sortspace.reserve(indices.size());

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

private:
    struct DenseSupplement {
        std::vector<size_t> reverse_mapping;
        static constexpr bool sparse = false;
    };

    struct SparseSupplement {
        std::vector<std::pair<size_t, size_t> > mapping_duplicates; 
        std::vector<IDX> mapping_duplicates_pool; 
        static constexpr bool sparse = true;
    };

    template<bool SPARSE>
    using ConditionalSupplement = typename std::conditional<SPARSE, SparseSupplement, DenseSupplement>::type;

    struct DenseBase2 : public DenseBase, public DenseSupplement {
        DenseBase2(size_t bufsize, const WorkspaceOptions& opt, DenseSupplement host) : DenseBase(bufsize, opt), DenseSupplement(std::move(host)) {}
    };

    struct SparseBase2 : public SparseBase, public SparseSupplement {
        SparseBase2(size_t bufsize, const WorkspaceOptions& opt, SparseSupplement host) : SparseBase(bufsize, opt), SparseSupplement(std::move(host)) {}
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
        AlongBlockWorkspace(size_t s, size_t l, size_t bufsize, const WorkspaceOptions& opt, ConditionalSupplement<SPARSE> host) : 
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
            ConditionalSupplement<SPARSE> temp;
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

private:
    template<bool WORKROW, bool SPARSE>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, WORKROW, SPARSE>, public ConditionalBase2<SPARSE> {
        AlongIndexWorkspace(std::vector<IDX> i, size_t bufsize, const WorkspaceOptions& opt, ConditionalSupplement<SPARSE> host) : 
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
            ConditionalSupplement<SPARSE> temp;
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

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 * This will automatically dispatch to `DelayedSubsetSortedUnique`, `DelayedSubsetUnique`, `DelayedSubsetSorted` or `DelayedSubset`, depending on the values in `idx`.
 *
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 * @tparam V Vector containing the subset indices, to be automatically deducted.
 *
 * @param p Pointer to a `Matrix`.
 * @param idx Instance of the index vector.
 *
 * @return A pointer to a `DelayedSubset` instance.
 */
template<int MARGIN, class MAT, class V>
std::shared_ptr<MAT> make_DelayedSubset(std::shared_ptr<MAT> p, V idx) {
    bool is_unsorted = false;
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        if (i) {
            if (idx[i] < idx[i-1]) {
                is_unsorted = true;
                break;
            }
        }
    }

    if (!is_unsorted) {
        bool has_duplicates = false;
        for (size_t i = 0, end = idx.size(); i < end; ++i) {
            if (i) {
                if (idx[i] == idx[i-1]) {
                    has_duplicates = true;
                    break;
                }
            }
        }

        if (!has_duplicates) {
            return std::shared_ptr<MAT>(
                new DelayedSubsetSortedUnique<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), std::move(idx), false)
            );
        } else {
            return std::shared_ptr<MAT>(
                new DelayedSubsetSorted<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), std::move(idx), false)
            );
        }
    }

    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<V>()[0])>::type>::type V_type;
    std::vector<std::pair<V_type, size_t> > collected;
    collected.reserve(idx.size());
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        collected.emplace_back(idx[i], i);
    }
    std::sort(collected.begin(), collected.end());

    bool has_duplicates = false;
    for (size_t i = 1, end = collected.size(); i < end; ++i) {
        if (collected[i].first == collected[i-1].first) {
            has_duplicates = true;
            break;
        }
    }

    if (!has_duplicates) {
        return std::shared_ptr<MAT>(
            new DelayedSubsetUnique<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), collected, std::move(idx))
        );
    } else {
        return std::shared_ptr<MAT>(
            new DelayedSubset<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), collected, std::move(idx))
        );
    }
}

}

#endif
