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
    DelayedSubsetSortedUnique(std::shared_ptr<const Matrix<T, IDX> > p, V idx, bool check = true) : mat(std::move(p)), original_indices(std::move(idx)) {
        if (check) {
            for (size_t i = 1, end = original_indices.size(); i < end; ++i) {
                if (original_indices[i] <= original_indices[i-1]) {
                    throw std::runtime_error("indices should be unique and sorted");
                }
            }
        }

        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
        mapping_single.resize(mapping_dim);
        for (IDX i = 0, end = original_indices.size(); i < end; ++i) {
            mapping_single[original_indices[i]] = i;
        }

        if constexpr(use_unique_and_sorted) {
            unique_and_sorted.insert(unique_and_sorted.end(), original_indices.begin(), original_indices.end());
            original_indices.clear(); // freeing up some memory.
        }
    }

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    V original_indices;

    static constexpr bool use_unique_and_sorted = !std::is_same<V, std::vector<IDX> >::value;
    std::vector<IDX> unique_and_sorted;
    std::vector<IDX> mapping_single;

    size_t host_length() const {
        if constexpr(use_unique_and_sorted) {
            return unique_and_sorted.size();
        } else {
            return original_indices.size();
        }
    }

    const IDX* host_indices() const {
        if constexpr(use_unique_and_sorted) {
            return unique_and_sorted.data();
        } else {
            return original_indices.data();
        }
    }

    IDX host_index(size_t x) const {
        if constexpr(use_unique_and_sorted) {
            return unique_and_sorted[x];
        } else {
            return original_indices[x];
        }
    }

public:
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return host_length();
        } else {
            return mat->nrow();
        }
    }

    size_t ncol() const {
        if constexpr(MARGIN==0) {
            return mat->ncol();
        } else {
            return host_length();
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

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongWorkspace : public Workspace<ROW> {
        AlongWorkspace(std::shared_ptr<IndexWorkspace<IDX, ROW> > i) : internal(std::move(i)) {}
        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace();
        } else {
            return std::shared_ptr<RowWorkspace>(new AlongWorkspace<true>(mat->new_row_workspace(host_length(), host_indices())));
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        if constexpr(MARGIN == 0) {
            return std::shared_ptr<ColumnWorkspace>(new AlongWorkspace<false>(mat->new_column_workspace(host_length(), host_indices())));
        } else {
            return mat->new_column_workspace();
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(host_index(r), buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return mat->row(r, buffer, wptr->internal.get());
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return mat->column(c, buffer, wptr->internal.get());
        } else {
            return mat->column(host_index(c), buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(host_index(r), out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            auto raw = mat->sparse_row(r, out_values, out_indices, wptr->internal.get(), sorted);
            return remap_indices(raw, out_indices);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            auto raw = mat->sparse_column(c, out_values, out_indices, wptr->internal.get(), sorted);
            return remap_indices(raw, out_indices);
        } else {
            return mat->sparse_column(host_index(c), out_values, out_indices, work, sorted);
        }
    }

private:
    SparseRange<T, IDX> remap_indices(const SparseRange<T, IDX>& raw, IDX* ibuffer) const {
        auto originali = ibuffer;
        for (size_t i = 0; i < raw.number; ++i, ++ibuffer) {
            *ibuffer = mapping_single[raw.index[i]];
        }
        return SparseRange<T, IDX>(raw.number, raw.value, originali);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongBlockWorkspace : public BlockWorkspace<ROW> {
        AlongBlockWorkspace(size_t start, size_t length, std::shared_ptr<IndexWorkspace<IDX, ROW> > i) : BlockWorkspace<ROW>(start, length), internal(std::move(i)) {}
        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(start, length);
        } else {
            return std::shared_ptr<RowBlockWorkspace>(new AlongBlockWorkspace<true>(start, length, mat->new_row_workspace(length, host_indices() + start)));
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            return std::shared_ptr<ColumnBlockWorkspace>(new AlongBlockWorkspace<false>(start, length, mat->new_column_workspace(length, host_indices() + start)));
        } else {
            return mat->new_column_workspace(start, length);
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(host_index(r), buffer, work);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return mat->row(r, buffer, wptr->internal.get());
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return mat->column(c, buffer, wptr->internal.get());
        } else {
            return mat->column(host_index(c), buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(host_index(r), out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            auto raw = mat->sparse_row(r, out_values, out_indices, wptr->internal.get(), sorted);
            return remap_indices(raw, out_indices);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            auto raw = mat->sparse_column(c, out_values, out_indices, wptr->internal.get(), sorted);
            return remap_indices(raw, out_indices);
        } else {
            return mat->sparse_column(host_index(c), out_values, out_indices, work, sorted);
        }
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        AlongIndexWorkspace(size_t length, const IDX* subset) : IndexWorkspace<IDX, ROW>(length, subset) {}
        std::vector<IDX> local_unique_and_sorted;
        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(size_t length, const IDX* subset) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(length, subset);
        } else {
            auto ptr = new AlongIndexWorkspace<true>(length, subset);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);

            // We can use the local indices directly; if 'subset' is sorted and unique,
            // and 'indices' is sorted and unique, then 'local' must also be sorted and unique. 
            transplant_indices(*ptr, length, subset);
            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_row_workspace(local.size(), local.data());

            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(size_t length, const IDX* subset) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongIndexWorkspace<false>(length, subset);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, subset);
            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_column_workspace(local.size(), local.data());

            return output;
        } else {
            return mat->new_column_workspace(length, subset);
        }
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(host_index(r), buffer, work);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return mat->row(r, buffer, wptr->internal.get());
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return mat->column(c, buffer, wptr->internal.get());
        } else {
            return mat->column(host_index(c), buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(host_index(r), out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            auto raw = mat->sparse_row(r, out_values, out_indices, wptr->internal.get());
            return remap_indices(raw, out_indices);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            auto raw = mat->sparse_column(c, out_values, out_indices, wptr->internal.get());
            return remap_indices(raw, out_indices);
        } else {
            return mat->sparse_column(host_index(c), out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t length, const IDX* subset) const {
        work.local_unique_and_sorted.reserve(length);
        for (size_t i = 0; i < length; ++i) {
            work.local_unique_and_sorted.push_back(host_index(subset[i])); 
        }
    }
};

/**
 * @cond
 */
namespace subset_utils {

template<bool ROW, typename T, typename IDX, class InputWorkspace>
const T* extract_dense(const Matrix<T, IDX>* mat, size_t i, T* buffer, InputWorkspace* work, const std::vector<size_t>& rmapping) {
    const T* dump;
    if constexpr(ROW) {
        dump = mat->row(i, work->vbuffer.data(), work->internal.get());
    } else {
        dump = mat->column(i, work->vbuffer.data(), work->internal.get());
    }

    auto temp = buffer;
    for (auto i : rmapping) {
        *temp = dump[i];
        ++temp;
    } 

    return buffer;
}

template<bool ROW, typename T, typename IDX, class InputWorkspace>
SparseRange<T, IDX> extract_sparse(
    const Matrix<T, IDX>* mat, 
    size_t i, 
    T* vbuffer, 
    IDX* ibuffer, 
    InputWorkspace* work, 
    const std::vector<std::pair<size_t, size_t> >& dups, 
    const std::vector<IDX>& pool,
    bool sorted)
{
    if (work->ibuffer.empty()) {
        work->ibuffer.resize(work->vbuffer.size());
    }

    SparseRange<T, IDX> raw;
    if constexpr(ROW) {
        raw = mat->sparse_row(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), sorted);
    } else {
        raw = mat->sparse_column(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), sorted);
    }

    auto originalv = vbuffer;
    auto originali = ibuffer;
    size_t counter = 0;

    for (size_t i = 0; i < raw.number; ++i) {
        const auto& pool_pos = dups[raw.index[i]];
        size_t pool_end = pool_pos.first + pool_pos.second;
        for (size_t j = pool_pos.first; j < pool_end; ++j, ++counter, ++vbuffer, ++ibuffer) {
            *vbuffer = raw.value[i];
            *ibuffer = pool[j];
        }
    }

    return SparseRange<T, IDX>(counter, originalv, originali);
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

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongWorkspace : public Workspace<ROW> {
        AlongWorkspace(size_t n, std::shared_ptr<IndexWorkspace<IDX, ROW> > i) : vbuffer(n), internal(std::move(i)) {}
        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;
        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace();
        } else {
            size_t buffer_size = unique_and_sorted.size();
            return std::shared_ptr<RowWorkspace>(new AlongWorkspace<true>(buffer_size, mat->new_row_workspace(buffer_size, unique_and_sorted.data())));
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        if constexpr(MARGIN == 0) {
            size_t buffer_size = unique_and_sorted.size();
            return std::shared_ptr<ColumnWorkspace>(new AlongWorkspace<false>(buffer_size, mat->new_column_workspace(buffer_size, unique_and_sorted.data())));
        } else {
            return mat->new_column_workspace();
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return subset_utils::extract_sparse<true>(mat.get(), r, out_values, out_indices, wptr, mapping_duplicates, mapping_duplicates_pool, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return subset_utils::extract_sparse<false>(mat.get(), c, out_values, out_indices, wptr, mapping_duplicates, mapping_duplicates_pool, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongBlockWorkspace : public BlockWorkspace<ROW> {
        AlongBlockWorkspace(size_t start, size_t length) : BlockWorkspace<ROW>(start, length) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;
        std::vector<std::pair<size_t, size_t> > local_mapping_duplicates; 
        std::vector<IDX> local_mapping_duplicates_pool; 

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(start, length);
        } else {
            auto ptr = new AlongBlockWorkspace<true>(start, length);
            std::shared_ptr<RowBlockWorkspace> output(ptr);

            transplant_indices(*ptr, start, length);

            auto& local = ptr->local_unique_and_sorted;
            ptr->vbuffer.resize(local.size());
            ptr->internal = mat->new_row_workspace(local.size(), local.data());

            return output;
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongBlockWorkspace<false>(start, length);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);

            transplant_indices(*ptr, start, length);

            auto& local = ptr->local_unique_and_sorted;
            ptr->vbuffer.resize(local.size());
            ptr->internal = mat->new_column_workspace(local.size(), local.data());

            return output;
        } else {
            return mat->new_column_workspace(start, length);
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return subset_utils::extract_sparse<true>(mat.get(), r, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return subset_utils::extract_sparse<false>(mat.get(), c, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t start, size_t length) const {
        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
        work.local_mapping_duplicates.resize(mapping_dim);
        work.local_mapping_duplicates_pool.reserve(length);
        work.local_unique_and_sorted.reserve(length);
        work.local_reverse_mapping.reserve(length);

        size_t end = start + length;
        for (size_t i = start; i < end; ++i) {
            auto& range = work.local_mapping_duplicates[indices[i]];
            if (work.local_unique_and_sorted.empty() || indices[i] != work.local_unique_and_sorted.back()) {
                work.local_unique_and_sorted.push_back(indices[i]);
                range.first = work.local_mapping_duplicates_pool.size();
            }

            work.local_reverse_mapping.push_back(work.local_unique_and_sorted.size() - 1);
            work.local_mapping_duplicates_pool.push_back(i);
            ++range.second;
        }
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        AlongIndexWorkspace(size_t length, const IDX* indices) : IndexWorkspace<IDX, ROW>(length, indices) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;
        std::vector<std::pair<size_t, size_t> > local_mapping_duplicates; 
        std::vector<IDX> local_mapping_duplicates_pool; 

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(size_t length, const IDX* indices) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(length, indices);
        } else {
            auto ptr = new AlongIndexWorkspace<true>(length, indices);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, indices);

            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_row_workspace(local.size(), local.data());
            ptr->vbuffer.resize(local.size());

            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(size_t length, const IDX* indices) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongIndexWorkspace<false>(length, indices);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, indices);

            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_column_workspace(local.size(), local.data());
            ptr->vbuffer.resize(local.size());

            return output;
        } else {
            return mat->new_column_workspace(length, indices);
        }
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return subset_utils::extract_sparse<true>(mat.get(), r, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return subset_utils::extract_sparse<false>(mat.get(), c, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t length, const IDX* subset) const {
        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
        work.local_mapping_duplicates.resize(mapping_dim);
        work.local_mapping_duplicates_pool.reserve(indices.size());
        work.local_unique_and_sorted.reserve(length);
        work.local_reverse_mapping.reserve(length);

        for (size_t i = 0; i < length; ++i) {
            auto s = subset[i];
            auto& range = work.local_mapping_duplicates[indices[s]];
            if (work.local_unique_and_sorted.empty() || indices[s] != work.local_unique_and_sorted.back()) {
                work.local_unique_and_sorted.push_back(indices[s]);
                range.first = work.local_mapping_duplicates_pool.size();
            }
            work.local_reverse_mapping.push_back(work.local_unique_and_sorted.size() - 1);
            work.local_mapping_duplicates_pool.push_back(s);
            ++range.second;
        }
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

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongWorkspace : public Workspace<ROW> {
        AlongWorkspace(size_t n, std::shared_ptr<IndexWorkspace<IDX, ROW> > i) : vbuffer(n), internal(std::move(i)) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;

        std::vector<std::pair<IDX, T> > sortspace;

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace();
        } else {
            size_t buffer_size = unique_and_sorted.size();
            return std::shared_ptr<RowWorkspace>(new AlongWorkspace<true>(buffer_size, mat->new_row_workspace(buffer_size, unique_and_sorted.data())));
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        if constexpr(MARGIN == 0) {
            size_t buffer_size = unique_and_sorted.size();
            return std::shared_ptr<ColumnWorkspace>(new AlongWorkspace<false>(buffer_size, mat->new_column_workspace(buffer_size, unique_and_sorted.data())));
        } else {
            return mat->new_column_workspace();
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return extract_sparse<true>(r, out_values, out_indices, wptr, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return extract_sparse<false>(c, out_values, out_indices, wptr, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<bool ROW, class InputWorkspace>
    SparseRange<T, IDX> extract_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work, bool sorted) const {
        SparseRange<T, IDX> raw;

        if (work->ibuffer.empty()) {
            work->ibuffer.resize(work->vbuffer.size());
        }

        // no need for sorting if it's unsorted, as we'll need to sort again anyway.
        if constexpr(ROW) {
            raw = mat->sparse_row(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), false);
        } else {
            raw = mat->sparse_column(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), false);
        }

        if (!sorted) {
            auto originali = ibuffer;
            auto originalv = vbuffer;
            for (size_t i = 0; i < raw.number; ++i, ++ibuffer, ++vbuffer) {
                *ibuffer = mapping_single[raw.index[i]];
                *vbuffer = raw.value[i];
            }
            return SparseRange<T, IDX>(raw.number, originalv, originali);
        }

        auto& sortspace = work->sortspace;
        sortspace.clear();
        sortspace.reserve(raw.number);
        for (size_t i = 0; i < raw.number; ++i) {
            sortspace.emplace_back(mapping_single[raw.index[i]], raw.value[i]);
        }

        auto originalv = vbuffer;
        auto originali = ibuffer;

        std::sort(sortspace.begin(), sortspace.end());
        for (const auto& x : sortspace) {
            *vbuffer = x.second;
            *ibuffer = x.first;
            ++vbuffer;
            ++ibuffer;
        }

        return SparseRange<T, IDX>(raw.number, originalv, originali);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongBlockWorkspace : public BlockWorkspace<ROW> {
        AlongBlockWorkspace(size_t start, size_t length) : BlockWorkspace<ROW>(start, length) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;

        std::vector<std::pair<IDX, T> > sortspace;

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(start, length);
        } else {
            auto ptr = new AlongBlockWorkspace<true>(start, length);
            std::shared_ptr<RowBlockWorkspace> output(ptr);

            transplant_indices(*ptr, start, length);
            auto& local = ptr->local_unique_and_sorted;
            ptr->vbuffer.resize(local.size());
            ptr->internal = mat->new_row_workspace(local.size(), local.data());

            return output;
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongBlockWorkspace<false>(start, length);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);

            transplant_indices(*ptr, start, length);
            auto& local = ptr->local_unique_and_sorted;
            ptr->vbuffer.resize(local.size());
            ptr->internal = mat->new_column_workspace(local.size(), local.data());

            return output;
        } else {
            return mat->new_column_workspace(start, length);
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            return extract_sparse<true>(r, out_values, out_indices, static_cast<AlongBlockWorkspace<true>*>(work), sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return extract_sparse<false>(c, out_values, out_indices, static_cast<AlongBlockWorkspace<false>*>(work), sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work) const {
        auto& collected = work.sortspace;
        std::sort(collected.begin(), collected.end());

        work.local_unique_and_sorted.reserve(collected.size());
        work.local_reverse_mapping.resize(collected.size());
        for (size_t i = 0, end = collected.size(); i < end; ++i) {
            work.local_unique_and_sorted.push_back(collected[i].first);
            work.local_reverse_mapping[collected[i].second] = work.local_unique_and_sorted.size() - 1;
        }
    }

    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t start, size_t length) const {
        size_t end = start + length;

        auto& collected = work.sortspace;
        collected.reserve(length);
        for (size_t i = start; i < end; ++i) {
            collected.emplace_back(indices[i], i - start);
        }

        transplant_indices(work);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        AlongIndexWorkspace(size_t length, const IDX* indices) : IndexWorkspace<IDX, ROW>(length, indices) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;

        std::vector<std::pair<IDX, T> > sortspace;

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(size_t length, const IDX* subset) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(length, subset);
        } else {
            auto ptr = new AlongIndexWorkspace<true>(length, subset);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, subset);
            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_row_workspace(local.size(), local.data());
            ptr->vbuffer.resize(local.size());

            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(size_t length, const IDX* subset) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongIndexWorkspace<false>(length, subset);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, subset);
            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_column_workspace(local.size(), local.data());
            ptr->vbuffer.resize(local.size());

            return output;
        } else {
            return mat->new_column_workspace(length, subset);
        }
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return extract_sparse<true>(r, out_values, out_indices, wptr, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return extract_sparse<false>(c, out_values, out_indices, wptr, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t length, const IDX* subset) const {
        auto& collected = work.sortspace;
        collected.reserve(length);
        for (size_t i = 0; i < length; ++i) {
            collected.emplace_back(indices[subset[i]], i);
        }
        transplant_indices(work);
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

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongWorkspace : public Workspace<ROW> {
        AlongWorkspace(size_t n, std::shared_ptr<IndexWorkspace<IDX, ROW> > i) : vbuffer(n), internal(std::move(i)) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;
        std::vector<std::pair<IDX, T> > sortspace;

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace();
        } else {
            size_t buffer_size = unique_and_sorted.size();
            return std::shared_ptr<RowWorkspace>(new AlongWorkspace<true>(buffer_size, mat->new_row_workspace(buffer_size, unique_and_sorted.data())));
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        if constexpr(MARGIN == 0) {
            size_t buffer_size = unique_and_sorted.size();
            return std::shared_ptr<ColumnWorkspace>(new AlongWorkspace<false>(buffer_size, mat->new_column_workspace(buffer_size, unique_and_sorted.data())));
        } else {
            return mat->new_column_workspace();
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return extract_sparse<true>(r, out_values, out_indices, wptr, mapping_duplicates, mapping_duplicates_pool, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return extract_sparse<false>(c, out_values, out_indices, wptr, mapping_duplicates, mapping_duplicates_pool, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<bool ROW, class InputWorkspace>
    SparseRange<T, IDX> extract_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work, const std::vector<std::pair<size_t, size_t> >& dups, const std::vector<IDX>& pool, bool sorted) const {
        if (!sorted) {
            return subset_utils::extract_sparse<ROW>(mat.get(), i, vbuffer, ibuffer, work, dups, pool, sorted);
        }

        SparseRange<T, IDX> raw;
        if (work->ibuffer.empty()) {
            work->ibuffer.resize(work->vbuffer.size());
        }

        // no need for sorting if it's unsorted, as we'll need to sort again anyway.
        if constexpr(ROW) {
            raw = mat->sparse_row(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), false);
        } else {
            raw = mat->sparse_column(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), false);
        }

        auto& sortspace = work->sortspace;
        sortspace.clear();
        sortspace.reserve(indices.size());
        for (size_t i = 0; i < raw.number; ++i) {
            const auto& pool_pos = dups[raw.index[i]];
            size_t pool_end = pool_pos.first + pool_pos.second;
            for (size_t j = pool_pos.first; j < pool_end; ++j) {
                sortspace.emplace_back(pool[j], raw.value[i]);
            }
        }

        auto originalv = vbuffer;
        auto originali = ibuffer;

        std::sort(sortspace.begin(), sortspace.end());
        for (const auto& x : sortspace) {
            *vbuffer = x.second;
            *ibuffer = x.first;
            ++vbuffer;
            ++ibuffer;
        }

        return SparseRange<T, IDX>(sortspace.size(), originalv, originali);
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongBlockWorkspace : public BlockWorkspace<ROW> {
        AlongBlockWorkspace(size_t start, size_t length) : BlockWorkspace<ROW>(start, length) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;
        std::vector<std::pair<IDX, T> > sortspace;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;
        std::vector<std::pair<size_t, size_t> > local_mapping_duplicates; 
        std::vector<IDX> local_mapping_duplicates_pool; 

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(start, length);
        } else {
            auto ptr = new AlongBlockWorkspace<true>(start, length);
            std::shared_ptr<RowBlockWorkspace> output(ptr);

            transplant_indices(*ptr, start, length);
            auto& local = ptr->local_unique_and_sorted;
            ptr->vbuffer.resize(local.size());
            ptr->internal = mat->new_row_workspace(local.size(), local.data());

            return output;
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongBlockWorkspace<false>(start, length);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);

            transplant_indices(*ptr, start, length);
            auto& local = ptr->local_unique_and_sorted;
            ptr->vbuffer.resize(local.size());
            ptr->internal = mat->new_column_workspace(local.size(), local.data());

            return output;
        } else {
            return mat->new_column_workspace(start, length);
        }
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return extract_sparse<true>(r, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return extract_sparse<false>(c, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace, class Function>
    void transplant_indices(InputWorkspace& work, Function to_index) const {
        auto& collected = work.sortspace;
        std::sort(collected.begin(), collected.end());

        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();
        work.local_mapping_duplicates.resize(mapping_dim);
        work.local_mapping_duplicates_pool.reserve(collected.size());
        work.local_unique_and_sorted.reserve(collected.size());
        work.local_reverse_mapping.resize(collected.size());

        for (size_t i = 0, end = collected.size(); i < end; ++i) {
            const auto& current = collected[i];
            auto& range = work.local_mapping_duplicates[current.first];
            if (work.local_unique_and_sorted.empty() || current.first != work.local_unique_and_sorted.back()) {
                work.local_unique_and_sorted.push_back(current.first);
                range.first = work.local_mapping_duplicates_pool.size();
            }

            work.local_reverse_mapping[current.second] = work.local_unique_and_sorted.size() - 1;
            work.local_mapping_duplicates_pool.push_back(to_index(current.second));
            ++range.second;
        }
    }

    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t start, size_t length) const {
        size_t end = start + length;
        auto& collected = work.sortspace;
        collected.reserve(length);
        for (size_t i = start; i < end; ++i) {
            collected.emplace_back(indices[i], i - start);
        }
        transplant_indices(work, [&](size_t i) -> size_t { return i + start; });
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        AlongIndexWorkspace(size_t length, const IDX* indices) : IndexWorkspace<IDX, ROW>(length, indices) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;
        std::vector<std::pair<size_t, size_t> > local_mapping_duplicates; 
        std::vector<IDX> local_mapping_duplicates_pool; 

        std::vector<std::pair<IDX, T> > sortspace;

        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(size_t length, const IDX* indices) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(length, indices);
        } else {
            auto ptr = new AlongIndexWorkspace<true>(length, indices);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, indices);

            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_row_workspace(local.size(), local.data());
            ptr->vbuffer.resize(local.size());

            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(size_t length, const IDX* indices) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongIndexWorkspace<false>(length, indices);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);

            transplant_indices(*ptr, length, indices);

            auto& local = ptr->local_unique_and_sorted;
            ptr->internal = mat->new_column_workspace(local.size(), local.data());
            ptr->vbuffer.resize(local.size());

            return output;
        } else {
            return mat->new_column_workspace(length, indices);
        }
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return subset_utils::extract_dense<true>(mat.get(), r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return subset_utils::extract_dense<false>(mat.get(), c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongIndexWorkspace<true>*>(work);
            return extract_sparse<true>(r, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return extract_sparse<false>(c, out_values, out_indices, wptr, wptr->local_mapping_duplicates, wptr->local_mapping_duplicates_pool, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class InputWorkspace>
    void transplant_indices(InputWorkspace& work, size_t length, const IDX* subset) const {
        auto& collected = work.sortspace;
        collected.reserve(length);
        for (size_t i = 0; i < length; ++i) {
            collected.emplace_back(indices[subset[i]], i);
        }
        transplant_indices(work, [&](size_t i) -> size_t { return subset[i]; });
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
