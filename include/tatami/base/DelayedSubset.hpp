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
 * @brief Delayed subsetting of a matrix.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix.
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
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `MARGIN = 0`) or columns (if `MARGIN = 1`).
     */
    DelayedSubset(std::shared_ptr<const Matrix<T, IDX> > p, V idx) : mat(std::move(p)), indices(std::move(idx)) {
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            if (i) {
                if (indices[i] < indices[i-1]) {
                    is_unsorted = true;
                    break;
                }
            }
        }

        size_t mapping_dim = MARGIN == 0 ? mat->nrow() : mat->ncol();

        if (!is_unsorted) {
            for (size_t i = 0, end = indices.size(); i < end; ++i) {
                if (i) {
                    if (indices[i] == indices[i-1]) {
                        has_duplicates = true;
                        break;
                    }
                }
            }

            if (!has_duplicates) {
                if constexpr(use_unique_and_sorted) {
                    unique_and_sorted.insert(unique_and_sorted.end(), indices.begin(), indices.end());
                }

                mapping_single.resize(mapping_dim);
                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    mapping_single[indices[i]] = i;
                }

            } else {
                unique_and_sorted.reserve(indices.size());
                reverse_mapping.reserve(indices.size());
                mapping_duplicates.resize(mapping_dim);
                mapping_duplicates_pool.reserve(indices.size());

                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    auto& range = mapping_duplicates[indices[i]];
                    if (unique_and_sorted.empty() || indices[i] != unique_and_sorted.back()) {
                        unique_and_sorted.push_back(indices[i]);
                        range.first = mapping_duplicates_pool.size();
                    }
                    mapping_duplicates_pool.push_back(i);
                    reverse_mapping.push_back(unique_and_sorted.size() - 1);
                    ++range.second;
                }
            }

        } else {
            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(indices.size());
            for (size_t i = 0, end = indices.size(); i < end; ++i) {
                collected.emplace_back(indices[i], i);
            }
            std::sort(collected.begin(), collected.end());

            for (size_t i = 0, end = collected.size(); i < end; ++i) {
                if (i) {
                    if (collected[i].first == collected[i-1].first) {
                        has_duplicates = true;
                        break;
                    }
                }
            }

            unique_and_sorted.reserve(indices.size());
            reverse_mapping.resize(indices.size());

            if (!has_duplicates) {
                for (size_t i = 0, end = collected.size(); i < end; ++i) {
                    const auto& current = collected[i];
                    unique_and_sorted.push_back(current.first);
                    reverse_mapping[current.second] = unique_and_sorted.size() - 1;
                }

                mapping_single.resize(mapping_dim);
                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    mapping_single[indices[i]] = i;
                }

            } else {
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
        }
      
        return;
    }

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    V indices;
    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<V>()[0])>::type>::type V_type;

    bool is_unsorted = false;
    std::vector<size_t> reverse_mapping;

    static constexpr bool use_unique_and_sorted = !std::is_same<V, std::vector<IDX> >::value;
    std::vector<IDX> unique_and_sorted;

    bool has_duplicates = false;
    std::vector<size_t> mapping_single;
    std::vector<std::pair<size_t, size_t> > mapping_duplicates; // holds (position, size) in the pool.
    std::vector<size_t> mapping_duplicates_pool; 

    size_t host_length() const {
        if constexpr(use_unique_and_sorted) {
            return unique_and_sorted.size();
        } else {
            return indices.size();
        }
    }

    const IDX* host_indices() const {
        if constexpr(use_unique_and_sorted) {
            return unique_and_sorted.data();
        } else {
            return indices.data();
        }
    }

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
            size_t buffer_size = (is_unsorted || has_duplicates ? unique_and_sorted.size() : 0);
            return std::shared_ptr<RowWorkspace>(new AlongWorkspace<true>(buffer_size, mat->new_row_workspace(host_length(), host_indices())));
        }
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        if constexpr(MARGIN == 0) {
            size_t buffer_size = (is_unsorted || has_duplicates ? unique_and_sorted.size() : 0);
            return std::shared_ptr<ColumnWorkspace>(new AlongWorkspace<false>(buffer_size, mat->new_column_workspace(host_length(), host_indices())));
        } else {
            return mat->new_column_workspace();
        }
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return mat->row(indices[r], buffer, work);
        } else {
            auto wptr = static_cast<AlongWorkspace<true>*>(work);
            return expand_dense<true>(r, buffer, wptr, reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongWorkspace<false>*>(work);
            return expand_dense<false>(c, buffer, wptr, reverse_mapping);
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
    const T* expand_dense(size_t i, T* buffer, InputWorkspace* work, const std::vector<size_t>& rmapping) const {
        if (!is_unsorted && !has_duplicates) {
            if constexpr(ROW) {
                return mat->row(i, buffer, work->internal.get());
            } else {
                return mat->column(i, buffer, work->internal.get());
            }
        } 

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

    SparseRange<T, IDX> extract_sparse_unique(SparseRange<T, IDX> raw, T* vbuffer, IDX* ibuffer) const {
        auto originalv = vbuffer;
        auto originali = ibuffer;
        for (size_t i = 0; i < raw.number; ++i, ++vbuffer, ++ibuffer) {
            *vbuffer = raw.value[i];
            *ibuffer = mapping_single[raw.index[i]];
        }
        return SparseRange<T, IDX>(raw.number, originalv, originali);
    }

    SparseRange<T, IDX> extract_sparse_duplicates(SparseRange<T, IDX> raw, T* vbuffer, IDX* ibuffer) const {
        auto originalv = vbuffer;
        auto originali = ibuffer;
        size_t counter = 0;

        for (size_t i = 0; i < raw.number; ++i) {
            const auto& pool_pos = mapping_duplicates[raw.index[i]];
            size_t pool_end = pool_pos.first + pool_pos.second;
            for (size_t j = pool_pos.first; j < pool_end; ++j, ++counter, ++vbuffer, ++ibuffer) {
                *vbuffer = raw.value[i];
                *ibuffer = mapping_duplicates_pool[j];
            }
        }

        return SparseRange<T, IDX>(counter, originalv, originali);
    }

    template<bool ROW, class InputWorkspace>
    SparseRange<T, IDX> extract_sparse(size_t i, T* vbuffer, IDX* ibuffer, InputWorkspace* work, bool sorted) const {
        SparseRange<T, IDX> raw;
        if (!is_unsorted && !has_duplicates) {
            if constexpr(ROW) {
                raw = mat->sparse_row(i, vbuffer, ibuffer, work->internal.get(), sorted);
            } else {
                raw = mat->sparse_column(i, vbuffer, ibuffer, work->internal.get(), sorted);
            }
            return extract_sparse_unique(raw, vbuffer, ibuffer);
        }

        if (work->ibuffer.empty()) {
            work->ibuffer.resize(work->vbuffer.size());
        }

        // no need for sorting if it's unsorted, as we'll need to sort again anyway.
        bool raw_sorted = sorted && !is_unsorted;

        if constexpr(ROW) {
            raw = mat->sparse_row(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), raw_sorted);
        } else {
            raw = mat->sparse_column(i, work->vbuffer.data(), work->ibuffer.data(), work->internal.get(), raw_sorted);
        }

        if (!is_unsorted && has_duplicates) {
            return extract_sparse_duplicates(raw, vbuffer, ibuffer);
        }

        if (!sorted) {
            if (!has_duplicates) {
                return extract_sparse_unique(raw, vbuffer, ibuffer);
            } else {
                return extract_sparse_duplicates(raw, vbuffer, ibuffer);
            }
        }

        auto& sortspace = work->sortspace;
        sortspace.clear();

        if (!has_duplicates) {
            sortspace.reserve(raw.number);
            for (size_t i = 0; i < raw.number; ++i) {
                sortspace.emplace_back(mapping_single[raw.index[i]], raw.value[i]);
            }

        } else {
            sortspace.reserve(indices.size());
            for (size_t i = 0; i < raw.number; ++i) {
                const auto& pool_pos = mapping_duplicates[raw.index[i]];
                size_t pool_end = pool_pos.first + pool_pos.second;
                for (size_t j = pool_pos.first; j < pool_end; ++j) {
                    sortspace.emplace_back(mapping_duplicates_pool[j], raw.value[i]);
                }
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

            if (!is_unsorted && !has_duplicates) {
                ptr->internal = mat->new_row_workspace(length, host_indices() + start);
            } else {
                auto& local = ptr->local_unique_and_sorted;
                transplant_indices(local, ptr->local_reverse_mapping, start, length);
                ptr->vbuffer.resize(local.size());
                ptr->internal = mat->new_row_workspace(local.size(), local.data());
            }

            return output;
        }
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        if constexpr(MARGIN == 0) {
            auto ptr = new AlongBlockWorkspace<false>(start, length);
            std::shared_ptr<ColumnBlockWorkspace> output(ptr);

            if (!is_unsorted && !has_duplicates) {
                ptr->internal = mat->new_column_workspace(length, host_indices() + start);
            } else {
                auto& local = ptr->local_unique_and_sorted;
                transplant_indices(local, ptr->local_reverse_mapping, start, length);
                ptr->vbuffer.resize(local.size());
                ptr->internal = mat->new_column_workspace(local.size(), local.data());
            }

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
            return expand_dense<true>(r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return expand_dense<false>(c, buffer, wptr, wptr->local_reverse_mapping);
        } else {
            return mat->column(indices[c], buffer, work);
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(indices[r], out_values, out_indices, work, sorted);
        } else {
            auto wptr = static_cast<AlongBlockWorkspace<true>*>(work);
            return extract_sparse<true>(r, out_values, out_indices, wptr, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongBlockWorkspace<false>*>(work);
            return extract_sparse<false>(c, out_values, out_indices, wptr, sorted);
        } else {
            return mat->sparse_column(indices[c], out_values, out_indices, work, sorted);
        }
    }

private:
    template<class Sortspace>
    void transplant_indices(std::vector<IDX>& local_unique_and_sorted, std::vector<size_t>& local_reverse_mapping, Sortspace& collected) const {
        std::sort(collected.begin(), collected.end());
        local_reverse_mapping.resize(collected.size());

        if (!has_duplicates) {
            for (size_t i = 0, end = collected.size(); i < end; ++i) {
                local_unique_and_sorted.push_back(collected[i].first);
                local_reverse_mapping[collected[i].second] = local_unique_and_sorted.size() - 1;
            }
        } else {
            for (size_t i = 0, end = collected.size(); i < end; ++i) {
                const auto& current = collected[i];
                if (local_unique_and_sorted.empty() || current.first != local_unique_and_sorted.back()) {
                    local_unique_and_sorted.push_back(current.first);
                }
                local_reverse_mapping[current.second] = local_unique_and_sorted.size() - 1;
            }
        }
    }

    void transplant_indices(std::vector<IDX>& local_unique_and_sorted, std::vector<size_t>& local_reverse_mapping, size_t start, size_t length) const {
        local_unique_and_sorted.reserve(length);
        size_t end = start + length;

        if (!is_unsorted) {
            local_reverse_mapping.reserve(length);
            for (size_t i = start; i < end; ++i) {
                if (local_unique_and_sorted.empty() || indices[i] != local_unique_and_sorted.back()) {
                    local_unique_and_sorted.push_back(indices[i]);
                }
                local_reverse_mapping.push_back(local_unique_and_sorted.size() - 1);
            }

        } else {
            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(length);
            for (size_t i = start; i < end; ++i) {
                collected.emplace_back(indices[i], i - start);
            }
            transplant_indices(local_unique_and_sorted, local_reverse_mapping, collected);
        }
    }

public:
    /**
     * @cond
     */
    template<bool ROW>
    struct AlongIndexWorkspace : public IndexWorkspace<IDX, ROW> {
        AlongIndexWorkspace(size_t length, const IDX* indices, size_t n) : IndexWorkspace<IDX, ROW>(length, indices), vbuffer(n) {}

        std::vector<T> vbuffer;
        std::vector<IDX> ibuffer;
        std::vector<std::pair<IDX, T> > sortspace;

        std::vector<IDX> local_unique_and_sorted;
        std::vector<size_t> local_reverse_mapping;
        std::shared_ptr<IndexWorkspace<IDX, ROW> > internal;

    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(size_t length, const IDX* indices) const {
        if constexpr(MARGIN == 0) {
            return mat->new_row_workspace(length, indices);
        } else {
            size_t buffer_size = (is_unsorted || has_duplicates ? length : 0);
            auto ptr = new AlongIndexWorkspace<true>(length, indices, buffer_size);
            std::shared_ptr<RowIndexWorkspace<IDX> > output(ptr);

            auto& local = ptr->local_unique_and_sorted;
            transplant_indices(local, ptr->local_reverse_mapping, length, indices);
            ptr->internal = mat->new_row_workspace(local.size(), local.data());

            return output;
        }
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(size_t length, const IDX* indices) const {
        if constexpr(MARGIN == 0) {
            size_t buffer_size = (is_unsorted || has_duplicates ? length : 0);
            auto ptr = new AlongIndexWorkspace<false>(length, indices, buffer_size);
            std::shared_ptr<ColumnIndexWorkspace<IDX> > output(ptr);

            auto& local = ptr->local_unique_and_sorted;
            transplant_indices(local, ptr->local_reverse_mapping, length, indices);
            ptr->internal = mat->new_column_workspace(local.size(), local.data());

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
            return expand_dense<true>(r, buffer, wptr, wptr->local_reverse_mapping);
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            auto wptr = static_cast<AlongIndexWorkspace<false>*>(work);
            return expand_dense<false>(c, buffer, wptr, wptr->local_reverse_mapping);
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
    void transplant_indices(std::vector<IDX>& local_unique_and_sorted, std::vector<size_t>& local_reverse_mapping, size_t length, const IDX* subset) const {
        local_unique_and_sorted.reserve(length);

        if (!is_unsorted) {
            if (!has_duplicates) {
                for (size_t i = 0; i < length; ++i) {
                    local_unique_and_sorted.push_back(indices[subset[i]]); 
                }

            } else {
                local_reverse_mapping.reserve(length);
                for (size_t i = 0; i < length; ++i) {
                    auto s = subset[i];
                    if (local_unique_and_sorted.empty() || indices[s] != local_unique_and_sorted.back()) {
                        local_unique_and_sorted.push_back(indices[s]);
                    }
                    local_reverse_mapping.push_back(local_unique_and_sorted.size() - 1);
                }
            }

        } else {
            std::vector<std::pair<V_type, size_t> > collected;
            collected.reserve(length);
            for (size_t i = 0; i < length; ++i) {
                collected.emplace_back(indices[subset[i]], i);
            }
            transplant_indices(local_unique_and_sorted, local_reverse_mapping, collected);
        }
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
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
    return std::shared_ptr<MAT>(
        new DelayedSubset<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(
            p,
            std::move(idx)
        )
    );
}

}

#endif
