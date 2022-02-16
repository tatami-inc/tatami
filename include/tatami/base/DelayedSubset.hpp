#ifndef TATAMI_DELAYED_SUBSET
#define TATAMI_DELAYED_SUBSET

#include "Matrix.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubset.hpp
 *
 * Delayed subsetting, equivalent to the `DelayedSubset` class in the **DelayedArray** package.
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
    DelayedSubset(std::shared_ptr<const Matrix<T, IDX> > p, V idx) : 
        mat(p), 
        indices(std::move(idx)),
        reverse_indices(MARGIN==1 ? mat->ncol() : mat->nrow(), indices.size())
    {
        for (size_t i = 0; i < indices.size(); ++i) {
            if (i && indices[i] < indices[i-1]) {
                // unsorted, we give up.
                reverse_indices.clear();
                break;
            }

            auto& chosen = reverse_indices[indices[i]];
            if (chosen != indices.size()) {
                // duplicates, we give up.
                reverse_indices.clear(); 
                break;
            }

            chosen = i;
        }
       
        return;
    }

public:
    const T* row(size_t r, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        if constexpr(MARGIN==1) {
            subset_expanded<true>(r, buffer, start, end, work);
            return buffer;
        } else {
            return mat->row(indices[r], buffer, start, end, work);
        }
    }

    const T* column(size_t c, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        if constexpr(MARGIN==1) {
            return mat->column(indices[c], buffer, start, end, work);
        } else {
            subset_expanded<false>(c, buffer, start, end, work);
            return buffer;
        }
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

public:
    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(MARGIN==1) {
            auto total = subset_sparse<true>(r, out_values, out_indices, start, end, work, sorted);
            return SparseRange<T, IDX>(total, out_values, out_indices);
        } else {
            return mat->sparse_row(indices[r], out_values, out_indices, start, end, work, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(MARGIN==1) {
            return mat->sparse_column(indices[c], out_values, out_indices, start, end, work, sorted);
        } else {
            auto total = subset_sparse<false>(c, out_values, out_indices, start, end, work, sorted);
            return SparseRange<T, IDX>(total, out_values, out_indices);
        }
    }

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @return Number of rows after any subsetting is applied.
     */
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }
    
    /**
     * @return Number of columns after any subsetting is applied.
     */
    size_t ncol() const {
        if constexpr(MARGIN==0) {
            return mat->ncol();
        } else {
            return indices.size();
        }
    }

    /**
     * @return The sparsity status of the underlying (pre-subsetted) matrix.
     */
    bool sparse() const {
        return mat->sparse();
    }

    /**
     * @return Whether the underlying (pre-subsetted) matrix prefers row access.
     */
    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    /**
     * @param row Should a workspace be created for row-wise extraction?
     *
     * @return A null pointer or a shared pointer to a `Workspace` object, depending on whether `row` is equal to `MARGIN == 0`.
     */
    std::shared_ptr<Workspace> new_workspace(bool row) const {
        if (row == (MARGIN==0)) {
            return mat->new_workspace(row);
        } else {
            return std::shared_ptr<Workspace>(new SubsetWorkspace(mat.get(), indices, row));
        }
    }
private:
    struct SubsetWorkspace : public Workspace {
        template<class M>
        SubsetWorkspace(const M* ptr, const V& indices, bool row) :
            value_buffer(buffer_size(ptr, row)),
            index_buffer(value_buffer.size()),
            work(ptr->new_workspace(row)) 
        {
            update_last(0, indices.size(), indices);
            return;
        }

        void update_last(size_t start, size_t end, const V& indices) {
            bool changed = false;
            if (start != last_start.first) {
                last_start.first = start;
                changed = true;
            }
            if (end != last_end.first) {
                last_end.first = end;
                changed = true;
            }
            if (changed && start < end) {
                find_min_max(start, end, last_start.second, last_end.second, indices);
            }
            return;
        }

        static void find_min_max(size_t start, size_t end, size_t& min_index, size_t& max_index, const V& indices) {
            min_index = *std::min_element(indices.begin() + start, indices.begin() + end);
            max_index = *std::max_element(indices.begin() + start, indices.begin() + end) + 1;
        }

        template <class M>
        static size_t buffer_size(const M* ptr, bool row) {
            return row ? ptr->ncol() : ptr->nrow();
        }

        std::vector<T> value_buffer;
        std::vector<IDX> index_buffer;
        std::shared_ptr<Workspace> work;

        std::pair<size_t, size_t> last_start, last_end;
    };

private:
    std::shared_ptr<const Matrix<T, IDX> > mat;
    V indices;
    std::vector<IDX> reverse_indices;

    template<bool ROW>
    void subset_expanded(size_t r, T* buffer, size_t start, size_t end, Workspace* work) const {
        if (start >= end) {
            return;
        }

        if (work == NULL) {
            std::vector<T> xbuffer(SubsetWorkspace::buffer_size(mat.get(), ROW));
            size_t min_index, max_index;
            SubsetWorkspace::find_min_max(start, end, min_index, max_index, indices);
            subset_expanded_inner<ROW>(r, buffer, xbuffer.data(), start, end, min_index, max_index, NULL);
        } else {
            auto work0 = reinterpret_cast<SubsetWorkspace*>(work);
            work0->update_last(start, end, indices);
            subset_expanded_inner<ROW>(r, buffer, work0->value_buffer.data(), start, end, work0->last_start.second, work0->last_end.second, work0->work.get());
        }
    }

    template<bool ROW>
    void subset_expanded_inner(size_t r, T* buffer, T* inner_buffer, size_t start, size_t end, size_t min_index, size_t max_index, Workspace* work) const {
        const T* ptr = NULL;
        if constexpr(ROW) {
            ptr = mat->row(r, inner_buffer, min_index, max_index, work);
        } else { 
            ptr = mat->column(r, inner_buffer, min_index, max_index, work);
        }
        for (size_t i = start; i < end; ++i) {
            *buffer = ptr[indices[i] - min_index];
            ++buffer;
        }
        return;
    }

    template<bool ROW>
    size_t subset_sparse(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work, bool sorted) const {
        if (start >= end) {
            return 0;
        }

        if (work == NULL) {
            std::vector<T> xbuffer(SubsetWorkspace::buffer_size(mat.get(), ROW));
            std::vector<IDX> ibuffer(xbuffer.size());
            size_t min_index, max_index;
            SubsetWorkspace::find_min_max(start, end, min_index, max_index, indices);
            return subset_sparse_inner<ROW>(r, out_values, out_indices, xbuffer.data(), ibuffer.data(), start, end, min_index, max_index, NULL, sorted);
        } else {
            auto work0 = reinterpret_cast<SubsetWorkspace*>(work);
            work0->update_last(start, end, indices);
            return subset_sparse_inner<ROW>(r, out_values, out_indices, work0->value_buffer.data(), work0->index_buffer.data(), start, end, work0->last_start.second, work0->last_end.second, work0->work.get(), sorted);
        }
    }

    template<bool ROW>
    size_t subset_sparse_inner(size_t r, T* out_values, IDX* out_indices, T* inner_out_values, IDX* inner_out_indices, size_t start, size_t end, size_t min_index, size_t max_index, Workspace* work, bool sorted) const {
        if (reverse_indices.empty()) {
            // Has duplicates or is out-of-order... need to expand the sparse vector into an array for indexing.
            const T* ptr = NULL;
            if constexpr(ROW) {
                ptr = mat->row(r, inner_out_values, min_index, max_index, work);
            } else { 
                ptr = mat->column(r, inner_out_values, min_index, max_index, work);
            }

            auto copy = out_indices;
            for (size_t i = start; i < end; ++i) {
                auto val = ptr[indices[i] - min_index];
                if (val) { // note that this means that any deliberately inserted zero elements are lost.
                    *out_values = val;
                    ++out_values;
                    *out_indices = i;
                    ++out_indices;
                }
            }

            return static_cast<size_t>(out_indices - copy);

        } else {
            // No duplicates, this just involves doing the reverse look-up.
            SparseRange<T, IDX> range;
            if constexpr(ROW) {
                range = mat->sparse_row(r, inner_out_values, inner_out_indices, min_index, max_index, work, sorted);
            } else {
                range = mat->sparse_column(r, inner_out_values, inner_out_indices, min_index, max_index, work, sorted);
            }

            auto copy = out_indices;
            for (size_t i = 0; i < range.number; ++i, ++range.index, ++range.value) {
                auto final_idx = reverse_indices[*range.index];
                if (final_idx != static_cast<IDX>(indices.size())) {
                    *out_values = *range.value;
                    ++out_values;
                    *out_indices = final_idx;
                    ++out_indices;
                }
            }

            return static_cast<size_t>(out_indices - copy);
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
