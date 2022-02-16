#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "Matrix.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedBind.hpp
 *
 * Delayed combining, equivalent to the `DelayedAbind` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed combining of a matrix.
 *
 * Implements delayed combining by rows or columns of a matrix.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam MARGIN Dimension along which the combining is to occur.
 * If 0, the matrices are combined along the rows; if 1, the combining is applied along the columns.
 * @tparam T Type of matrix value.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX>
class DelayedBind : public Matrix<T, IDX> {
public:
    /**
     * @param ps Pointers to the matrices to be combined.
     */
    DelayedBind(std::vector<std::shared_ptr<const Matrix<T, IDX> > > ps) : mats(std::move(ps)), cumulative(mats.size()+1) {
        for (size_t i = 0; i < mats.size(); ++i) {
            if constexpr(MARGIN==0) {
                cumulative[i+1] = mats[i]->nrow();
            } else {
                cumulative[i+1] = mats[i]->ncol();
            }
            cumulative[i+1] += cumulative[i];
        }
    }

    /**
     * @param ps Pointers to the matrices to be combined.
     */
    DelayedBind(const std::vector<std::shared_ptr<Matrix<T, IDX> > >& ps) : mats(ps.begin(), ps.end()), cumulative(mats.size()+1) {
        for (size_t i = 0; i < mats.size(); ++i) {
            if constexpr(MARGIN==0) {
                cumulative[i+1] = mats[i]->nrow();
            } else {
                cumulative[i+1] = mats[i]->ncol();
            }
            cumulative[i+1] += cumulative[i];
        }
    }

public:
    const T* row(size_t r, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        if constexpr(MARGIN==1) {
            assemble_along_dimension<true>(r, buffer, start, end, work);
            return buffer;
        } else {
            return extract_one_dimension<true>(r, buffer, start, end, work);
        }
    }

    const T* column(size_t c, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        if constexpr(MARGIN==1) {
            return extract_one_dimension<false>(c, buffer, start, end, work);
        } else {
            assemble_along_dimension<false>(c, buffer, start, end, work);
            return buffer;
        }
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

private:
    template<bool ROW>
    const T* extract_one_dimension(size_t i, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        size_t chosen = std::upper_bound(cumulative.begin(), cumulative.end(), i) - cumulative.begin() - 1;

        if (work != nullptr) {
            auto work2 = static_cast<BindWorkspace*>(work);
            work = work2->workspaces[chosen].get();
        }

        if constexpr(ROW) {
            return mats[chosen]->row(i - cumulative[chosen], buffer, start, end, work);
        } else {
            return mats[chosen]->column(i - cumulative[chosen], buffer, start, end, work);
        }
    }

    template<bool ROW>
    void assemble_along_dimension(size_t i, T* buffer, size_t start, size_t end, Workspace* work=nullptr) const {
        BindWorkspace* work2 = NULL;
        if (work != nullptr) {
            work2 = static_cast<BindWorkspace*>(work);
        }
        
        size_t left = 0;
        if (start) {
            left = std::upper_bound(cumulative.begin(), cumulative.end(), start) - cumulative.begin() - 1;
        }

        size_t current = start;
        while (current < end) {
            size_t curstart = current - cumulative[left];
            size_t curend = std::min(cumulative[left + 1], end) - cumulative[left];

            Workspace* curwork = NULL;
            if (work2 != NULL) {
                curwork = work2->workspaces[left].get();
            }

            const T * ptr;
            if constexpr(ROW) {
                ptr = mats[left]->row(i, buffer, curstart, curend, curwork);
            } else {
                ptr = mats[left]->column(i, buffer, curstart, curend, curwork);
            }

            size_t filled = curend - curstart;
            if (buffer != ptr) {
                std::copy(ptr, ptr + filled, buffer);
            }
            current += filled;
            buffer += filled;
            ++left;
        }
    }

public:
    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(MARGIN==1) {
            return assemble_along_dimension_sparse<true>(r, out_values, out_indices, start, end, work, sorted);
        } else {
            return extract_one_dimension_sparse<true>(r, out_values, out_indices, start, end, work, sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        if constexpr(MARGIN==1) {
            return extract_one_dimension_sparse<false>(c, out_values, out_indices, start, end, work, sorted);
        } else {
            return assemble_along_dimension_sparse<false>(c, out_values, out_indices, start, end, work, sorted);
        }
    }

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

private:
    template<bool ROW>
    SparseRange<T, IDX> extract_one_dimension_sparse(size_t i, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        size_t chosen = std::upper_bound(cumulative.begin(), cumulative.end(), i) - cumulative.begin() - 1;

        if (work != nullptr) {
            auto work2 = static_cast<BindWorkspace*>(work);
            work = work2->workspaces[chosen].get();
        }

        if constexpr(ROW) {
            return mats[chosen]->sparse_row(i - cumulative[chosen], out_values, out_indices, start, end, work, sorted);
        } else {
            return mats[chosen]->sparse_column(i - cumulative[chosen], out_values, out_indices, start, end, work, sorted);
        }
    }

    template<bool ROW>
    SparseRange<T, IDX> assemble_along_dimension_sparse(size_t i, T* out_values, IDX* out_indices, size_t start, size_t end, Workspace* work=nullptr, bool sorted=true) const {
        BindWorkspace* work2 = NULL;
        if (work != nullptr) {
            work2 = static_cast<BindWorkspace*>(work);
        }
        
        size_t left = 0;
        if (start) {
            left = std::upper_bound(cumulative.begin(), cumulative.end(), start) - cumulative.begin() - 1;
        }

        size_t current = start;
        SparseRange<T, IDX> output(0, out_values, out_indices);
        auto& ntotal = output.number;

        while (current < end) {
            size_t curstart = current - cumulative[left];
            size_t curend = std::min(cumulative[left + 1], end) - cumulative[left];

            Workspace* curwork = NULL;
            if (work2 != NULL) {
                curwork = work2->workspaces[left].get();
            }

            SparseRange<T, IDX> range;
            if constexpr(ROW) {
                range = mats[left]->sparse_row(i, out_values, out_indices, curstart, curend, curwork);
            } else {
                range = mats[left]->sparse_column(i, out_values, out_indices, curstart, curend, curwork);
            }

            if (range.value != out_values) {
                std::copy(range.value, range.value + range.number, out_values);
            }
            if (range.index != out_indices) {
                std::copy(range.index, range.index + range.number, out_indices);
            }
            for (size_t i = 0; i < range.number; ++i) {
                out_indices[i] += cumulative[left];
            }
            ntotal += range.number;

            current += curend - curstart;
            out_indices += range.number;
            out_values += range.number;
            ++left;
        }

        return output;
    }


public:
    /**
     * @return Number of rows after any combining is applied.
     */
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return cumulative.back();
        } else {
            return mats.front()->nrow();
        }
    }
    
    /**
     * @return Number of columns after any combining is applied.
     */
    size_t ncol() const {
        if constexpr(MARGIN==0) {
            return mats.front()->ncol();
        } else {
            return cumulative.back();
        }
    }

public:
    struct BindWorkspace : public Workspace {
        BindWorkspace(std::vector<std::shared_ptr<Workspace> > w) : workspaces(std::move(w)) {}
        std::vector<std::shared_ptr<Workspace> > workspaces;
    };

    /**
     * @param row Should a workspace be created for row-wise extraction?
     *
     * @return A null pointer or a shared pointer to a `Workspace` object.
     */
    std::shared_ptr<Workspace> new_workspace(bool row) const {
        std::vector<std::shared_ptr<Workspace> > workspaces;
        for (const auto& x : mats) {
            workspaces.push_back(x->new_workspace(row));
        }           
        return std::shared_ptr<Workspace>(new BindWorkspace(std::move(workspaces)));
    }

    /**
     * @return The sparsity status of the underlying matrices.
     * If any individual matrix is not sparse, the combination is also considered to be non-sparse.
     */
    bool sparse() const {
        bool is_sparse = true;
        for (const auto& x : mats) {
            is_sparse &= x->sparse();
        }
        return is_sparse;
    }

    /**
     * @return Whether the underlying matrices prefers row access.
     *
     * The preference is defined as the access pattern that is favored by a majority of matrix elements,
     * according to `dimension_preference()`.
     */
    bool prefer_rows() const {
        auto dimpref = dimension_preference();
        return dimpref.first > dimpref.second;
    }

    /**
     * @return A `pair` containing the total number of matrix entries that prefer row-level access (`first`) or column-level access (`second`).
     *
     * The number of entries for each preference is obtained from each matrix, and the values are summed across all of matrices to be combined.
     * This aims to provide `prefer_rows()` with a more nuanced preference that accounts for potential differences in access patterns amongst the underlying matrices.
     */
    std::pair<double, double> dimension_preference() const {
        std::pair<double, double> output;
        for (const auto& x : mats) {
            auto current = x->dimension_preference();
            output.first += current.first;
            output.second += current.second;
        }
        return output;
    }

private:
    std::vector<std::shared_ptr<const Matrix<T, IDX> > > mats;
    std::vector<size_t> cumulative;
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MARGIN Dimension along which the combining is to occur.
 * If 0, matrices are combined along the rows; if 1, matrices are combined to the columns.
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 *
 * @param ps Pointers to `Matrix` objects.
 *
 * @return A pointer to a `DelayedBind` instance.
 */
template<int MARGIN, class MAT>
std::shared_ptr<MAT> make_DelayedBind(std::vector<std::shared_ptr<MAT> > ps) {
    return std::shared_ptr<MAT>(
        new DelayedBind<MARGIN, typename MAT::data_type, typename MAT::index_type>(std::move(ps))
    );
}

}

#endif
