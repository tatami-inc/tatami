#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "Matrix.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedBind.hpp
 *
 * @brief Delayed combining of multiple `tatami::Matrix` objects.
 *
 * This is equivalent to the `DelayedAbind` class in the **DelayedArray** package.
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

        set_constants();
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

        set_constants();
    }

private:
    std::vector<std::shared_ptr<const Matrix<T, IDX> > > mats;
    std::vector<size_t> cumulative;

    bool sparse_;
    std::pair<double, double> dimension_preference_;

    void set_constants() {
        sparse_ = true;
        for (const auto& x : mats) {
            if (!(x->sparse())) {
                sparse_ = false;
                break;
            }
        }

        dimension_preference_.first = 0;
        dimension_preference_.second = 0;
        for (const auto& x : mats) {
            auto current = x->dimension_preference();
            dimension_preference_.first += current.first;
            dimension_preference_.second += current.second;
        }
    }

public:
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return cumulative.back();
        } else {
            if (mats.empty()) {
                return 0;
            } else {
                return mats.front()->nrow();
            }
        }
    }

    size_t ncol() const {
        if constexpr(MARGIN==0) {
            if (mats.empty()) {
                return 0;
            } else {
                return mats.front()->ncol();
            }
        } else {
            return cumulative.back();
        }
    }

    bool sparse() const {
        return sparse_;
    }

    bool prefer_rows() const {
        const auto& dimpref = dimension_preference_;
        return dimpref.first > dimpref.second;
    }

    std::pair<double, double> dimension_preference() const {
        return dimension_preference_;
    }

    using Matrix<T, IDX>::column;

    using Matrix<T, IDX>::row;

    using Matrix<T, IDX>::sparse_column;

    using Matrix<T, IDX>::sparse_row;

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct PerpendicularWorkspace : public Workspace<WORKROW> {
        PerpendicularWorkspace(std::vector<std::shared_ptr<Workspace<WORKROW> > > w) : workspaces(std::move(w)) {}
        std::vector<std::shared_ptr<Workspace<WORKROW> > > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW>
    struct ParallelWorkspace : public Workspace<WORKROW> {
        ParallelWorkspace(std::vector<std::shared_ptr<Workspace<WORKROW> > > w) : workspaces(std::move(w)) {}
        std::vector<std::shared_ptr<Workspace<WORKROW> > > workspaces;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowWorkspace> new_row_workspace() const {
        return new_workspace<true>();
    }

    std::shared_ptr<ColumnWorkspace> new_column_workspace() const {
        return new_workspace<false>();
    }

    const T* row(size_t r, T* buffer, RowWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, buffer, static_cast<PerpendicularWorkspace<true>*>(work));
        } else {
            assemble_along_dimension_simple<true>(r, buffer, static_cast<ParallelWorkspace<true>*>(work));
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, ColumnWorkspace* work) const {
        if constexpr(MARGIN==0) {
            assemble_along_dimension_simple<false>(c, buffer, static_cast<ParallelWorkspace<false>*>(work));
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, static_cast<PerpendicularWorkspace<false>*>(work));
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, static_cast<PerpendicularWorkspace<true>*>(work), sorted);
        } else {
            return assemble_along_dimension_simple<true>(r, out_values, out_indices, static_cast<ParallelWorkspace<true>*>(work), sorted);
        }
    }

    SparseRange<T, IDX> simple_column(size_t c, T* out_values, IDX* out_indices, ColumnWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return assemble_along_dimension_simple<false>(c, out_values, out_indices, static_cast<ParallelWorkspace<false>*>(work), sorted);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, static_cast<PerpendicularWorkspace<false>*>(work), sorted);
        }
    }

private:
    template<bool WORKROW>
    std::shared_ptr<Workspace<WORKROW> > new_workspace() const {
        std::vector<std::shared_ptr<Workspace<WORKROW> > > workspaces;
        workspaces.reserve(mats.size());
        for (const auto& x : mats) {
            if constexpr(WORKROW) {
                workspaces.push_back(x->new_row_workspace());
            } else {
                workspaces.push_back(x->new_column_workspace());
            }
        }

        if constexpr((MARGIN == 0) == WORKROW) {
            return std::shared_ptr<Workspace<WORKROW> >(new PerpendicularWorkspace<WORKROW>(std::move(workspaces)));
        } else {
            return std::shared_ptr<Workspace<WORKROW> >(new ParallelWorkspace<WORKROW>(std::move(workspaces)));
        }
    }

    size_t choose_segment(size_t i) const {
        return std::upper_bound(cumulative.begin(), cumulative.end(), i) - cumulative.begin() - 1;
    }

    template<class PerpendicularWorkspace_> 
    size_t choose_segment(size_t i, PerpendicularWorkspace_* work) const {
        if (cumulative[work->last_segment] > i) {
            if (work->last_segment && cumulative[work->last_segment - 1] <= i) {
                --(work->last_segment);
            } else {
                work->last_segment = choose_segment(i);
            }
        } else if (cumulative[work->last_segment + 1] <= i) {
            if (work->last_segment + 2 < cumulative.size() && cumulative[work->last_segment + 2] > i) {
                ++(work->last_segment);
            } else {
                work->last_segment = choose_segment(i);
            }
        }
        return work->last_segment;
    }

    template<bool WORKROW, class PerpendicularWorkspace_>
    const T* extract_one_dimension(size_t i, T* buffer, PerpendicularWorkspace_* work) const {
        size_t chosen = choose_segment(i, work);
        auto work2 = work->workspaces[chosen].get();
        if constexpr(WORKROW) {
            return mats[chosen]->row(i - cumulative[chosen], buffer, work2);
        } else {
            return mats[chosen]->column(i - cumulative[chosen], buffer, work2);
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    void assemble_along_dimension_simple(size_t i, T* buffer, ParallelWorkspace_* work) const {
        for (size_t m = 0; m < mats.size(); ++m) {
            auto curwork = work->workspaces[m].get();
            const auto& curmat = mats[m];
            if constexpr(WORKROW) {
                curmat->row_copy(i, buffer, curwork);
                buffer += curmat->ncol();
            } else {
                curmat->column_copy(i, buffer, curwork);
                buffer += curmat->nrow();
            }
        }
    }

    template<bool WORKROW, class AlongWorkspace>
    SparseRange<T, IDX> extract_one_dimension(size_t i, T* out_values, IDX* out_indices, AlongWorkspace* work, bool sorted=true) const {
        size_t chosen = choose_segment(i, work);
        auto work2 = work->workspaces[chosen].get();
        if constexpr(WORKROW) {
            return mats[chosen]->sparse_row(i - cumulative[chosen], out_values, out_indices, work2, sorted);
        } else {
            return mats[chosen]->sparse_column(i - cumulative[chosen], out_values, out_indices, work2, sorted);
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    SparseRange<T, IDX> assemble_along_dimension_simple(size_t i, T* out_values, IDX* out_indices, ParallelWorkspace_* work, bool sorted=true) const {
        size_t total = 0;
        auto originali = out_indices;
        auto originalv = out_values;

        for (size_t m = 0; m < mats.size(); ++m) {
            auto curwork = work->workspaces[m].get();
            const auto& curmat = mats[m];

            SparseRange<T, IDX> found;
            if constexpr(WORKROW) {
                found = curmat->sparse_row_copy(i, out_values, out_indices, curwork, SPARSE_COPY_BOTH, sorted);
            } else {
                found = curmat->sparse_column_copy(i, out_values, out_indices, curwork, SPARSE_COPY_BOTH, sorted);
            }

            for (size_t i = 0; i < found.number; ++i) {
                out_indices[i] += cumulative[m];
            }
            out_values += found.number;
            out_indices += found.number;
            total += found.number;
        }

        return SparseRange<T, IDX>(total, originalv, originali);
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct PerpendicularBlockWorkspace : public BlockWorkspace<WORKROW> {
        PerpendicularBlockWorkspace(size_t s, size_t l, std::vector<std::shared_ptr<BlockWorkspace<WORKROW> > > w) : details(s, l), workspaces(std::move(w)) {}

        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }

        std::vector<std::shared_ptr<BlockWorkspace<WORKROW> > > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW>
    struct ParallelBlockWorkspace : public BlockWorkspace<WORKROW> {
        ParallelBlockWorkspace(size_t s, size_t l, std::vector<std::shared_ptr<BlockWorkspace<WORKROW> > > w, std::vector<size_t> k) : details(s, l), workspaces(std::move(w)), kept(std::move(k)) {}

        std::pair<size_t, size_t> details;
        const std::pair<size_t, size_t>& block() const { return details; }

        std::vector<std::shared_ptr<BlockWorkspace<WORKROW> > > workspaces;
        std::vector<size_t> kept;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowBlockWorkspace> new_row_workspace(size_t start, size_t length) const {
        return new_workspace<true>(start, length);
    }

    std::shared_ptr<ColumnBlockWorkspace> new_column_workspace(size_t start, size_t length) const {
        return new_workspace<false>(start, length);
    }

    const T* row(size_t r, T* buffer, RowBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, buffer, static_cast<PerpendicularBlockWorkspace<true>*>(work));
        } else {
            assemble_along_dimension_complex<true>(r, buffer, static_cast<ParallelBlockWorkspace<true>*>(work));
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, ColumnBlockWorkspace* work) const {
        if constexpr(MARGIN==0) {
            assemble_along_dimension_complex<false>(c, buffer, static_cast<ParallelBlockWorkspace<false>*>(work));
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, static_cast<PerpendicularBlockWorkspace<false>*>(work));
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, static_cast<PerpendicularBlockWorkspace<true>*>(work), sorted);
        } else {
            return assemble_along_dimension_complex<true>(r, out_values, out_indices, static_cast<ParallelBlockWorkspace<true>*>(work), sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnBlockWorkspace* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return assemble_along_dimension_complex<false>(c, out_values, out_indices, static_cast<ParallelBlockWorkspace<false>*>(work), sorted);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, static_cast<PerpendicularBlockWorkspace<false>*>(work), sorted);
        }
    }

private:
    template<bool WORKROW>
    std::shared_ptr<BlockWorkspace<WORKROW> > new_workspace(size_t start, size_t length) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            std::vector<std::shared_ptr<BlockWorkspace<WORKROW> > > workspaces;
            workspaces.reserve(mats.size());

            for (const auto& x : mats) {
                if constexpr(WORKROW) {
                    workspaces.push_back(x->new_row_workspace(start, length));
                } else {
                    workspaces.push_back(x->new_column_workspace(start, length));
                }
            }

            return std::shared_ptr<BlockWorkspace<WORKROW> >(new PerpendicularBlockWorkspace<WORKROW>(start, length, std::move(workspaces)));

        } else {
            std::vector<std::shared_ptr<BlockWorkspace<WORKROW> > > workspaces;
            std::vector<size_t> kept;

            if (length) {
                workspaces.reserve(mats.size());
                kept.reserve(mats.size());

                size_t index = std::upper_bound(cumulative.begin(), cumulative.end(), start) - cumulative.begin() - 1; // finding the first matrix.
                size_t actual_start = start - cumulative[index];
                size_t end = start + length;

                for (; index < mats.size(); ++index) {
                    size_t actual_end = end;
                    bool not_final = (end > cumulative[index + 1]);
                    if (not_final) {
                        actual_end = cumulative[index + 1];
                    }
                    actual_end -= cumulative[index];

                    size_t len = actual_end - actual_start;
                    const auto& x = mats[index];
                    if constexpr(WORKROW) {
                        workspaces.push_back(x->new_row_workspace(actual_start, len));
                    } else {
                        workspaces.push_back(x->new_column_workspace(actual_start, len));
                    }
                    kept.push_back(index);

                    if (!not_final) {
                        break;
                    }
                    actual_start = 0;
                }
            }

            return std::shared_ptr<BlockWorkspace<WORKROW> >(new ParallelBlockWorkspace<WORKROW>(start, length, std::move(workspaces), std::move(kept)));
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    void assemble_along_dimension_complex(size_t i, T* buffer, ParallelWorkspace_* work) const {
        for (size_t j = 0, end = work->kept.size(); j < end; ++j) {
            auto curwork = work->workspaces[j].get();
            size_t m = work->kept[j];
            const auto& curmat = mats[m];
            if constexpr(WORKROW) {
                curmat->row_copy(i, buffer, curwork);
                buffer += curwork->length();
            } else {
                curmat->column_copy(i, buffer, curwork);
                buffer += curwork->length();
            }
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    SparseRange<T, IDX> assemble_along_dimension_complex(size_t i, T* out_values, IDX* out_indices, ParallelWorkspace_* work, bool sorted=true) const {
        size_t total = 0;
        auto originali = out_indices;
        auto originalv = out_values;

        for (size_t j = 0, end = work->kept.size(); j < end; ++j) {
            auto curwork = work->workspaces[j].get();
            size_t m = work->kept[j];
            const auto& curmat = mats[m];

            SparseRange<T, IDX> found;
            if constexpr(WORKROW) {
                found = curmat->sparse_row_copy(i, out_values, out_indices, curwork, SPARSE_COPY_BOTH, sorted);
            } else {
                found = curmat->sparse_column_copy(i, out_values, out_indices, curwork, SPARSE_COPY_BOTH, sorted);
            }

            for (size_t k = 0; k < found.number; ++k) {
                out_indices[k] += cumulative[m];
            }
            out_values += found.number;
            out_indices += found.number;
            total += found.number;
        }

        return SparseRange<T, IDX>(total, originalv, originali);
    }

public:
    /**
     * @cond
     */
    template<bool WORKROW>
    struct PerpendicularIndexWorkspace : public IndexWorkspace<IDX, WORKROW> {
        PerpendicularIndexWorkspace() = default;

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { 
            if (!workspaces.empty()) {
                return workspaces.front()->indices(); 
            } else {
                return indices_;
            }
        }

        std::vector<std::shared_ptr<IndexWorkspace<IDX, WORKROW> > > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW>
    struct ParallelIndexWorkspace : public IndexWorkspace<IDX, WORKROW> {
        ParallelIndexWorkspace(std::vector<IDX> subset) : indices_(std::move(subset)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }

        std::vector<std::shared_ptr<IndexWorkspace<IDX, WORKROW> > > workspaces;
        std::vector<size_t> kept;
    };
    /**
     * @endcond
     */

    std::shared_ptr<RowIndexWorkspace<IDX> > new_row_workspace(std::vector<IDX> subset) const {
        return new_workspace<true>(std::move(subset));
    }

    std::shared_ptr<ColumnIndexWorkspace<IDX> > new_column_workspace(std::vector<IDX> subset) const {
        return new_workspace<false>(std::move(subset));
    }

    const T* row(size_t r, T* buffer, RowIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, buffer, static_cast<PerpendicularIndexWorkspace<true>*>(work));
        } else {
            assemble_along_dimension_complex<true>(r, buffer, static_cast<ParallelIndexWorkspace<true>*>(work));
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, ColumnIndexWorkspace<IDX>* work) const {
        if constexpr(MARGIN==0) {
            assemble_along_dimension_complex<false>(c, buffer, static_cast<ParallelIndexWorkspace<false>*>(work));
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, static_cast<PerpendicularIndexWorkspace<false>*>(work));
        }
    }

    SparseRange<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, RowIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, static_cast<PerpendicularIndexWorkspace<true>*>(work), sorted);
        } else {
            return assemble_along_dimension_complex<true>(r, out_values, out_indices, static_cast<ParallelIndexWorkspace<true>*>(work), sorted);
        }
    }

    SparseRange<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, ColumnIndexWorkspace<IDX>* work, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return assemble_along_dimension_complex<false>(c, out_values, out_indices, static_cast<ParallelIndexWorkspace<false>*>(work), sorted);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, static_cast<PerpendicularIndexWorkspace<false>*>(work), sorted);
        }
    }

private:
    template<bool WORKROW>
    std::shared_ptr<IndexWorkspace<IDX, WORKROW> > new_workspace(std::vector<IDX> i) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            auto ptr = new PerpendicularIndexWorkspace<WORKROW>;
            std::shared_ptr<IndexWorkspace<IDX, WORKROW> > output(ptr);
            ptr->workspaces.reserve(mats.size());

            for (const auto& x : mats) {
                if constexpr(WORKROW) {
                    ptr->workspaces.push_back(x->new_row_workspace(i)); // deliberate copies here.
                } else {
                    ptr->workspaces.push_back(x->new_column_workspace(i));
                }
            }

            if (mats.empty()) { // make sure we can provide the indices if there are no matrices.
                ptr->indices_ = std::move(i);
            }
            return output;

        } else {
            auto ptr = new ParallelIndexWorkspace<WORKROW>(std::move(i));
            std::shared_ptr<IndexWorkspace<IDX, WORKROW> > output(ptr);
            const auto& subset = ptr->indices_;
            size_t length = subset.size();

            if (!subset.empty()) {
                auto& workspaces = ptr->workspaces;
                auto& kept = ptr->kept;

                workspaces.reserve(mats.size());
                kept.reserve(mats.size());

                size_t index = std::upper_bound(cumulative.begin(), cumulative.end(), subset[0]) - cumulative.begin() - 1; // finding the first matrix.
                size_t counter = 0;

                for (; index < mats.size(); ++index) {
                    IDX lower = cumulative[index];
                    IDX upper = cumulative[index + 1];

                    std::vector<IDX> curslice;
                    while (counter < length && subset[counter] < upper) {
                        curslice.push_back(subset[counter] - lower);
                        ++counter;
                    }

                    const auto& x = mats[index];
                    if constexpr(WORKROW) {
                        workspaces.push_back(x->new_row_workspace(std::move(curslice)));
                    } else {
                        workspaces.push_back(x->new_column_workspace(std::move(curslice)));
                    }
                    kept.push_back(index);

                    if (counter == length) {
                        break;
                    }
                }
            }

            return output;
        }
    }
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
