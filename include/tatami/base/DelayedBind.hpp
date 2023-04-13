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

    using Matrix<T, IDX>::dense_row_workspace;

    using Matrix<T, IDX>::dense_column_workspace;

    using Matrix<T, IDX>::sparse_row_workspace;

    using Matrix<T, IDX>::sparse_column_workspace;

private:
    struct DenseBase {
        DenseBase(const WorkspaceOptions& opt) {}
    };

    struct SparseBase {
        SparseBase(const WorkspaceOptions& opt) : extract_mode(opt.sparse_extract_mode) {}
        SparseExtractMode extract_mode;

        // It's always sorted anyway, no need to consider 'sorted' here.
    };

    template<class Parent>
    using ConditionalBase = typename std::conditional<Parent::sparse, SparseBase, DenseBase>::type;

    template<bool WORKROW, template<bool> class ParentWorkspace, class Parent = ParentWorkspace<WORKROW> >
    struct PerpendicularWorkspace__ : public Parent, public ConditionalBase<Parent> {
        PerpendicularWorkspace__(const WorkspaceOptions& opt, std::vector<std::shared_ptr<Parent> > w) : 
            ConditionalBase<Parent>(opt), workspaces(std::move(w)) {}

        std::vector<std::shared_ptr<Parent> > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW, template<bool> class ParentWorkspace, class Parent = ParentWorkspace<WORKROW> >
    struct ParallelWorkspace__ : public Parent, public ConditionalBase<Parent> {
        ParallelWorkspace__(const WorkspaceOptions& opt, std::vector<std::shared_ptr<Parent> > w) : 
            ConditionalBase<Parent>(opt), workspaces(std::move(w)) {}

        std::vector<std::shared_ptr<Parent> > workspaces;
    };

    template<bool WORKROW, template<bool> class ParentWorkspace>
    static auto quick_cast(ParentWorkspace<WORKROW>* wrk) {
        if constexpr((MARGIN == 0) == WORKROW) {
            return static_cast<PerpendicularWorkspace__<WORKROW, ParentWorkspace>*>(wrk);
        } else {
            return static_cast<ParallelWorkspace__<WORKROW, ParentWorkspace>*>(wrk);
        }
    }

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, DenseWorkspace>(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, DenseWorkspace>(opt);
    }

    const T* row(size_t r, T* buffer, DenseRowWorkspace* work) const {
        auto wptr = quick_cast<true, DenseWorkspace>(work);
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, buffer, wptr);
        } else {
            assemble_along_dimension_simple<true>(r, buffer, wptr);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, DenseColumnWorkspace* work) const {
        auto wptr = quick_cast<false, DenseWorkspace>(work);
        if constexpr(MARGIN==0) {
            assemble_along_dimension_simple<false>(c, buffer, wptr);
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, wptr);
        }
    }

    std::shared_ptr<SparseRowWorkspace> sparse_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, SparseWorkspace>(opt);
    }

    std::shared_ptr<SparseColumnWorkspace> sparse_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, SparseWorkspace>(opt);
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowWorkspace* work) const {
        auto wptr = quick_cast<true, SparseWorkspace>(work);
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, wptr);
        } else {
            return assemble_along_dimension_simple<true>(r, out_values, out_indices, wptr);
        }
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnWorkspace* work) const {
        auto wptr = quick_cast<false, SparseWorkspace>(work);
        if constexpr(MARGIN==0) {
            return assemble_along_dimension_simple<false>(c, out_values, out_indices, wptr);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, wptr);
        }
    }

private:
    template<bool WORKROW, template<bool> class ParentWorkspace, class Parent = ParentWorkspace<WORKROW> >
    std::shared_ptr<Parent> create_new_workspace(const WorkspaceOptions& opt) const {
        std::vector<std::shared_ptr<Parent> > workspaces;
        workspaces.reserve(mats.size());
        for (const auto& x : mats) {
            workspaces.push_back(new_workspace<WORKROW, Parent::sparse>(x.get(), opt));
        }

        if constexpr((MARGIN == 0) == WORKROW) {
            return std::shared_ptr<Parent>(new PerpendicularWorkspace__<WORKROW, ParentWorkspace>(opt, std::move(workspaces)));
        } else {
            return std::shared_ptr<Parent>(new ParallelWorkspace__<WORKROW, ParentWorkspace>(opt, std::move(workspaces)));
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
    SparseRange<T, IDX> extract_one_dimension(size_t i, T* out_values, IDX* out_indices, AlongWorkspace* work) const {
        size_t chosen = choose_segment(i, work);
        auto work2 = work->workspaces[chosen].get();
        if constexpr(WORKROW) {
            return mats[chosen]->row(i - cumulative[chosen], out_values, out_indices, work2);
        } else {
            return mats[chosen]->column(i - cumulative[chosen], out_values, out_indices, work2);
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    SparseRange<T, IDX> assemble_along_dimension_simple(size_t i, T* out_values, IDX* out_indices, ParallelWorkspace_* work) const {
        nullify_sparse_extract_pointers(work->extract_mode, out_values, out_indices);
        auto originali = out_indices;
        auto originalv = out_values;
        size_t total = 0;

        for (size_t m = 0; m < mats.size(); ++m) {
            auto curwork = work->workspaces[m].get();
            const auto& curmat = mats[m];

            SparseRange<T, IDX> found;
            if constexpr(WORKROW) {
                found = curmat->row_copy(i, out_values, out_indices, curwork);
            } else {
                found = curmat->column_copy(i, out_values, out_indices, curwork); 
            }

            if (out_indices) {
                for (size_t i = 0; i < found.number; ++i) {
                    out_indices[i] += cumulative[m];
                }
                out_indices += found.number;
            }

            if (out_values) {
                out_values += found.number;
            }

            total += found.number;
        }

        return SparseRange<T, IDX>(total, originalv, originali);
    }

private:
    template<bool WORKROW, template<bool> class ParentBlockWorkspace, class Parent = ParentBlockWorkspace<WORKROW> >
    struct PerpendicularBlockWorkspace__ : public Parent, public ConditionalBase<Parent> {
        PerpendicularBlockWorkspace__(size_t s, size_t l, const WorkspaceOptions& opt, std::vector<std::shared_ptr<Parent> > w) : 
            Parent(s, l), ConditionalBase<Parent>(opt), workspaces(std::move(w)) {}

        std::vector<std::shared_ptr<Parent> > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW, template<bool> class ParentBlockWorkspace, class Parent = ParentBlockWorkspace<WORKROW> >
    struct ParallelBlockWorkspace__ : public ParentBlockWorkspace<WORKROW>, public ConditionalBase<Parent> {
        ParallelBlockWorkspace__(size_t s, size_t l, const WorkspaceOptions& opt, std::vector<std::shared_ptr<Parent> > w, std::vector<size_t> k) : 
            Parent(s, l), ConditionalBase<Parent>(opt), workspaces(std::move(w)), kept(std::move(k)) {}

        std::vector<std::shared_ptr<Parent> > workspaces;
        std::vector<size_t> kept;
    };

    template<bool WORKROW, template<bool> class ParentBlockWorkspace>
    static auto quick_block_cast(ParentBlockWorkspace<WORKROW>* wrk) {
        if constexpr((MARGIN == 0) == WORKROW) {
            return static_cast<PerpendicularBlockWorkspace__<WORKROW, ParentBlockWorkspace>*>(wrk);
        } else {
            return static_cast<ParallelBlockWorkspace__<WORKROW, ParentBlockWorkspace>*>(wrk);
        }
    }

public:
    std::shared_ptr<DenseRowBlockWorkspace> dense_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, DenseBlockWorkspace>(start, length, opt);
    }

    std::shared_ptr<DenseColumnBlockWorkspace> dense_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, DenseBlockWorkspace>(start, length, opt);
    }

    const T* row(size_t r, T* buffer, DenseRowBlockWorkspace* work) const {
        auto wptr = quick_block_cast<true, DenseBlockWorkspace>(work);
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, buffer, wptr);
        } else {
            assemble_along_dimension_complex<true>(r, buffer, wptr);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, DenseColumnBlockWorkspace* work) const {
        auto wptr = quick_block_cast<false, DenseBlockWorkspace>(work);
        if constexpr(MARGIN==0) {
            assemble_along_dimension_complex<false>(c, buffer, wptr);
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, wptr);
        }
    }

    std::shared_ptr<SparseRowBlockWorkspace> sparse_row_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, SparseBlockWorkspace>(start, length, opt);
    }

    std::shared_ptr<SparseColumnBlockWorkspace> sparse_column_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, SparseBlockWorkspace>(start, length, opt);
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowBlockWorkspace* work) const {
        auto wptr = quick_block_cast<true, SparseBlockWorkspace>(work);
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, wptr);
        } else {
            return assemble_along_dimension_complex<true>(r, out_values, out_indices, wptr);
        }
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnBlockWorkspace* work) const {
        auto wptr = quick_block_cast<false, SparseBlockWorkspace>(work);
        if constexpr(MARGIN==0) {
            return assemble_along_dimension_complex<false>(c, out_values, out_indices, wptr);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, wptr);
        }
    }

private:
    template<bool WORKROW, template<bool> class ParentBlockWorkspace, class Parent = ParentBlockWorkspace<WORKROW> >
    std::shared_ptr<Parent> create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            std::vector<std::shared_ptr<Parent> > workspaces;
            workspaces.reserve(mats.size());
            for (const auto& x : mats) {
                workspaces.push_back(new_workspace<WORKROW, Parent::sparse>(x.get(), start, length, opt));
            }
            return std::shared_ptr<Parent>(new PerpendicularBlockWorkspace__<WORKROW, ParentBlockWorkspace>(start, length, opt, std::move(workspaces)));

        } else {
            std::vector<std::shared_ptr<Parent> > workspaces;
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
                    workspaces.push_back(new_workspace<WORKROW, Parent::sparse>(x.get(), actual_start, len, opt));
                    kept.push_back(index);

                    if (!not_final) {
                        break;
                    }
                    actual_start = 0;
                }
            }

            return std::shared_ptr<Parent>(new ParallelBlockWorkspace__<WORKROW, ParentBlockWorkspace>(start, length, opt, std::move(workspaces), std::move(kept)));
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
                buffer += curwork->length;
            } else {
                curmat->column_copy(i, buffer, curwork);
                buffer += curwork->length;
            }
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    SparseRange<T, IDX> assemble_along_dimension_complex(size_t i, T* out_values, IDX* out_indices, ParallelWorkspace_* work) const {
        size_t total = 0;
        auto originali = out_indices;
        auto originalv = out_values;

        for (size_t j = 0, end = work->kept.size(); j < end; ++j) {
            auto curwork = work->workspaces[j].get();
            size_t m = work->kept[j];
            const auto& curmat = mats[m];

            SparseRange<T, IDX> found;
            if constexpr(WORKROW) {
                found = curmat->row_copy(i, out_values, out_indices, curwork);
            } else {
                found = curmat->column_copy(i, out_values, out_indices, curwork);
            }

            if (out_indices) {
                for (size_t k = 0; k < found.number; ++k) {
                    out_indices[k] += cumulative[m];
                }
                out_indices += found.number;
            }

            if (out_values) {
                out_values += found.number;
            }

            total += found.number;
        }

        return SparseRange<T, IDX>(total, originalv, originali);
    }

private:
    template<bool WORKROW, template<typename, bool> class ParentIndexWorkspace, class Parent = ParentIndexWorkspace<IDX, WORKROW> >
    struct PerpendicularIndexWorkspace__ : public Parent, public ConditionalBase<Parent> {
        PerpendicularIndexWorkspace__(size_t l, const WorkspaceOptions& opt) : Parent(l), ConditionalBase<Parent>(opt) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { 
            if (!workspaces.empty()) {
                return workspaces.front()->indices(); 
            } else {
                return indices_;
            }
        }

        std::vector<std::shared_ptr<Parent> > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW,  template<typename, bool> class ParentIndexWorkspace, class Parent = ParentIndexWorkspace<IDX, WORKROW> >
    struct ParallelIndexWorkspace__ : public Parent, public ConditionalBase<Parent> {
        ParallelIndexWorkspace__(std::vector<IDX> subset, const WorkspaceOptions& opt) : Parent(subset.size()), ConditionalBase<Parent>(opt), indices_(std::move(subset)) {}

        std::vector<IDX> indices_;
        const std::vector<IDX>& indices() const { return indices_; }

        std::vector<std::shared_ptr<Parent> > workspaces;
        std::vector<size_t> kept;
    };

    template<bool WORKROW, template<typename, bool> class ParentIndexWorkspace>
    static auto quick_index_cast(ParentIndexWorkspace<IDX, WORKROW>* wrk) {
        if constexpr((MARGIN == 0) == WORKROW) {
            return static_cast<PerpendicularIndexWorkspace__<WORKROW, ParentIndexWorkspace>*>(wrk);
        } else {
            return static_cast<ParallelIndexWorkspace__<WORKROW, ParentIndexWorkspace>*>(wrk);
        }
    }

public:
    std::shared_ptr<DenseRowIndexWorkspace<IDX> > dense_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, DenseIndexWorkspace>(std::move(subset), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<IDX> > dense_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, DenseIndexWorkspace>(std::move(subset), opt);
    }

    const T* row(size_t r, T* buffer, DenseRowIndexWorkspace<IDX>* work) const {
        auto wptr = quick_index_cast<true, DenseIndexWorkspace>(work); 
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, buffer, wptr);
        } else {
            assemble_along_dimension_complex<true>(r, buffer, wptr);
            return buffer;
        }
    }

    const T* column(size_t c, T* buffer, DenseColumnIndexWorkspace<IDX>* work) const {
        auto wptr = quick_index_cast<false, DenseIndexWorkspace>(work); 
        if constexpr(MARGIN==0) {
            assemble_along_dimension_complex<false>(c, buffer, wptr);
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, wptr);
        }
    }

    std::shared_ptr<SparseRowIndexWorkspace<IDX> > sparse_row_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, SparseIndexWorkspace>(std::move(subset), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<IDX> > sparse_column_workspace(std::vector<IDX> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, SparseIndexWorkspace>(std::move(subset), opt);
    }

    SparseRange<T, IDX> row(size_t r, T* out_values, IDX* out_indices, SparseRowIndexWorkspace<IDX>* work) const {
        auto wptr = quick_index_cast<true, SparseIndexWorkspace>(work); 
        if constexpr(MARGIN==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, wptr);
        } else {
            return assemble_along_dimension_complex<true>(r, out_values, out_indices, wptr);
        }
    }

    SparseRange<T, IDX> column(size_t c, T* out_values, IDX* out_indices, SparseColumnIndexWorkspace<IDX>* work) const {
        auto wptr = quick_index_cast<false, SparseIndexWorkspace>(work); 
        if constexpr(MARGIN==0) {
            return assemble_along_dimension_complex<false>(c, out_values, out_indices, wptr);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, wptr);
        }
    }

private:
    template<bool WORKROW, template<typename, bool> class ParentIndexWorkspace, class Parent = ParentIndexWorkspace<IDX, WORKROW> >
    std::shared_ptr<Parent> create_new_workspace(std::vector<IDX> i, const WorkspaceOptions& opt) const {
        if constexpr((MARGIN == 0) == WORKROW) {
            auto ptr = new PerpendicularIndexWorkspace__<WORKROW, ParentIndexWorkspace>(i.size(), opt);
            std::shared_ptr<Parent> output(ptr);
            ptr->workspaces.reserve(mats.size());

            for (const auto& x : mats) {
                ptr->workspaces.push_back(new_workspace<WORKROW, Parent::sparse>(x.get(), i, opt)); // deliberate copies here.
            }

            if (mats.empty()) { // make sure we can provide the indices if there are no matrices.
                ptr->indices_ = std::move(i);
            }
            return output;

        } else {
            auto ptr = new ParallelIndexWorkspace__<WORKROW, ParentIndexWorkspace>(std::move(i), opt);
            std::shared_ptr<Parent> output(ptr);
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
                    workspaces.push_back(new_workspace<WORKROW, Parent::sparse>(x.get(), std::move(curslice), opt));
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
