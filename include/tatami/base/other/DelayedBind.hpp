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
 * @tparam margin_ Dimension along which the combining is to occur.
 * If 0, the matrices are combined along the rows; if 1, the combining is applied along the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Data_ Type of index value.
 */
template<int margin_, typename Value_, typename Index_>
class DelayedBind : public Matrix<Value_, Index_> {
public:
    /**
     * @param ps Pointers to the matrices to be combined.
     */
    DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > ps) : mats(std::move(ps)), cumulative(mats.size()+1) {
        for (size_t i = 0; i < mats.size(); ++i) {
            if constexpr(margin_==0) {
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
    DelayedBind(const std::vector<std::shared_ptr<Matrix<Value_, Index_> > >& ps) : mats(ps.begin(), ps.end()), cumulative(mats.size()+1) {
        for (size_t i = 0; i < mats.size(); ++i) {
            if constexpr(margin_==0) {
                cumulative[i+1] = mats[i]->nrow();
            } else {
                cumulative[i+1] = mats[i]->ncol();
            }
            cumulative[i+1] += cumulative[i];
        }

        set_constants();
    }

private:
    std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > mats;
    std::vector<Index_> cumulative;

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
    Index_ nrow() const {
        if constexpr(margin_==0) {
            return cumulative.back();
        } else {
            if (mats.empty()) {
                return 0;
            } else {
                return mats.front()->nrow();
            }
        }
    }

    Index_ ncol() const {
        if constexpr(margin_==0) {
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

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

    /****************************************
     ********** Parallel extractor **********
     ****************************************/
private:
    template<DimensionSelectionType selection_, bool use_start_, bool sparse_>
    struct ParallelExtractor : public StandardExtractor<selection_, use_start_, sparse_, Value_, Index_> {
        ParallelExtractor(const DelayedBind* p, IterationOptions<Index_>& iopt, ExtractionOptions<Index_>& eopt) : 
            StandardExtractor<selection_, use_start_, sparse_, Value_, Index_>(eopt),
            parent(p)
        {
            constexpr bool accrow_ = margin_ != 0;
            workspaces.reserve(p->mats.size());

            // We need to create copies of 'eopt', so we flush the indices to keep it cheap.
            // This information should already have been moved to the StandardExtractor anyway.
            eopt.selection.indices.clear(); 

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                for (const auto& m : mats) {
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), iopt, eopt));
                }

            } else if constexpr(selection_ == DimensionSelection::BLOCK) {
                size_t index = std::upper_bound(cumulative.begin(), cumulative.end(), this->extracted_block) - cumulative.begin() - 1; // finding the first matrix.
                Index_ actual_start = this->extracted_block - cumulative[index];
                Index_ end = this->extracted_block + this->extracted_length;
                if constexpr(sparse_) {
                    kept.reserve(mats.size());
                }

                for (; index < mats.size(); ++index) {
                    bool not_final = (end > cumulative[index + 1]);

                    auto ecopy = eopt;
                    if (not_final && actual_start == 0) {
                        // Switch to full extraction, which may be more optimal for the underlying matrix.
                        ecopy.selection.type = DimensionSelectionType::FULL;
                        workspaces.push_back(new_extractor<accrow_, sparse_>(mats[index].get(), iopt, std::move(ecopy)));
                    } else {
                        Index_ actual_end = end - cumulative[index];
                        ecopy.block_start = actual_start;
                        ecopy.block_end = actual_end - actual_start;
                        workspaces.push_back(new_extractor<accrow_, sparse_>(mats[index].get(), iopt, std::move(ecopy)));
                    }

                    if constexpr(sparse_) {
                        kept.push_back(index);
                    }

                    if (!not_final) {
                        break;
                    }
                    actual_start = 0;
                }

            } else {
                Index_ length = this->extracted_length;

                if (length) {
                    auto subset = this->quick_extract_index();
                    size_t index = std::upper_bound(cumulative.begin(), cumulative.end(), subset[0]) - cumulative.begin() - 1; // finding the first matrix.
                    if constexpr(sparse_) {
                        kept.reserve(mats.size());
                    }

                    Index_ counter = 0;
                    for (; index < mats.size(); ++index) {
                        Index_ lower = cumulative[index];
                        Index_ upper = cumulative[index + 1];

                        auto ecopy = eopt;
                        ecopy.index_start = NULL;
                        auto& curslice = ecopy.selection.indices;
                        while (counter < length && subset[counter] < upper) {
                            curslice.push_back(subset[counter] - lower);
                            ++counter;
                        }

                        if (curslice.size()) {
                            const auto& x = mats[index];
                            workspaces.push_back(new_extractor<accrow_ == (margin_ == 0), sparse_>(mats[index].get(), iopt, std::move(ecopy)));
                            if constexpr(sparse_) {
                                kept.push_back(index);
                            }
                        }

                        if (counter == length) {
                            break;
                        }
                    }
                }
            }
        }

    protected:
        const DelayedBind* parent;
        std::vector<std::unique_ptr<Extractor<sparse_, Value_, Index_> > > workspaces;
        typename std::conditional<sparse && selection_ != DimensionSelectionType::FULL, std::vector<size_t>, bool> kept;
    };

private:
    template<DimensionSelectionType selection_, bool use_start_>
    struct DenseParallelExtractor : public ParallelExtractor<selection_, use_start_, false> {
        DenseParallelExtractor(const DelayedBind* p, IterationOptions<Index_>& iopt, ExtractionOptions<Index_>& eopt) :
            ParallelExtractor<selection_, use_start_, false>(p, iopt, eopt)
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto copy = buffer;
            for (auto& w : workspaces) {
                w->fetch_copy(i, copy);
                copy += w->extracted_length;
            }
            return buffer;
        }
    };

    template<DimensionSelectionType selection_, bool use_start_>
    struct SparseParallelExtractor : public ParallelExtractor<selection_, use_start_, true> {
        SparseParallelExtractor(const DelayedBind* p, IterationOptions<Index_>& iopt, ExtractionOptions<Index_>& eopt) :
            ParallelExtractor<selection_, use_start_, true>(p, iopt, eopt), 
            needs_value(eopt.sparse_extract_value), 
            needs_index(eopt.sparse_extract_index)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto vcopy = vbuffer;
            auto icopy = ibuffer;
            Index_ total = 0;
            size_t counter = 0;

            for (auto& w : workspaces) {
                auto raw = w->fetch_copy(i, vcopy, icopy);
                total += raw.number;

                if (needs_value) {
                    vcopy += raw.number;
                }

                if (needs_index) {
                    Index_ adjustment;
                    if constexpr(selection_ == DimensionSelectionType::FULL) {
                        adjustment = cumulative[counter];
                    } else {
                        adjustment = cumulative[keep][counter];
                    }

                    for (Index_ j = 0; j < raw.number; ++j) {
                        icopy[j] += adjustment;
                    }

                    icopy += raw.number;
                }
            }

            if (needs_value) {
                vbuffer = NULL;
            }
            if (needs_index) {
                ibuffer = NULL;
            }

            return SparseRange<Value_, Index_>(total, vbuffer, ibuffer);
        }

    protected:
        bool needs_value;
        bool needs_index;
    };

    /*********************************************
     ********** Perpendicular extractor **********
     *********************************************/
private:
    template<DimensionSelectionType selection_, bool use_start_, bool sparse_>
    struct PerpendicularExtractor : public StandardExtractor<selection_, use_start_, sparse_, Value_, Index_> {
        PerpendicularExtractor(const DelayedBind* p, IterationOptions<Index_>& iopt, ExtractionOptions<Index_>& eopt) : 
            StandardExtractor<selection_, use_start_, sparse_, Value_, Index_>(eopt),
            parent(p)
    {
        constexpr bool accrow_ = margin_ == 0;
        workspaces.reserve(p->mats.size());

        // We need to create copies of 'eopt', so we flush the indices to keep it cheap.
        // This information should already have been moved to the StandardExtractor anyway.
        eopt.selection.indices.clear(); 

        if constexpr(selection_ != DimensionSelectionType::INDEX) {
            for (const auto& m : mats) {
                workspaces.push_back(new_extractor<accrow_, sparse_, Value_, Index_>(m.get(), iopt, eopt));
            }
        } else {
            // Pass along pointers instead; lifetime is guaranteed?
            auto ptr = this->quick_extracted_index();
            eopt.selection.indices.insert(eopt.selection.indices.end(), ptr, ptr + this->extracted_length);
            for (const auto& m : mats) {
                workspaces.push_back(new_extractor<accrow_, sparse_, Value_, Index_>(m.get(), iopt, eopt));
            }
        }
   }

    protected:
        const DelayedBind* parent;
        std::vector<std::unique_ptr<Extractor<sparse_, Value_, Index_> > > workspaces;
    };

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

public:
    std::shared_ptr<DenseRowWorkspace> dense_row_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<true, DenseWorkspace>(opt);
    }

    std::shared_ptr<DenseColumnWorkspace> dense_column_workspace(const WorkspaceOptions& opt) const {
        return create_new_workspace<false, DenseWorkspace>(opt);
    }

    const Value_* row(size_t r, Value_* buffer, DenseRowWorkspace* work) const {
        auto wptr = quick_cast<true, DenseWorkspace>(work);
        if constexpr(margin_==0) {
            return extract_one_dimension<true>(r, buffer, wptr);
        } else {
            assemble_along_dimension_simple<true>(r, buffer, wptr);
            return buffer;
        }
    }

    const Value_* column(size_t c, Value_* buffer, DenseColumnWorkspace* work) const {
        auto wptr = quick_cast<false, DenseWorkspace>(work);
        if constexpr(margin_==0) {
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

    SparseRange<Value_, Index_> row(size_t r, Value_* out_values, Index_* out_indices, SparseRowWorkspace* work) const {
        auto wptr = quick_cast<true, SparseWorkspace>(work);
        if constexpr(margin_==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, wptr);
        } else {
            return assemble_along_dimension_simple<true>(r, out_values, out_indices, wptr);
        }
    }

    SparseRange<Value_, Index_> column(size_t c, Value_* out_values, Index_* out_indices, SparseColumnWorkspace* work) const {
        auto wptr = quick_cast<false, SparseWorkspace>(work);
        if constexpr(margin_==0) {
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

        if constexpr((margin_ == 0) == WORKROW) {
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
    const Value_* extract_one_dimension(size_t i, Value_* buffer, PerpendicularWorkspace_* work) const {
        size_t chosen = choose_segment(i, work);
        auto work2 = work->workspaces[chosen].get();
        if constexpr(WORKROW) {
            return mats[chosen]->row(i - cumulative[chosen], buffer, work2);
        } else {
            return mats[chosen]->column(i - cumulative[chosen], buffer, work2);
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    void assemble_along_dimension_simple(size_t i, Value_* buffer, ParallelWorkspace_* work) const {
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
    SparseRange<Value_, Index_> extract_one_dimension(size_t i, Value_* out_values, Index_* out_indices, AlongWorkspace* work) const {
        size_t chosen = choose_segment(i, work);
        auto work2 = work->workspaces[chosen].get();
        if constexpr(WORKROW) {
            return mats[chosen]->row(i - cumulative[chosen], out_values, out_indices, work2);
        } else {
            return mats[chosen]->column(i - cumulative[chosen], out_values, out_indices, work2);
        }
    }

    template<bool WORKROW, class ParallelWorkspace_>
    SparseRange<Value_, Index_> assemble_along_dimension_simple(size_t i, Value_* out_values, Index_* out_indices, ParallelWorkspace_* work) const {
        nullify_sparse_extract_pointers(work->extract_mode, out_values, out_indices);
        auto originali = out_indices;
        auto originalv = out_values;
        size_t total = 0;

        for (size_t m = 0; m < mats.size(); ++m) {
            auto curwork = work->workspaces[m].get();
            const auto& curmat = mats[m];

            SparseRange<Value_, Index_> found;
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

        return SparseRange<Value_, Index_>(total, originalv, originali);
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
        if constexpr((margin_ == 0) == WORKROW) {
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

    const Value_* row(size_t r, Value_* buffer, DenseRowBlockWorkspace* work) const {
        auto wptr = quick_block_cast<true, DenseBlockWorkspace>(work);
        if constexpr(margin_==0) {
            return extract_one_dimension<true>(r, buffer, wptr);
        } else {
            assemble_along_dimension_complex<true>(r, buffer, wptr);
            return buffer;
        }
    }

    const Value_* column(size_t c, Value_* buffer, DenseColumnBlockWorkspace* work) const {
        auto wptr = quick_block_cast<false, DenseBlockWorkspace>(work);
        if constexpr(margin_==0) {
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

    SparseRange<Value_, Index_> row(size_t r, Value_* out_values, Index_* out_indices, SparseRowBlockWorkspace* work) const {
        auto wptr = quick_block_cast<true, SparseBlockWorkspace>(work);
        if constexpr(margin_==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, wptr);
        } else {
            return assemble_along_dimension_complex<true>(r, out_values, out_indices, wptr);
        }
    }

    SparseRange<Value_, Index_> column(size_t c, Value_* out_values, Index_* out_indices, SparseColumnBlockWorkspace* work) const {
        auto wptr = quick_block_cast<false, SparseBlockWorkspace>(work);
        if constexpr(margin_==0) {
            return assemble_along_dimension_complex<false>(c, out_values, out_indices, wptr);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, wptr);
        }
    }

private:
    template<bool WORKROW, template<bool> class ParentBlockWorkspace, class Parent = ParentBlockWorkspace<WORKROW> >
    std::shared_ptr<Parent> create_new_workspace(size_t start, size_t length, const WorkspaceOptions& opt) const {
        if constexpr((margin_ == 0) == WORKROW) {
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
    void assemble_along_dimension_complex(size_t i, Value_* buffer, ParallelWorkspace_* work) const {
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
    SparseRange<Value_, Index_> assemble_along_dimension_complex(size_t i, Value_* out_values, Index_* out_indices, ParallelWorkspace_* work) const {
        size_t total = 0;
        auto originali = out_indices;
        auto originalv = out_values;

        for (size_t j = 0, end = work->kept.size(); j < end; ++j) {
            auto curwork = work->workspaces[j].get();
            size_t m = work->kept[j];
            const auto& curmat = mats[m];

            SparseRange<Value_, Index_> found;
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

        return SparseRange<Value_, Index_>(total, originalv, originali);
    }

private:
    template<bool WORKROW, template<typename, bool> class ParentIndexWorkspace, class Parent = ParentIndexWorkspace<Index_, WORKROW> >
    struct PerpendicularIndexWorkspace__ : public Parent, public ConditionalBase<Parent> {
        PerpendicularIndexWorkspace__(size_t l, const WorkspaceOptions& opt) : Parent(l), ConditionalBase<Parent>(opt) {}

        std::vector<Index_> indices_;
        const std::vector<Index_>& indices() const { 
            if (!workspaces.empty()) {
                return workspaces.front()->indices(); 
            } else {
                return indices_;
            }
        }

        std::vector<std::shared_ptr<Parent> > workspaces;
        size_t last_segment = 0;
    };

    template<bool WORKROW,  template<typename, bool> class ParentIndexWorkspace, class Parent = ParentIndexWorkspace<Index_, WORKROW> >
    struct ParallelIndexWorkspace__ : public Parent, public ConditionalBase<Parent> {
        ParallelIndexWorkspace__(std::vector<Index_> subset, const WorkspaceOptions& opt) : Parent(subset.size()), ConditionalBase<Parent>(opt), indices_(std::move(subset)) {}

        std::vector<Index_> indices_;
        const std::vector<Index_>& indices() const { return indices_; }

        std::vector<std::shared_ptr<Parent> > workspaces;
        std::vector<size_t> kept;
    };

    template<bool WORKROW, template<typename, bool> class ParentIndexWorkspace>
    static auto quick_index_cast(ParentIndexWorkspace<Index_, WORKROW>* wrk) {
        if constexpr((margin_ == 0) == WORKROW) {
            return static_cast<PerpendicularIndexWorkspace__<WORKROW, ParentIndexWorkspace>*>(wrk);
        } else {
            return static_cast<ParallelIndexWorkspace__<WORKROW, ParentIndexWorkspace>*>(wrk);
        }
    }

public:
    std::shared_ptr<DenseRowIndexWorkspace<Index_> > dense_row_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, DenseIndexWorkspace>(std::move(subset), opt);
    }

    std::shared_ptr<DenseColumnIndexWorkspace<Index_> > dense_column_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, DenseIndexWorkspace>(std::move(subset), opt);
    }

    const Value_* row(size_t r, Value_* buffer, DenseRowIndexWorkspace<Index_>* work) const {
        auto wptr = quick_index_cast<true, DenseIndexWorkspace>(work); 
        if constexpr(margin_==0) {
            return extract_one_dimension<true>(r, buffer, wptr);
        } else {
            assemble_along_dimension_complex<true>(r, buffer, wptr);
            return buffer;
        }
    }

    const Value_* column(size_t c, Value_* buffer, DenseColumnIndexWorkspace<Index_>* work) const {
        auto wptr = quick_index_cast<false, DenseIndexWorkspace>(work); 
        if constexpr(margin_==0) {
            assemble_along_dimension_complex<false>(c, buffer, wptr);
            return buffer;
        } else {
            return extract_one_dimension<false>(c, buffer, wptr);
        }
    }

    std::shared_ptr<SparseRowIndexWorkspace<Index_> > sparse_row_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<true, SparseIndexWorkspace>(std::move(subset), opt);
    }

    std::shared_ptr<SparseColumnIndexWorkspace<Index_> > sparse_column_workspace(std::vector<Index_> subset, const WorkspaceOptions& opt) const {
        return create_new_workspace<false, SparseIndexWorkspace>(std::move(subset), opt);
    }

    SparseRange<Value_, Index_> row(size_t r, Value_* out_values, Index_* out_indices, SparseRowIndexWorkspace<Index_>* work) const {
        auto wptr = quick_index_cast<true, SparseIndexWorkspace>(work); 
        if constexpr(margin_==0) {
            return extract_one_dimension<true>(r, out_values, out_indices, wptr);
        } else {
            return assemble_along_dimension_complex<true>(r, out_values, out_indices, wptr);
        }
    }

    SparseRange<Value_, Index_> column(size_t c, Value_* out_values, Index_* out_indices, SparseColumnIndexWorkspace<Index_>* work) const {
        auto wptr = quick_index_cast<false, SparseIndexWorkspace>(work); 
        if constexpr(margin_==0) {
            return assemble_along_dimension_complex<false>(c, out_values, out_indices, wptr);
        } else {
            return extract_one_dimension<false>(c, out_values, out_indices, wptr);
        }
    }

private:
    template<bool WORKROW, template<typename, bool> class ParentIndexWorkspace, class Parent = ParentIndexWorkspace<Index_, WORKROW> >
    std::shared_ptr<Parent> create_new_workspace(std::vector<Index_> i, const WorkspaceOptions& opt) const {
        if constexpr((margin_ == 0) == WORKROW) {
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
                    Index_ lower = cumulative[index];
                    Index_ upper = cumulative[index + 1];

                    std::vector<Index_> curslice;
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
 * @tparam margin_ Dimension along which the combining is to occur.
 * If 0, matrices are combined along the rows; if 1, matrices are combined to the columns.
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 *
 * @param ps Pointers to `Matrix` objects.
 *
 * @return A pointer to a `DelayedBind` instance.
 */
template<int margin_, class MAT>
std::shared_ptr<MAT> make_DelayedBind(std::vector<std::shared_ptr<MAT> > ps) {
    return std::shared_ptr<MAT>(
        new DelayedBind<margin_, typename MAT::data_type, typename MAT::index_type>(std::move(ps))
    );
}

}

#endif
