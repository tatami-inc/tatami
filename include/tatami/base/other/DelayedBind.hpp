#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "../Matrix.hpp"
#include "../utils.hpp"
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
        size_t sofar = 0;
        for (size_t i = 0; i < mats.size(); ++i) {
            auto& current = mats[i];
            auto& cum = cumulative[sofar + 1];

            if constexpr(margin_==0) {
                auto nr = current->nrow();
                if (nr == 0) {
                    continue;
                }
                cum += nr;

            } else {
                auto nc = current->ncol();
                if (nc == 0) {
                    continue;
                }
                cum += nc;
            }

            cum += cumulative[sofar];
            if (sofar != i) {
                mats[sofar] = std::move(current);
            }
            ++sofar;

        }

        if (sofar != mats.size()) {
            mats.resize(sofar);
            cumulative.resize(sofar + 1);
        }

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

    /**
     * @param ps Pointers to the matrices to be combined.
     */
    DelayedBind(const std::vector<std::shared_ptr<Matrix<Value_, Index_> > >& ps) : DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >(ps.begin(), ps.end())) {}

private:
    std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > mats;
    std::vector<Index_> cumulative;

    bool sparse_;
    std::pair<double, double> dimension_preference_;

private:
    Index_ internal_nrow() const {
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

    Index_ internal_ncol() const {
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

public:
    Index_ nrow() const {
        return internal_nrow();
    }

    Index_ ncol() const {
        return internal_ncol();
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
    template<DimensionSelectionType selection_, bool sparse_>
    struct ParallelExtractor : public Extractor<selection_, sparse_, Value_, Index_> {
        ParallelExtractor(const DelayedBind* p, const Options<Index_>& opt) : parent(p) {
            constexpr bool accrow_ = margin_ != 0;
            workspaces.reserve(p->mats.size());

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? p->internal_ncol() : p->internal_nrow());
                for (const auto& m : parent->mats) {
                    // TODO: manage opt.selection, manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), opt));
                }
            }
        }

        ParallelExtractor(const DelayedBind* p, const Options<Index_>& opt, Index_ bs, Index_ bl) : parent(p) {
            constexpr bool accrow_ = margin_ != 0;
            size_t nmats = parent->mats.size();
            workspaces.reserve(nmats);

            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;

                const auto& cumulative = this->parent->cumulative;
                size_t index = std::upper_bound(cumulative.begin(), cumulative.end(), this->block_start) - cumulative.begin() - 1; // finding the first matrix.
                Index_ actual_start = this->block_start - cumulative[index];
                Index_ end = this->block_start + this->block_length;
                if constexpr(sparse_) {
                    kept.reserve(nmats);
                }

                for (; index < nmats; ++index) {
                    bool not_final = (end > cumulative[index + 1]);
                    Index_ actual_end = (not_final ? cumulative[index + 1] : end) - cumulative[index];

                    // TODO: manage opt.selection, manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(this->parent->mats[index].get(), actual_start, actual_end - actual_start, opt));

                    if constexpr(sparse_) {
                        kept.push_back(index);
                    }

                    if (!not_final) {
                        break;
                    }
                    actual_start = 0;
                }
            }
        }

        ParallelExtractor(const DelayedBind* p, const Options<Index_>& opt, const Index_* is, size_t il) : parent(p) {
            constexpr bool accrow_ = margin_ != 0;
            size_t nmats = parent->mats.size();
            workspaces.reserve(nmats);

            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = il;

                if (il) {
                    const auto& cumulative = this->parent->cumulative;
                    size_t index = std::upper_bound(cumulative.begin(), cumulative.end(), is[0]) - cumulative.begin() - 1; // finding the first matrix.
                    if constexpr(sparse_) {
                        kept.reserve(nmats);
                    }

                    size_t counter = 0;
                    indices.reserve(il);
                    for (; index < nmats; ++index) {
                        Index_ lower = cumulative[index];
                        Index_ upper = cumulative[index + 1];

                        size_t old_counter = counter;
                        while (counter < il && is[counter] < upper) {
                            indices.push_back(is[counter] - lower);
                            ++counter;
                        }

                        if (old_counter != counter) {
                            auto indices_start = indices.data() + old_counter;

                            // TODO: manage opt.selection, manage opt.access.sequencer!
                            workspaces.push_back(new_extractor<accrow_, sparse_>(this->parent->mats[index].get(), indices_start, counter - old_counter, opt));
                            if constexpr(sparse_) {
                                kept.push_back(index);
                            }

                            // Actually setting the indices properly.
                            std::copy(is + old_counter, is + counter, indices_start);
                        }

                        if (counter == il) {
                            break;
                        }
                    }
                }
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

    protected:
        const DelayedBind* parent;
        std::vector<std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > > workspaces;
        typename std::conditional<sparse_ && selection_ != DimensionSelectionType::FULL, std::vector<size_t>, bool>::type kept;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
    };

private:
    template<DimensionSelectionType selection_>
    struct DenseParallelExtractor : public ParallelExtractor<selection_, false> {
        template<typename ... Args_>
        DenseParallelExtractor(const DelayedBind* p, const Options<Index_>& opt, Args_... args) : ParallelExtractor<selection_, false>(p, opt, args...) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto copy = buffer;
            for (auto& w : this->workspaces) {
                w->fetch_copy(i, copy);
                copy += extracted_length<selection_, Index_>(*w);
            }
            return buffer;
        }
    };

    template<DimensionSelectionType selection_>
    struct SparseParallelExtractor : public ParallelExtractor<selection_, true> {
        template<typename ... Args_>
        SparseParallelExtractor(const DelayedBind* p, const Options<Index_>& opt, Args_... args) : 
            ParallelExtractor<selection_, true>(p, opt, args...),
            needs_value(opt.sparse.extract_value), 
            needs_index(opt.sparse.extract_index)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto vcopy = vbuffer;
            auto icopy = ibuffer;
            Index_ total = 0;
            size_t counter = 0;

            for (auto& w : this->workspaces) {
                auto raw = w->fetch_copy(i, vcopy, icopy);
                total += raw.number;

                if (needs_value) {
                    vcopy += raw.number;
                }

                if (needs_index) {
                    Index_ adjustment;
                    if constexpr(selection_ == DimensionSelectionType::FULL) {
                        adjustment = this->parent->cumulative[counter];
                    } else {
                        adjustment = this->parent->cumulative[this->kept[counter]];
                    }

                    for (Index_ j = 0; j < raw.number; ++j) {
                        icopy[j] += adjustment;
                    }

                    icopy += raw.number;
                }

                ++counter;
            }

            if (!needs_value) {
                vbuffer = NULL;
            }
            if (!needs_index) {
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
    template<DimensionSelectionType selection_, bool sparse_>
    struct PerpendicularExtractor : public Extractor<selection_, sparse_, Value_, Index_> {
        PerpendicularExtractor(const DelayedBind* p, const Options<Index_>& opt) : parent(p) {
            constexpr bool accrow_ = margin_ == 0;
            workspaces.reserve(parent->mats.size());

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? p->internal_ncol() : p->internal_nrow());
                for (const auto& m : parent->mats) {
                    // TODO: manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), opt));
                }
            }
        }

        PerpendicularExtractor(const DelayedBind* p, const Options<Index_>& opt, Index_ bs, Index_ bl) : parent(p) {
            constexpr bool accrow_ = margin_ == 0;
            workspaces.reserve(p->mats.size());

            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
                for (const auto& m : parent->mats) {
                    // TODO: manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), bs, bl, opt));
                }
            }
        }

        PerpendicularExtractor(const DelayedBind* p, const Options<Index_>& opt, const Index_* is, size_t il) : parent(p) {
            constexpr bool accrow_ = margin_ == 0;
            workspaces.reserve(p->mats.size());

            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = il;
                for (const auto& m : parent->mats) {
                    // TODO: manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), is, il, opt));
                }

                if (workspaces.empty()) {
                    indices = std::vector<Index_>(is, is + il);
                }
            }
        }

        const Index_* index_start() const {
            if (workspaces.empty()) {
                return indices.data();
            } else {
                return workspaces.front()->index_start();
            }
        }

    protected:
        const DelayedBind* parent;
        std::vector<std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > > workspaces;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
        size_t last_segment = 0;

        size_t choose_segment(size_t i) {
            if (parent->cumulative[last_segment] > i) {
                if (last_segment && parent->cumulative[last_segment - 1] <= i) {
                    --last_segment;
                } else {
                    last_segment = choose_segment_raw(i);
                }
            } else if (parent->cumulative[last_segment + 1] <= i) {
                if (last_segment + 2 < parent->cumulative.size() && parent->cumulative[last_segment + 2] > i) {
                    ++last_segment;
                } else {
                    last_segment = choose_segment_raw(i);
                }
            }
            return last_segment;
        }

    private:
        size_t choose_segment_raw(size_t i) const {
            return std::upper_bound(parent->cumulative.begin(), parent->cumulative.end(), i) - parent->cumulative.begin() - 1;
        }
    };

private:
    template<DimensionSelectionType selection_>
    struct DensePerpendicularExtractor : public PerpendicularExtractor<selection_, false> {
        template<typename ... Args_>
        DensePerpendicularExtractor(const DelayedBind* p, const Options<Index_>& opt, Args_... args) : PerpendicularExtractor<selection_, false>(p, opt, args...) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            size_t chosen = this->choose_segment(i);
            return this->workspaces[chosen]->fetch(i - this->parent->cumulative[chosen], buffer);
        }
    };

    template<DimensionSelectionType selection_>
    struct SparsePerpendicularExtractor : public PerpendicularExtractor<selection_, true> {
        template<typename ... Args_>
        SparsePerpendicularExtractor(const DelayedBind* p, const Options<Index_>& opt, Args_... args) : PerpendicularExtractor<selection_, true>(p, opt, args...) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            size_t chosen = this->choose_segment(i);
            return this->workspaces[chosen]->fetch(i - this->parent->cumulative[chosen], vbuffer, ibuffer);
        }
    };

    /************************************
     ********** Public methods **********
     ************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options<Index_>& opt, Args... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(sparse_) {
            if constexpr(accrow_ == (margin_ == 0)) {
                output.reset(new SparsePerpendicularExtractor<selection_>(this, opt, args...));
            } else {
                output.reset(new SparseParallelExtractor<selection_>(this, opt, args...));
            }
        } else {
            if constexpr(accrow_ == (margin_ == 0)) {
                output.reset(new DensePerpendicularExtractor<selection_>(this, opt, args...));
            } else {
                output.reset(new DenseParallelExtractor<selection_>(this, opt, args...));
            }
        }

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options<Index_>& opt) const {
        return populate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const {
        return populate<true, DimensionSelectionType::INDEX, false>(opt, index_start, index_length);
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options<Index_>& opt) const {
        return populate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const {
        return populate<false, DimensionSelectionType::INDEX, false>(opt, index_start, index_length);
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options<Index_>& opt) const {
        return populate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const {
        return populate<true, DimensionSelectionType::INDEX, true>(opt, index_start, index_length);
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options<Index_>& opt) const {
        return populate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(const Index_* index_start, size_t index_length, const Options<Index_>& opt) const {
        return populate<false, DimensionSelectionType::INDEX, true>(opt, index_start, index_length);
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam margin_ Dimension along which the combining is to occur.
 * If 0, matrices are combined along the rows; if 1, matrices are combined to the columns.
 * @tparam Matrix_ A realized `Matrix` class, possibly const. 
 *
 * @param ps Pointers to `Matrix` objects.
 *
 * @return A pointer to a `DelayedBind` instance.
 */
template<int margin_, class Matrix_>
std::shared_ptr<Matrix_> make_DelayedBind(std::vector<std::shared_ptr<Matrix_> > ps) {
    return std::shared_ptr<Matrix_>(new DelayedBind<margin_, typename Matrix_::value_type, typename Matrix_::index_type>(std::move(ps)));
}

}

#endif
