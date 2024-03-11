#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"

#include <algorithm>
#include <memory>
#include <array>
#include <deque>

/**
 * @file DelayedBind.hpp
 *
 * @brief Delayed combining of multiple `tatami::Matrix` objects.
 *
 * This is equivalent to the `DelayedAbind` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedBind_internal {

template<typename Value_, typename Index_>
struct ParallelFullDense {
    ParallelFullDense(const std::vector<Index_>& cum) : cumulative(cum) {}

protected:
    const std::vector<Index_>& cumulative;

    template<class Store_>
    void fetch_base(Value_* buffer, Store_ store) {
        for (size_t x = 0, end = cumulative.size() - 1; x < end; ++x) {
            Index_ last = cumulative[x];
            auto copy = buffer + last;
            auto ptr = store(x, copy);
            copy_n(ptr, cumulative[x + 1] - last, copy);
        }
    }
};

template<typename Value_, typename Index_>
struct MyopicParallelFullDense : public MyopicDenseExtractor<Value_, Index_>, public ParallelFullDense<Value_, Index_> {
    MyopicParallelFullBase(const std::vector<Index_>& cum, std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > in) : 
        ParallelFullDense<Value_, Index_>(cum), internal(std::move(in)) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        this->fetch_base(buffer, [&](size_t x, Value_* copy) -> Value_* { return internal[x]->fetch(i, copy); });
        return buffer;
    }

    Index_ number() const {
        return this->cumulative.back();
    }

private:
    std::vector<std::unique_ptr<MyopicDenseExtractor_<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_>
struct OracularParallelFullDense : public OracularDenseExtractor<Value_, Index_>, public ParallelFullDense<Value_, Index_> {
    OracularParallelFullBase(const std::vector<Index_>& cum, std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > in, std::shared_ptr<Oracle<Index_> > ora) : 
        ParallelFullDense<Value_, Index_>(cum), internal(std::move(in)) 
    {
        if (internal.empty()) {
            oracle = std::move(ora);
        }
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        if (internal.empty()) {
            i = oracle->get(used);
            ++used;
        } else {
            this->fetch_base(buffer, [&](size_t x, Value_* copy) -> Value_* { return internal[x]->fetch(i, copy); });
        }
        return buffer;
    }

    Index_ number() const {
        return this->cumulative.back();
    }

private:
    std::vector<std::unique_ptr<OracularDenseExtractor_<Value_, Index_> > > internal;
    std::shared_ptr<Oracle<Index_> > oracle;
    size_t used = 0;
};

template<typename Value_, typename Index_>
struct ParallelFullSparse {
    ParallelFullSparse(const std::vector<Index_>& cum, const Options& opt) : 
        cumulative(cum), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {}

protected:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;

    template<class Store_>
    SparseRange<Value_, Index_> fetch_base(Value_* vbuffer, Index_* ibuffer, Store_ store) {
        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        Index_ accumulated = 0;

        for (size_t x = 0, end = cumulative.size() - 1; x < end; ++x) {
            auto range = store(x, vcopy, icopy);
            accumulated += range.number;

            if (needs_value) {
                copy_n(range.value, range.number, vcopy);
                vcopy += range.number;
            }
            if (needs_index) {
                auto offset = cumulative[x];
                for (Index_ y = 0; y < range.number; ++y) {
                    icopy[y] = range.index[y] + offset;
                }
                icopy += range.number;
            }
        }

        return SparseRange<Value_, Index_>(accumulated, (needs_value ? vbuffer : NULL), (needs_index ? ibuffer : NULL));
    }
};

template<typename Value_, typename Index_>
struct MyopicParallelFullSparse : public MyopicSparseExtractor<Value_, Index_>, public ParallelFullSparse<Value_, Index_> {
    MyopicParallelFullBase(const std::vector<Index_>& cum, std::vector<std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > > in) : 
        ParallelFullSparse<Value_, Index_>(cum), internal(std::move(in)) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return this->fetch_base(vbuffer, ibuffer, [&](size_t x, Value_* vcopy, Index_* icopy) -> Value_* { return internal[x]->fetch(i, vcopy, icopy); });
    }

    Index_ number() const {
        return this->cumulative.back();
    }

private:
    std::vector<std::unique_ptr<MyopicSparseExtractor_<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_>
struct OracularParallelFullSparse : public OracularSparseExtractor<Value_, Index_>, public ParallelFullSparse<Value_, Index_> {
    OracularParallelFullBase(const std::vector<Index_>& cum, std::vector<std::unique_ptr<OracularSparseExtractor<Value_, Index_> > > in, std::shared_ptr<Oracle<Index_> > ora) : 
        ParallelFullSparse<Value_, Index_>(cum), internal(std::move(in)) 
    {
        if (internal.empty()) {
            oracle = std::move(ora);
        }
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        if (internal.empty()) {
            i = oracle->get(used);
            ++used;
        } 
        return this->fetch_base(vbuffer, ibuffer, [&](size_t x, Value_* vcopy, Index_* icopy) -> Value_* { return internal[x]->fetch(i, vcopy, icopy); });
    }

    Index_ number() const {
        return this->cumulative.back();
    }

private:
    std::vector<std::unique_ptr<OracularSparseExtractor_<Value_, Index_> > > internal;
    std::shared_ptr<Oracle<Index_> > oracle;
    size_t used = 0;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed combining of a matrix.
 *
 * Implements delayed combining by rows or columns of a matrix.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the combining is to occur.
 * If 0, the matrices are combined along the rows; if 1, the combining is applied along the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
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

        double denom = 0;
        for (const auto& x : mats) {
            double total = x->nrow() * x->ncol();
            denom += total;
            sparse_prop += total * x->sparse_proportion();
            row_prop += total * x->prefer_rows_proportion();
        }
        if (denom) {
            sparse_prop /= denom;
            row_prop /= denom;
        }

        for (int d = 0; d < 2; ++d) {
            stored_uses_oracle[d] = false;
            for (const auto& x : mats) {
                if (x->uses_oracle(d)) {
                    stored_uses_oracle[d] = true;
                    break;
                }
            }
        }
    }

    /**
     * @param ps Pointers to the matrices to be combined.
     */
    DelayedBind(const std::vector<std::shared_ptr<Matrix<Value_, Index_> > >& ps) : DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >(ps.begin(), ps.end())) {}

private:
    std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > mats;
    std::vector<Index_> cumulative;

    double sparse_prop = 0, row_prop = 0;
    std::array<bool, 2> stored_uses_oracle;

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
        return sparse_prop > 0.5;
    }

    double sparse_proportion() const {
        return sparse_prop;
    }

    bool prefer_rows() const {
        return row_prop > 0.5;
    }

    double prefer_rows_proportion() const {
        return row_prop;
    }

    bool uses_oracle(bool row) const {
        return stored_uses_oracle[row];
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;


private:
    template<DimensionSelectionType selection_>
    struct DenseParallelExtractor : public ParallelExtractor<selection_, false> {
        template<typename ... Args_>
        DenseParallelExtractor(const DelayedBind* p, const Options& opt, Args_&& ... args) : ParallelExtractor<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

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
        SparseParallelExtractor(const DelayedBind* p, const Options& opt, Args_&& ... args) : 
            ParallelExtractor<selection_, true>(p, opt, std::forward<Args_>(args)...),
            needs_value(opt.sparse_extract_value), 
            needs_index(opt.sparse_extract_index)
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

                    if (adjustment) {
                        for (Index_ j = 0; j < raw.number; ++j) {
                            icopy[j] += adjustment;
                        }
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
        static constexpr bool accrow_ = margin_ == 0;

        PerpendicularExtractor(const DelayedBind* p, const Options& opt) : parent(p) {
            workspaces.reserve(parent->mats.size());

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? p->internal_ncol() : p->internal_nrow());
                for (const auto& m : parent->mats) {
                    // TODO: manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), opt));
                }
            }
        }

        PerpendicularExtractor(const DelayedBind* p, const Options& opt, Index_ bs, Index_ bl) : parent(p) {
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

        PerpendicularExtractor(const DelayedBind* p, const Options& opt, std::vector<Index_> idx) : parent(p) {
            workspaces.reserve(p->mats.size());

            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                this->index_length = idx.size();
                for (const auto& m : parent->mats) {
                    // TODO: manage opt.access.sequencer!
                    workspaces.push_back(new_extractor<accrow_, sparse_>(m.get(), idx, opt)); // copy of 'idx' is deliberate here.
                }

                if (workspaces.empty()) {
                    indices = std::move(idx);
                }
            }
        }

    protected:
        const DelayedBind* parent;
        std::vector<std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > > workspaces;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
        size_t last_segment = 0;

    private:
        static size_t choose_segment_raw(Index_ i, const std::vector<Index_>& cumulative) {
            return std::upper_bound(cumulative.begin(), cumulative.end(), i) - cumulative.begin() - 1;
        }

        static void choose_segment(Index_ i, size_t& last_segment, const std::vector<Index_>& cumulative) {
            if (cumulative[last_segment] > i) {
                if (last_segment && cumulative[last_segment - 1] <= i) {
                    --last_segment;
                } else {
                    last_segment = choose_segment_raw(i, cumulative);
                }
            } else if (cumulative[last_segment + 1] <= i) {
                if (last_segment + 2 < cumulative.size() && cumulative[last_segment + 2] > i) {
                    ++last_segment;
                } else {
                    last_segment = choose_segment_raw(i, cumulative);
                }
            }
            return;
        }

    protected:
        size_t choose_segment(Index_ i) {
            choose_segment(i, last_segment, parent->cumulative);
            return last_segment;
        }

    public:
        const Index_* index_start() const {
            if (workspaces.empty()) {
                return indices.data();
            } else {
                return workspaces.front()->index_start();
            }
        }

    private:
        struct ParentOracle {
            ParentOracle(std::unique_ptr<Oracle<Index_> > o, std::vector<unsigned char> u, const std::vector<Index_>* c) : 
                source(std::move(o)), streams(u.size()), used(std::move(u)), cumulative(c) {}

            size_t fill(size_t id, Index_* buffer, size_t number) {
                auto& stream = streams[id];

                // Trying to top it up, but we won't try too hard, as we might
                // end up doing lots of predictions to fill up the stream for a
                // bound matrix with very few elements.
                if (stream.size() < number) {
                    const auto& cum = *cumulative;
                    bool unfound = true;
                    do {
                        size_t filled = source->predict(buffer, number);
                        if (!filled) {
                            break;
                        }

                        for (size_t f = 0; f < filled; ++f) {
                            PerpendicularExtractor::choose_segment(buffer[f], chosen_segment, cum);
                            if (used[chosen_segment]) {
                                streams[chosen_segment].push_back(buffer[f] - cum[chosen_segment]);
                                if (chosen_segment == id) {
                                    unfound = true;
                                }
                            }
                        }
                    } while (unfound);

                    if (stream.size() < number) {
                        number = stream.size();
                    }
                }

                if (number) {
                    std::copy(stream.begin(), stream.begin() + number, buffer);
                    stream.erase(stream.begin(), stream.begin() + number);
                }
                return number;
            }
        private:
            std::unique_ptr<Oracle<Index_> > source;
            std::vector<std::deque<Index_> > streams;
            std::vector<unsigned char> used;
            const std::vector<Index_>* cumulative;
            size_t chosen_segment = 0;
        };

        struct ChildOracle : public Oracle<Index_> {
            ChildOracle(ParentOracle* o, size_t i) : parent(o), id(i) {}
            size_t predict(Index_* buffer, size_t number) {
                return parent->fill(id, buffer, number);
            }
        private:
            ParentOracle* parent;
            size_t id;
        };

        std::unique_ptr<ParentOracle> parent_oracle;

    public:
        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            // Figure out how many of these need an oracle.
            std::vector<size_t> need_oracles;
            size_t nmats = parent->mats.size();
            need_oracles.reserve(nmats);

            for (size_t m = 0; m < nmats; ++m) {
                if (parent->mats[m]->uses_oracle(accrow_)) {
                    need_oracles.push_back(m);
                }
            }

            // Now we create child oracles that reference the parent;
            // or, if only one inner matrix needs an oracle, we pass it directly.
            if (!need_oracles.empty()) {
                std::vector<unsigned char> used(nmats);
                for (auto n : need_oracles) {
                    used[n] = 1;
                }
                parent_oracle.reset(new ParentOracle(std::move(o), std::move(used), &(parent->cumulative))); 
                for (auto n : need_oracles) {
                    workspaces[n]->set_oracle(std::make_unique<ChildOracle>(parent_oracle.get(), n));
                }
            }
        }
    };

private:
    template<DimensionSelectionType selection_>
    struct DensePerpendicularExtractor : public PerpendicularExtractor<selection_, false> {
        template<typename ... Args_>
        DensePerpendicularExtractor(const DelayedBind* p, const Options& opt, Args_&& ... args) : PerpendicularExtractor<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            size_t chosen = this->choose_segment(i);
            return this->workspaces[chosen]->fetch(i - this->parent->cumulative[chosen], buffer);
        }
    };

    template<DimensionSelectionType selection_>
    struct SparsePerpendicularExtractor : public PerpendicularExtractor<selection_, true> {
        template<typename ... Args_>
        SparsePerpendicularExtractor(const DelayedBind* p, const Options& opt, Args_&& ... args) : PerpendicularExtractor<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            size_t chosen = this->choose_segment(i);
            return this->workspaces[chosen]->fetch(i - this->parent->cumulative[chosen], vbuffer, ibuffer);
        }
    };

    /************************************
     ********** Public methods **********
     ************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options& opt, Args_&&... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(sparse_) {
            if constexpr(accrow_ == (margin_ == 0)) {
                output.reset(new SparsePerpendicularExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            } else {
                output.reset(new SparseParallelExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            }
        } else {
            if constexpr(accrow_ == (margin_ == 0)) {
                output.reset(new DensePerpendicularExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            } else {
                output.reset(new DenseParallelExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            }
        }

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam margin_ Dimension along which the combining is to occur.
 * If 0, matrices are combined along the rows; if 1, matrices are combined to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 *
 * @param ps Pointers to `Matrix` objects.
 *
 * @return A pointer to a `DelayedBind` instance.
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<Matrix<Value_, Index_> > > ps) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBind<margin_, Value_, Index_>(std::move(ps)));
}

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam margin_ Dimension along which the combining is to occur.
 * If 0, matrices are combined along the rows; if 1, matrices are combined to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 *
 * @param ps Pointers to `const` `Matrix` objects.
 *
 * @return A pointer to a `DelayedBind` instance.
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > ps) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBind<margin_, Value_, Index_>(std::move(ps)));
}

}

#endif
