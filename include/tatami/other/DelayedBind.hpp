#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "../base/Matrix.hpp"
#include "../utils/new_extractor.hpp"
#include "../utils/ConsecutiveOracle.hpp"
#include "../utils/FixedOracle.hpp"
#include "../utils/PseudoOracularExtractor.hpp"
#include "../utils/copy.hpp"

#include <numeric>
#include <algorithm>
#include <memory>
#include <array>
#include <type_traits>

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

/**********************
 *** Dense parallel ***
 **********************/

template<typename Index_>
size_t find_segment(const std::vector<Index_>& cumulative, Index_ target) {
    return std::upper_bound(cumulative.begin(), cumulative.end(), target) - cumulative.begin() - 1;
}

template<typename Index_, class Initialize_>
size_t initialize_parallel_block(
    const std::vector<Index_>& cumulative, 
    Index_ block_start, 
    Index_ block_length, 
    Initialize_ init) 
{
    size_t start_index = find_segment(cumulative, block_start); // finding the first matrix.
    Index_ actual_start = block_start - cumulative[start_index];
    Index_ block_end = block_start + block_length;

    for (auto index = start_index, nmats = cumulative.size() - 1; index < nmats; ++index) {
        bool not_final = (block_end > cumulative[index + 1]);
        Index_ actual_end = (not_final ? cumulative[index + 1] : block_end) - cumulative[index];
        init(index, actual_start, actual_end - actual_start);
        if (!not_final) {
            break;
        }
        actual_start = 0;
    }

    return start_index;
}

template<typename Index_, class Initialize_>
void initialize_parallel_index(
    const std::vector<Index_>& cumulative, 
    const std::vector<Index_>& indices, 
    Initialize_ init) 
{
    if (indices.empty()) {
        return;
    } 
    size_t start_index = find_segment(cumulative, indices.front());

    Index_ counter = 0, il = indices.size();
    for (size_t index = start_index, nmats = cumulative.size() -1; index < nmats; ++index) {
        Index_ lower = cumulative[index];
        Index_ upper = cumulative[index + 1];

        auto slice_ptr = std::make_shared<std::vector<Index_> >();
        auto& curslice = *slice_ptr;
        while (counter < il && indices[counter] < upper) {
            curslice.push_back(indices[counter] - lower);
            ++counter;
        }

        if (!curslice.empty()) {
            init(index, std::move(slice_ptr));
        }

        if (counter == il) {
            break;
        }
    }
}

template<bool oracle_, typename Value_, typename Index_>
struct ParallelDense : public DenseExtractor<oracle_, Value_, Index_> {
    ParallelDense(
        const std::vector<Index_>&, // Not used, just provided for consistency with other constructors.
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt)
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        for (const auto& m : mats) {
            count.emplace_back(row ? m->ncol() : m->nrow());
            internal.emplace_back(new_extractor<false, oracle_>(m.get(), row, oracle, opt));
        }
    }

    ParallelDense(
        const std::vector<Index_>& cumulative, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt)
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        initialize_parallel_block(
            cumulative, 
            block_start, 
            block_length,
            [&](size_t i, Index_ s, Index_ l) {
                count.emplace_back(l);
                internal.emplace_back(new_extractor<false, oracle_>(mats[i].get(), row, oracle, s, l, opt));
            }
        );
    }

    ParallelDense(
        const std::vector<Index_>& cumulative, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt)
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        initialize_parallel_index(
            cumulative, 
            *indices_ptr,
            [&](size_t i, VectorPtr<Index_> idx) {
                count.emplace_back(idx->size());
                internal.emplace_back(new_extractor<false, oracle_>(mats[i].get(), row, oracle, std::move(idx), opt));
            }
        );
    }

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto copy = buffer;
        for (size_t x = 0, end = count.size(); x < end; ++x) {
            auto ptr = internal[x]->fetch(i, copy); 
            auto num = count[x];
            copy_n(ptr, num, copy);
            copy += num;
        }
        return buffer;
    }

private:
    std::vector<std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > > internal;
    std::vector<Index_> count;
};

/***********************
 *** Sparse parallel ***
 ***********************/

template<bool oracle_, typename Value_, typename Index_>
struct ParallelFullSparse : public SparseExtractor<oracle_, Value_, Index_> {
    ParallelFullSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt) : 
        cumulative(cum), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index)
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(new_extractor<true, oracle_>(m.get(), row, oracle, opt));
        }
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        Index_ accumulated = 0;

        for (size_t x = 0, end = cumulative.size() - 1; x < end; ++x) {
            auto range = internal[x]->fetch(i, vcopy, icopy); 
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

private:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > > internal;
};

template<bool oracle_, typename Value_, typename Index_>
struct ParallelBlockSparse : public SparseExtractor<oracle_, Value_, Index_> {
    ParallelBlockSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        cumulative(cum), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {
        internal.reserve(mats.size());
        which_start = initialize_parallel_block(
            cumulative, 
            block_start, 
            block_length,
            [&](size_t i, Index_ s, Index_ l) {
                internal.emplace_back(new_extractor<true, oracle_>(mats[i].get(), row, oracle, s, l, opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        Index_ count = 0;

        for (size_t x = 0, end = internal.size(); x < end; ++x) {
            auto range = internal[x]->fetch(i, vcopy, icopy);
            count += range.number;
            if (needs_value) {
                copy_n(range.value, range.number, vcopy);
                vcopy += range.number;
            }
            if (needs_index) {
                Index_ offset = cumulative[x + which_start];
                for (Index_ y = 0; y < range.number; ++y) {
                    icopy[y] = range.index[y] + offset;
                }
                icopy += range.number;
            }
        }

        return SparseRange<Value_, Index_>(count, (needs_value ? vbuffer : NULL), (needs_index ? ibuffer : NULL));
    }

private:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > > internal;
    size_t which_start;
};

template<bool oracle_, typename Value_, typename Index_>
struct ParallelIndexSparse : public SparseExtractor<oracle_, Value_, Index_> {
    ParallelIndexSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt) : 
        cumulative(cum),
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {
        internal.reserve(mats.size());
        which.reserve(mats.size());
        initialize_parallel_index(
            cumulative, 
            *indices_ptr,
            [&](size_t i, VectorPtr<Index_> idx) {
                which.emplace_back(i);
                internal.emplace_back(new_extractor<true, oracle_>(mats[i].get(), row, oracle, std::move(idx), opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        Index_ count = 0;

        for (size_t x = 0, end = which.size(); x < end; ++x) {
            auto range = internal[x]->fetch(i, vcopy, icopy);
            count += range.number;
            if (needs_value) {
                copy_n(range.value, range.number, vcopy);
                vcopy += range.number;
            }

            if (needs_index) {
                Index_ offset = cumulative[which[x]];
                for (Index_ y = 0; y < range.number; ++y) {
                    icopy[y] = range.index[y] + offset;
                }
                icopy += range.number;
            }
        }

        return SparseRange<Value_, Index_>(count, (needs_value ? vbuffer : NULL), (needs_index ? ibuffer : NULL));
    }

private:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > > internal;
    std::vector<size_t> which;
};

/*********************
 *** Perpendicular ***
 *********************/

template<typename Index_>
struct ChooseSegment {
    ChooseSegment(const std::vector<Index_>& cum) : cumulative(cum) {}

    const std::vector<Index_>& cumulative;
private:
    size_t last_segment = 0;

public:
    size_t choose_segment(Index_ i) {
        if (cumulative[last_segment] > i) {
            if (last_segment && cumulative[last_segment - 1] <= i) {
                --last_segment;
            } else {
                last_segment = find_segment(cumulative, i);
            }
        } else if (cumulative[last_segment + 1] <= i) {
            if (last_segment + 2 < cumulative.size() && cumulative[last_segment + 2] > i) {
                ++last_segment;
            } else {
                last_segment = find_segment(cumulative, i);
            }
        }
        return last_segment;
    }
};

template<typename Value_, typename Index_>
struct MyopicPerpendicularDense : public MyopicDenseExtractor<Value_, Index_> {
    template<typename ... Args_>
    MyopicPerpendicularDense(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row, 
        const Args_& ... args) : 
        chooser(cum)
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(m->dense(row, args...));
        }
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        size_t chosen = chooser.choose_segment(i);
        return internal[chosen]->fetch(i - chooser.cumulative[chosen], buffer);
    }

private:
    ChooseSegment<Index_> chooser;
    std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_>
struct MyopicPerpendicularSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<typename ... Args_>
    MyopicPerpendicularSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        const Args_& ... args) : 
        chooser(cum)
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(m->sparse(row, args...));
        }
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        size_t chosen = chooser.choose_segment(i);
        return internal[chosen]->fetch(i - chooser.cumulative[chosen], vbuffer, ibuffer);
    }

private:
    ChooseSegment<Index_> chooser;
    std::vector<std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > > internal;
};

template<typename Index_, class Initialize_>
void initialize_perp_oracular(
    const std::vector<Index_>& cumulative, 
    const Oracle<Index_>* oracle, 
    std::vector<size_t>& chosen, 
    Initialize_ init) 
{
    size_t ntotal = oracle->total();
    chosen.reserve(ntotal);
    ChooseSegment<Index_> chooser(cumulative);

    struct Predictions {
        bool consecutive = true;
        Index_ start = 0;
        Index_ number = 0;
        std::vector<Index_> predictions;

        void add(Index_ p) {
            if (consecutive) {
                if (number == 0) {
                    start = p;
                    number = 1;
                    return;
                }
                if (number + start == p) {
                    ++number;
                    return;
                }
                consecutive = false;
                predictions.resize(number);
                std::iota(predictions.begin(), predictions.end(), start);
            }

            predictions.push_back(p);
        }
    };

    size_t nmats = cumulative.size() - 1;
    std::vector<Predictions> predictions(nmats);
    for (size_t i = 0; i < ntotal; ++i) {
        auto prediction = oracle->get(i);
        size_t choice = chooser.choose_segment(prediction);
        chosen.push_back(choice);
        predictions[choice].add(prediction - cumulative[choice]);
    }

    for (size_t x = 0; x < nmats; ++x) {
        auto& current = predictions[x];
        if (current.consecutive) {
            if (current.number) {
                init(x, std::make_shared<ConsecutiveOracle<Index_> >(current.start, current.number));
            }
        } else {
            if (!current.predictions.empty()) {
                init(x, std::make_shared<FixedVectorOracle<Index_> >(std::move(current.predictions)));
            }
        }
    }
}

template<typename Value_, typename Index_>
struct OracularPerpendicularDense : public OracularDenseExtractor<Value_, Index_> {
    template<typename ... Args_>
    OracularPerpendicularDense(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        std::shared_ptr<const Oracle<Index_> > ora, 
        const Args_& ... args) : 
        cumulative(cum)
    {
        internal.resize(mats.size());
        initialize_perp_oracular(
            cum,
            ora.get(),
            segments,
            [&](size_t x, std::shared_ptr<const Oracle<Index_> > subora) {
                internal[x] = mats[x]->dense(row, std::move(subora), args...);
            }
        );
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto chosen = segments[used];
        auto output = internal[chosen]->fetch(i, buffer);
        ++used;
        return output;
    }

private:
    const std::vector<Index_>& cumulative;
    std::vector<size_t> segments;
    std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > internal;
    size_t used = 0;
};

template<typename Value_, typename Index_>
struct OracularPerpendicularSparse : public OracularSparseExtractor<Value_, Index_> {
    template<typename ... Args_>
    OracularPerpendicularSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        bool row,
        std::shared_ptr<const Oracle<Index_> > ora, 
        const Args_& ... args) : 
        cumulative(cum)
    {
        internal.resize(mats.size());
        initialize_perp_oracular(
            cum,
            ora.get(),
            segments,
            [&](size_t x, std::shared_ptr<const Oracle<Index_> > subora) {
                internal[x] = mats[x]->sparse(row, std::move(subora), args...);
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto chosen = segments[used];
        auto output = internal[chosen]->fetch(i, vbuffer, ibuffer);
        ++used;
        return output;
    }

private:
    const std::vector<Index_>& cumulative;
    std::vector<size_t> segments;
    std::vector<std::unique_ptr<OracularSparseExtractor<Value_, Index_> > > internal;
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
     * All matrices to be combined should have the same number of columns (if `margin_ == 0`) or rows (otherwise).
     */
    DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > ps) : mats(std::move(ps)), cumulative(mats.size()+1) {
        size_t sofar = 0;
        for (size_t i = 0, nmats = mats.size(); i < nmats; ++i) {
            auto& current = mats[i];
            Index_ primary, secondary;
            if constexpr(margin_ == 0) {
                primary = current->nrow();
                secondary = current->ncol();
            } else {
                primary = current->ncol();
                secondary = current->nrow();
            }

            if (i == 0) {
                otherdim = secondary;
            } else if (otherdim != secondary) {
                throw std::runtime_error("all 'mats' should have the same number of " + (margin_ == 0 ? std::string("columns") : std::string("rows")));
            }

            // Removing the matrices that don't contribute anything,
            // so we don't have to deal with their overhead.
            if (primary > 0) {
                if (sofar != i) {
                    mats[sofar] = std::move(current);
                }
                cumulative[sofar + 1] = cumulative[sofar] + primary;
                ++sofar;
            }
        }

        cumulative.resize(sofar + 1);
        mats.resize(sofar);

        double denom = 0;
        for (const auto& x : mats) {
            double total = static_cast<double>(x->nrow()) * static_cast<double>(x->ncol());
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
     * All matrices to be combined should have the same number of columns (if `margin_ == 0`) or rows (otherwise).
     */
    DelayedBind(const std::vector<std::shared_ptr<Matrix<Value_, Index_> > >& ps) : DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >(ps.begin(), ps.end())) {}

private:
    std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > mats;
    Index_ otherdim = 0;
    std::vector<Index_> cumulative;

    double sparse_prop = 0, row_prop = 0;
    std::array<bool, 2> stored_uses_oracle;

public:
    Index_ nrow() const {
        if constexpr(margin_==0) {
            return cumulative.back();
        } else {
            return otherdim;
        }
    }

    Index_ ncol() const {
        if constexpr(margin_==0) {
            return otherdim;
        } else {
            return cumulative.back();
        }
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

    /**********************************
     ********** Myopic dense **********
     **********************************/
public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->dense(row, opt);
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(cumulative, mats, row, opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<false, Value_, Index_> >(cumulative, mats, row, false, opt);
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->dense(row, block_start, block_length, opt);
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(cumulative, mats, row, block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<false, Value_, Index_> >(cumulative, mats, row, false, block_start, block_length, opt);
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->dense(row, std::move(indices_ptr), opt);
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(cumulative, mats, row, std::move(indices_ptr), opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<false, Value_, Index_> >(cumulative, mats, row, false, std::move(indices_ptr), opt);
        }
    }

    /***********************************
     ********** Myopic sparse **********
     ***********************************/
private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->sparse(row, opt);
        } else  if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(cumulative, mats, row, opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelFullSparse<false, Value_, Index_> >(cumulative, mats, row, false, opt);
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->sparse(row, block_start, block_length, opt);
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(cumulative, mats, row, block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelBlockSparse<false, Value_, Index_> >(cumulative, mats, row, false, block_start, block_length, opt);
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->sparse(row, std::move(indices_ptr), opt);
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(cumulative, mats, row, std::move(indices_ptr), opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelIndexSparse<false, Value_, Index_> >(cumulative, mats, row, false, std::move(indices_ptr), opt);
        }
    }

    /************************************
     ********** Oracular dense **********
     ************************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->dense(row, std::move(oracle), opt);
        } else if (!stored_uses_oracle[row]) {
            return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, opt));
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(cumulative, mats, row, std::move(oracle), opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<true, Value_, Index_> >(cumulative, mats, row, std::move(oracle), opt);
        }
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->dense(row, std::move(oracle), block_start, block_length, opt);
        } else if (!stored_uses_oracle[row]) {
            return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, block_start, block_length, opt));
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(cumulative, mats, row, std::move(oracle), block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<true, Value_, Index_> >(cumulative, mats, row, std::move(oracle), block_start, block_length, opt);
        }
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->dense(row, std::move(oracle), std::move(indices_ptr), opt);
        } else if (!stored_uses_oracle[row]) {
            return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, std::move(indices_ptr), opt));
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(cumulative, mats, row, std::move(oracle), std::move(indices_ptr), opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<true, Value_, Index_> >(cumulative, mats, row, std::move(oracle), std::move(indices_ptr), opt);
        }
    }

    /*************************************
     ********** Oracular sparse **********
     *************************************/
private:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->sparse(row, std::move(oracle), opt);
        } else if (!stored_uses_oracle[row]) {
            return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, opt));
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(cumulative, mats, row, std::move(oracle), opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelFullSparse<true, Value_, Index_> >(cumulative, mats, row, std::move(oracle), opt);
        }
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->sparse(row, std::move(oracle), block_start, block_length, opt);
        } else if (!stored_uses_oracle[row]) {
            return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, block_start, block_length, opt));
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(cumulative, mats, row, std::move(oracle), block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelBlockSparse<true, Value_, Index_> >(cumulative, mats, row, std::move(oracle), block_start, block_length, opt);
        }
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        if (mats.size() == 1) {
            return mats[0]->sparse(row, std::move(oracle), std::move(indices_ptr), opt);
        } else if (!stored_uses_oracle[row]) {
            return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, std::move(indices_ptr), opt));
        } else if (row == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(cumulative, mats, row, std::move(oracle), std::move(indices_ptr), opt);
        } else {
            return std::make_unique<DelayedBind_internal::ParallelIndexSparse<true, Value_, Index_> >(cumulative, mats, row, std::move(oracle), std::move(indices_ptr), opt);
        }
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
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > ps) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBind<margin_, Value_, Index_>(std::move(ps)));
}

/**
 * @cond
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<Matrix<Value_, Index_> > > ps) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBind<margin_, Value_, Index_>(std::move(ps)));
}
/**
 * @endcond
 */

}

#endif
