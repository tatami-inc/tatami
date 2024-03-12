#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "../base/Matrix.hpp"
#include "../base/new_extractor.hpp"
#include "../utils/ConsecutiveOracle.hpp"
#include "../utils/FixedOracle.hpp"
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

/*********************
 *** Full parallel ***
 *********************/

template<typename Value_, typename Index_, class Store_>
void fetch_dense_parallel_full(const std::vector<Index_>& cumulative, Value_* buffer, Store_ store) {
    for (size_t x = 0, end = cumulative.size() - 1; x < end; ++x) {
        Index_ last = cumulative[x];
        auto copy = buffer + last;
        auto ptr = store(x, copy);
        copy_n(ptr, cumulative[x + 1] - last, copy);
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelFullDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicParallelFullDense(
        const std::vector<Index_>& cum,
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>, 
        const Options& opt) : 
        cumulative(cum) 
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(new_extractor<row_, false>(m.get(), opt));
        }
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        fetch_dense_parallel_full(
            cumulative, 
            buffer, 
            [&](size_t x, Value_* copy) -> const Value_* { 
                return internal[x]->fetch(i, copy); 
            }
        );
        return buffer;
    }

    Index_ number() const {
        return cumulative.back();
    }

private:
    const std::vector<Index_>& cumulative;
    std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_>
struct OracularParallelFullDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_>
    OracularParallelFullDense(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        std::shared_ptr<Oracle<Index_> > ora, 
        const Options& opt) : 
        cumulative(cum) 
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(new_extractor<row_, false>(m.get(), ora, opt));
        }
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        // We can assume that 'internal[x]->fetch()' is called at least once,
        // otherwise we should have created an OracularParallelDenseEmpty.
        fetch_dense_parallel_full(
            cumulative, 
            buffer, 
            [&](size_t x, Value_* copy) -> const Value_* {
                return internal[x]->fetch(i, copy);
            }
        );
        return buffer;
    }

    Index_ number() const {
        return cumulative.back();
    }

private:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_, class Store_>
SparseRange<Value_, Index_> fetch_sparse_full(
    const std::vector<Index_>& cumulative, 
    bool needs_value, 
    bool needs_index, 
    Value_* vbuffer, 
    Index_* ibuffer, 
    Store_ store) 
{
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

    return SparseRange<Value_, Index_>(
        accumulated, 
        (needs_value ? vbuffer : NULL), 
        (needs_index ? ibuffer : NULL)
    );
}

template<typename Value_, typename Index_>
struct MyopicParallelFullSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicParallelFullSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        const Options& opt) : 
        cumulative(cum), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index)
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(new_extractor<row_, true>(m.get(), opt));
        }
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return fetch_sparse_full(
            cumulative, 
            needs_value, 
            needs_index, 
            vbuffer, 
            ibuffer, 
            [&](size_t x, Value_* vcopy, Index_* icopy) -> SparseRange<Value_, Index_> {
                return internal[x]->fetch(i, vcopy, icopy); 
            }
        );
    }

    Index_ number() const {
        return cumulative.back();
    }

private:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_>
struct OracularParallelFullSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_>
    OracularParallelFullSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        std::shared_ptr<Oracle<Index_> > ora,
        const Options& opt) :
        cumulative(cum), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index)
    {
        internal.reserve(mats.size());
        for (const auto& m : mats) {
            internal.emplace_back(new_extractor<row_, true>(m.get(), ora, opt));
        }
    }

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        // We can assume that 'internal[x]->fetch()' is called at least once,
        // otherwise we should have created an OracularParallelSparseEmpty.
        return fetch_sparse_full(
            cumulative,
            needs_value,
            needs_index,
            vbuffer, 
            ibuffer, 
            [&](size_t x, Value_* vcopy, Index_* icopy) -> SparseRange<Value_, Index_> { 
                return internal[x]->fetch(i, vcopy, icopy); 
            }
        );
    }

    Index_ number() const {
        return cumulative.back();
    }

private:
    const std::vector<Index_>& cumulative;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<OracularSparseExtractor<Value_, Index_> > > internal;
};

/**********************
 *** Block parallel ***
 **********************/

template<typename Index_>
size_t find_segment(const std::vector<Index_>& cumulative, Index_ target) {
    return std::upper_bound(cumulative.begin(), cumulative.end(), target) - cumulative.begin() - 1;
}

template<typename Value_, typename Index_, class Initialize_>
size_t initialize_parallel_block(
    const std::vector<Index_>& cumulative, 
    const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
    Index_ block_start, 
    Index_ block_length, 
    Initialize_ init) 
{
    size_t start_index = find_segment(cumulative, block_start); // finding the first matrix.
    Index_ actual_start = block_start - cumulative[start_index];
    Index_ block_end = block_start + block_length;

    size_t nmats = mats.size();
    for (auto index = start_index; index < nmats; ++index) {
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

template<typename Value_, typename Index_, class Store_>
void fetch_dense_parallel_block(const std::vector<Index_>& count, Value_* buffer, Store_ store) {
    auto copy = buffer;
    for (size_t x = 0, end = count.size(); x < end; ++x) {
        auto ptr = store(x, copy);
        auto j = count[x];
        copy_n(ptr, j, copy);
        copy += j;
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicParallelBlockDense(
        const std::vector<Index_>& cumulative, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        block_length(block_length) 
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        initialize_parallel_block(
            cumulative, 
            mats, 
            block_start, 
            block_length,
            [&](size_t i, Index_ s, Index_ l) {
                count.emplace_back(l);
                internal.emplace_back(new_extractor<row_, false>(mats[i].get(), s, l, opt));
            }
        );
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        fetch_dense_parallel_block(
            count, 
            buffer, 
            [&](size_t x, Value_* copy) -> const Value_* { 
                return internal[x]->fetch(i, copy); 
            }
        );
        return buffer;
    }

    Index_ number() const {
        return block_length;
    }

private:
    Index_ block_length;
    std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > internal;
    std::vector<Index_> count;
};

template<typename Value_, typename Index_>
struct OracularParallelBlockDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_>
    OracularParallelBlockDense(
        const std::vector<Index_>& cumulative, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        std::shared_ptr<Oracle<Index_> > ora, 
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        block_length(block_length) 
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        initialize_parallel_block(
            cumulative, 
            mats, 
            block_start, 
            block_length,
            [&](size_t i, Index_ s, Index_ l) {
                count.emplace_back(l);
                internal.emplace_back(new_extractor<row_, false>(mats[i].get(), ora, s, l, opt));
            }
        );
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        // We can assume that 'internal[x]->fetch()' is called at least once,
        // otherwise we should have created an OracularParallelDenseEmpty.
        fetch_dense_parallel_block(
            count, 
            buffer, 
            [&](size_t x, Value_* copy) -> const Value_* { 
                return internal[x]->fetch(i, copy); 
            }
        );
        return buffer;
    }

    Index_ number() const {
        return block_length;
    }

private:
    Index_ block_length;
    std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > internal;
    std::vector<Index_> count;
    std::shared_ptr<Oracle<Index_> > oracle;
    size_t used = 0;
};

template<typename Value_, typename Index_, class Store_>
SparseRange<Value_, Index_> fetch_sparse_parallel_block(
    const std::vector<Index_>& cumulative, 
    size_t which_start, 
    size_t num, 
    bool needs_value, 
    bool needs_index, 
    Value_* vbuffer, 
    Index_* ibuffer, 
    Store_ store) 
{
    auto vcopy = vbuffer;
    auto icopy = ibuffer;
    Index_ count = 0;

    for (size_t x = 0; x < num; ++x) {
        auto range = store(x, vcopy, icopy);
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

    return SparseRange<Value_, Index_>(
        count,
        (needs_value ? vbuffer : NULL),
        (needs_index ? ibuffer : NULL)
    );
}

template<typename Value_, typename Index_>
struct MyopicParallelBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicParallelBlockSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        cumulative(cum), 
        block_length(block_length), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {
        internal.reserve(mats.size());
        which_start = initialize_parallel_block(
            cumulative, 
            mats, 
            block_start, 
            block_length,
            [&](size_t i, Index_ s, Index_ l) {
                internal.emplace_back(new_extractor<row_, true>(mats[i].get(), s, l, opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return fetch_sparse_parallel_block(
            cumulative, 
            which_start,
            internal.size(),
            needs_value, 
            needs_index, 
            vbuffer, 
            ibuffer, 
            [&](size_t x, Value_* vcopy, Index_* icopy) -> SparseRange<Value_, Index_> {
                return internal[x]->fetch(i, vcopy, icopy);
            }
        );
    }

    Index_ number() const {
        return block_length;
    }

private:
    const std::vector<Index_>& cumulative;
    Index_ block_length;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > > internal;
    size_t which_start;
};

template<typename Value_, typename Index_>
struct OracularParallelBlockSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_>
    OracularParallelBlockSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        std::shared_ptr<Oracle<Index_> > ora, 
        Index_ block_start, 
        Index_ block_length, 
        const Options& opt) : 
        cumulative(cum), 
        block_length(block_length), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index)
    {
        internal.reserve(mats.size());
        which_start = initialize_parallel_block(
            cumulative, 
            mats, 
            block_start, 
            block_length,
            [&](size_t i, Index_ s, Index_ l) {
                internal.emplace_back(new_extractor<row_, true>(mats[i].get(), ora, s, l, opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        // We can assume that 'internal[x]->fetch()' is called at least once,
        // otherwise we should have created an OracularParallelSparseEmpty.
        return fetch_sparse_parallel_block(
            cumulative, 
            which_start,
            internal.size(),
            needs_value, 
            needs_index, 
            vbuffer, 
            ibuffer, 
            [&](size_t x, Value_* vcopy, Index_* icopy) -> SparseRange<Value_, Index_> {
                return internal[x]->fetch(i, vcopy, icopy);
            }
        );
    }

    Index_ number() const {
        return block_length;
    }

private:
    const std::vector<Index_>& cumulative;
    Index_ block_length;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<OracularSparseExtractor<Value_, Index_> > > internal;
    size_t which_start;
    std::shared_ptr<Oracle<Index_> > oracle;
    size_t used = 0;
};

/**********************
 *** Index parallel ***
 **********************/

template<typename Value_, typename Index_, class Initialize_>
void initialize_parallel_index(
    const std::vector<Index_>& cumulative, 
    const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
    const std::vector<Index_>& indices, 
    Initialize_ init) 
{
    if (indices.empty()) {
        return;
    } 
    size_t start_index = find_segment(cumulative, indices.front());
    size_t nmats = cumulative.size();

    Index_ counter = 0, il = indices.size();
    for (size_t index = 0; index < nmats; ++index) {
        Index_ lower = cumulative[index];
        Index_ upper = cumulative[index + 1];

        std::vector<Index_> curslice;
        while (counter < il && indices[counter] < upper) {
            curslice.push_back(indices[counter] - lower);
            ++counter;
        }

        if (!curslice.empty()) {
            init(index, std::move(curslice));
        }

        if (counter == il) {
            break;
        }
    }
}

template<typename Value_, typename Index_, class Store_>
void fetch_dense_parallel_index(const std::vector<Index_>& count, Value_* buffer, Store_ store) {
    auto copy = buffer;
    for (size_t x = 0, end = count.size(); x < end; ++x) {
        auto ptr = store(x, copy);
        auto j = count[x];
        copy_n(ptr, j, copy);
        copy += j;
    }
}

template<typename Value_, typename Index_>
struct MyopicParallelIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicParallelIndexDense(
        const std::vector<Index_>& cumulative, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        const std::vector<Index_>& indices,
        const Options& opt) :
        index_length(indices.size()) 
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        initialize_parallel_index(
            cumulative, 
            mats, 
            indices,
            [&](size_t i, std::vector<Index_> idx) {
                count.emplace_back(idx.size());
                internal.emplace_back(new_extractor<row_, false>(mats[i].get(), std::move(idx), opt));
            }
        );
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        fetch_dense_parallel_index(
            count, 
            buffer, 
            [&](size_t x, Value_* copy) -> const Value_* { 
                return internal[x]->fetch(i, copy); 
            }
        );
        return buffer;
    }

    Index_ number() const {
        return index_length;
    }

private:
    Index_ index_length;
    std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > internal;
    std::vector<Index_> count;
};

template<typename Value_, typename Index_>
struct OracularParallelIndexDense : public OracularDenseExtractor<Value_, Index_> {
    template<bool row_>
    OracularParallelIndexDense(
        const std::vector<Index_>& cumulative, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        std::shared_ptr<Oracle<Index_> > ora, 
        const std::vector<Index_>& indices,
        const Options& opt) : 
        index_length(indices.size()) 
    {
        internal.reserve(mats.size());
        count.reserve(mats.size());
        initialize_parallel_index(
            cumulative, 
            mats, 
            indices,
            [&](size_t i, std::vector<Index_> idx) {
                count.emplace_back(idx.size());
                internal.emplace_back(new_extractor<row_, false>(mats[i].get(), ora, std::move(idx), opt));
            }
        );
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        // We can assume that 'internal[x]->fetch()' is called at least once,
        // otherwise we should have created an OracularParallelDenseEmpty.
        fetch_dense_parallel_index(
            count, 
            buffer, 
            [&](size_t x, Value_* copy) -> const Value_* { 
                return internal[x]->fetch(i, copy);
            }
        );
        return buffer;
    }

    Index_ number() const {
        return index_length;
    }

private:
    Index_ index_length;
    std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > internal;
    std::vector<Index_> count;
};

template<typename Value_, typename Index_, class Store_>
SparseRange<Value_, Index_> fetch_sparse_parallel_index(
    const std::vector<Index_>& cumulative, 
    const std::vector<size_t>& which, 
    bool needs_value, 
    bool needs_index, 
    Value_* vbuffer, 
    Index_* ibuffer, 
    Store_ store) 
{
    auto vcopy = vbuffer;
    auto icopy = ibuffer;
    Index_ count = 0;

    for (size_t x = 0, end = which.size(); x < end; ++x) {
        auto range = store(x, vcopy, icopy);
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

    return SparseRange<Value_, Index_>(
        count,
        (needs_value ? vbuffer : NULL),
        (needs_index ? ibuffer : NULL)
    );
}

template<typename Value_, typename Index_>
struct MyopicParallelIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_>
    MyopicParallelIndexSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        const std::vector<Index_>& indices,
        const Options& opt) : 
        cumulative(cum),
        index_length(indices.size()), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {
        internal.reserve(mats.size());
        which.reserve(mats.size());
        initialize_parallel_index(
            cumulative, 
            mats, 
            indices,
            [&](size_t i, std::vector<Index_> idx) {
                which.emplace_back(i);
                internal.emplace_back(new_extractor<row_, true>(mats[i].get(), std::move(idx), opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return fetch_sparse_parallel_index(
            cumulative,
            which, 
            needs_value, 
            needs_index, 
            vbuffer, 
            ibuffer, 
            [&](size_t x, Value_* vcopy, Index_* icopy) -> SparseRange<Value_, Index_> {
                return internal[x]->fetch(i, vcopy, icopy);
            }
        );
    }

    Index_ number() const {
        return index_length;
    }

private:
    const std::vector<Index_>& cumulative;
    Index_ index_length;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > > internal;
    std::vector<size_t> which;
};

template<typename Value_, typename Index_>
struct OracularParallelIndexSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_>
    OracularParallelIndexSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>,
        std::shared_ptr<Oracle<Index_> > ora, 
        const std::vector<Index_>& indices,
        const Options& opt) : 
        cumulative(cum),
        index_length(indices.size()), 
        needs_value(opt.sparse_extract_value), 
        needs_index(opt.sparse_extract_index) 
    {
        internal.reserve(mats.size());
        which.reserve(mats.size());
        initialize_parallel_index(
            cumulative, 
            mats, 
            indices,
            [&](size_t i, std::vector<Index_> idx) {
                which.emplace_back(i);
                internal.emplace_back(new_extractor<row_, true>(mats[i].get(), ora, std::move(idx), opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        // We can assume that 'internal[x]->fetch()' is called at least once,
        // otherwise we should have created an OracularParallelSparseEmpty.
        return fetch_sparse_parallel_index(
            cumulative,
            which, 
            needs_value, 
            needs_index, 
            vbuffer, 
            ibuffer, 
            [&](size_t x, Value_* vcopy, Index_* icopy) -> SparseRange<Value_, Index_> {
                return internal[x]->fetch(i, vcopy, icopy);
            }
        );
    }

    Index_ number() const {
        return index_length;
    }

private:
    const std::vector<Index_>& cumulative;
    Index_ index_length;
    bool needs_value, needs_index;
    std::vector<std::unique_ptr<OracularSparseExtractor<Value_, Index_> > > internal;
    std::vector<size_t> which;
};

/**********************
 *** Parallel empty ***
 **********************/

template<typename Value_, typename Index_>
struct OracularParallelDenseEmpty : public OracularDenseExtractor<Value_, Index_> {
    OracularParallelDenseEmpty(std::shared_ptr<Oracle<Index_> > ora, const Options&) : oracle(std::move(ora)) {}

    const Value_* fetch(Index_& i, Value_* buffer) {
        i = oracle->get(used);
        ++used;
        return buffer;
    }

    Index_ number() const {
        return 0;
    }
public:
    std::shared_ptr<Oracle<Index_> > oracle;
    size_t used = 0;
};

template<typename Value_, typename Index_>
struct OracularParallelSparseEmpty : public OracularSparseExtractor<Value_, Index_> {
    OracularParallelSparseEmpty(std::shared_ptr<Oracle<Index_> > ora, const Options& opt) : 
        oracle(std::move(ora)), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {}

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        i = oracle->get(used);
        ++used;
        return SparseRange<Value_, Index_>(0, (needs_value ? vbuffer : NULL), (needs_index ? ibuffer : NULL)); 
    }

    Index_ number() const {
        return 0;
    }
public:
    std::shared_ptr<Oracle<Index_> > oracle;
    bool needs_value, needs_index;
    size_t used = 0;
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
    template<bool row_, typename ... Args_>
    MyopicPerpendicularDense(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>, 
        Index_ extent,
        const Args_& ... args) : 
        chooser(cum), 
        extracted_number(extent)
    {
        internal.reserve(mats.size());
        for (auto& m : mats) {
            internal.emplace_back(new_extractor<row_, false>(m.get(), args...));
        }
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        size_t chosen = chooser.choose_segment(i);
        return internal[chosen]->fetch(i - chooser.cumulative[chosen], buffer);
    }

    Index_ number() const {
        return extracted_number;
    }

private:
    ChooseSegment<Index_> chooser;
    Index_ extracted_number;
    std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > internal;
};

template<typename Value_, typename Index_>
struct MyopicPerpendicularSparse : public MyopicSparseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    MyopicPerpendicularSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>, 
        Index_ extent,
        const Args_& ... args) : 
        chooser(cum), 
        extracted_number(extent)
    {
        internal.reserve(mats.size());
        for (auto& m : mats) {
            internal.emplace_back(new_extractor<row_, true>(m.get(), args...));
        }
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        size_t chosen = chooser.choose_segment(i);
        return internal[chosen]->fetch(i - chooser.cumulative[chosen], vbuffer, ibuffer);
    }

    Index_ number() const {
        return extracted_number;
    }

private:
    ChooseSegment<Index_> chooser;
    Index_ extracted_number;
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
    template<bool row_, typename ... Args_>
    OracularPerpendicularDense(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>, 
        Index_ extent,
        std::shared_ptr<Oracle<Index_> > ora, 
        const Args_& ... args) : 
        cumulative(cum), 
        extracted_number(extent) 
    {
        internal.resize(mats.size());
        initialize_perp_oracular(
            cum,
            ora.get(),
            segments,
            [&](size_t x, std::shared_ptr<Oracle<Index_> > subora) {
                internal[x] = new_extractor<row_, false>(mats[x].get(), std::move(subora), args...);
            }
        );
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        auto chosen = segments[used];
        auto output = internal[chosen]->fetch(i, buffer);
        i += cumulative[chosen];
        ++used;
        return output;
    }

    Index_ number() const {
        return extracted_number;
    }

private:
    const std::vector<Index_>& cumulative;
    Index_ extracted_number;
    std::vector<size_t> segments;
    std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > internal;
    size_t used = 0;
};

template<typename Value_, typename Index_>
struct OracularPerpendicularSparse : public OracularSparseExtractor<Value_, Index_> {
    template<bool row_, typename ... Args_>
    OracularPerpendicularSparse(
        const std::vector<Index_>& cum, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& mats, 
        std::integral_constant<bool, row_>, 
        Index_ extent, 
        std::shared_ptr<Oracle<Index_> > ora, 
        const Args_& ... args) : 
        cumulative(cum), 
        extracted_number(extent) 
    {
        internal.resize(mats.size());
        initialize_perp_oracular(
            cum,
            ora.get(),
            segments,
            [&](size_t x, std::shared_ptr<Oracle<Index_> > subora) {
                internal[x] = new_extractor<row_, true>(mats[x].get(), std::move(subora), args...);
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        auto chosen = segments[used];
        auto output = internal[chosen]->fetch(i, vbuffer, ibuffer);
        i += cumulative[chosen];
        ++used;
        return output;
    }

    Index_ number() const {
        return extracted_number;
    }

private:
    const std::vector<Index_>& cumulative;
    Index_ extracted_number;
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
private:
    template<bool accrow_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_internal(const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(cumulative, mats, flag, otherdim, opt);
        } else {
            return std::make_unique<DelayedBind_internal::MyopicParallelFullDense<Value_, Index_> >(cumulative, mats, flag, opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_internal(Index_ block_start, Index_ block_length, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(cumulative, mats, flag, block_length, block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::MyopicParallelBlockDense<Value_, Index_> >(cumulative, mats, flag, block_start, block_length, opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_internal(std::vector<Index_> indices, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            Index_ extent = indices.size();
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(cumulative, mats, flag, extent, std::move(indices), opt);
        } else {
            return std::make_unique<DelayedBind_internal::MyopicParallelIndexDense<Value_, Index_> >(cumulative, mats, flag, std::move(indices), opt);
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return dense_internal<true>(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<true>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return dense_internal<true>(std::move(indices), opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return dense_internal<false>(opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<false>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return dense_internal<false>(std::move(indices), opt);
    }

    /***********************************
     ********** Myopic sparse **********
     ***********************************/
private:
    template<bool accrow_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_internal(const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(cumulative, mats, flag, otherdim, opt);
        } else {
            return std::make_unique<DelayedBind_internal::MyopicParallelFullSparse<Value_, Index_> >(cumulative, mats, flag, opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_internal(Index_ block_start, Index_ block_length, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(cumulative, mats, flag, block_length, block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::MyopicParallelBlockSparse<Value_, Index_> >(cumulative, mats, flag, block_start, block_length, opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_internal(std::vector<Index_> indices, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            Index_ extent = indices.size();
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(cumulative, mats, flag, extent, std::move(indices), opt);
        } else {
            return std::make_unique<DelayedBind_internal::MyopicParallelIndexSparse<Value_, Index_> >(cumulative, mats, flag, std::move(indices), opt);
        }
    }

public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return sparse_internal<true>(opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<true>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return sparse_internal<true>(std::move(indices), opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return sparse_internal<false>(opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<false>(block_start, block_length, opt);
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return sparse_internal<false>(std::move(indices), opt);
    }

    /************************************
     ********** Oracular dense **********
     ************************************/
private:
    template<bool accrow_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_internal(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(cumulative, mats, flag, otherdim, std::move(oracle), opt);
        } else if (!mats.empty()) {
            return std::make_unique<DelayedBind_internal::OracularParallelFullDense<Value_, Index_> >(cumulative, mats, flag, std::move(oracle), opt);
        } else {
            return std::make_unique<DelayedBind_internal::OracularParallelDenseEmpty<Value_, Index_> >(std::move(oracle), opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_internal(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(cumulative, mats, flag, block_length, std::move(oracle), block_start, block_length, opt);
        } else if (block_length > 0) {
            return std::make_unique<DelayedBind_internal::OracularParallelBlockDense<Value_, Index_> >(cumulative, mats, flag, std::move(oracle), block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::OracularParallelDenseEmpty<Value_, Index_> >(std::move(oracle), opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_internal(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            Index_ extent = indices.size();
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(cumulative, mats, flag, extent, std::move(oracle), std::move(indices), opt);
        } else if (!indices.empty()) {
            return std::make_unique<DelayedBind_internal::OracularParallelIndexDense<Value_, Index_> >(cumulative, mats, flag, std::move(oracle), std::move(indices), opt);
        } else {
            return std::make_unique<DelayedBind_internal::OracularParallelDenseEmpty<Value_, Index_> >(std::move(oracle), opt);
        }
    }

public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return dense_internal<true>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<true>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return dense_internal<true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return dense_internal<false>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return dense_internal<false>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return dense_internal<false>(std::move(oracle), std::move(indices), opt);
    }

    /*************************************
     ********** Oracular sparse **********
     *************************************/
private:
    template<bool accrow_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_internal(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(cumulative, mats, flag, otherdim, std::move(oracle), opt);
        } else if (!mats.empty()) {
            return std::make_unique<DelayedBind_internal::OracularParallelFullSparse<Value_, Index_> >(cumulative, mats, flag, std::move(oracle), opt);
        } else {
            return std::make_unique<DelayedBind_internal::OracularParallelSparseEmpty<Value_, Index_> >(std::move(oracle), opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_internal(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(cumulative, mats, flag, block_length, std::move(oracle), block_start, block_length, opt);
        } else if (block_length > 0) {
            return std::make_unique<DelayedBind_internal::OracularParallelBlockSparse<Value_, Index_> >(cumulative, mats, flag, std::move(oracle), block_start, block_length, opt);
        } else {
            return std::make_unique<DelayedBind_internal::OracularParallelSparseEmpty<Value_, Index_> >(std::move(oracle), opt);
        }
    }

    template<bool accrow_>
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_internal(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        std::integral_constant<bool, accrow_> flag;
        if constexpr(accrow_ == (margin_ == 0)) {
            Index_ extent = indices.size();
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(cumulative, mats, flag, extent, std::move(oracle), std::move(indices), opt);
        } else if (!indices.empty()) {
            return std::make_unique<DelayedBind_internal::OracularParallelIndexSparse<Value_, Index_> >(cumulative, mats, flag, std::move(oracle), std::move(indices), opt);
        } else {
            return std::make_unique<DelayedBind_internal::OracularParallelSparseEmpty<Value_, Index_> >(std::move(oracle), opt);
        }
    }

public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return sparse_internal<true>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<true>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return sparse_internal<true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return sparse_internal<false>(std::move(oracle), opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return sparse_internal<false>(std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return sparse_internal<false>(std::move(oracle), std::move(indices), opt);
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
