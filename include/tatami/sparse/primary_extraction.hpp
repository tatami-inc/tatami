#ifndef TATAMI_SPARSE_PRIMARY_EXTRACTION_HPP
#define TATAMI_SPARSE_PRIMARY_EXTRACTION_HPP

#include "utils.hpp"
#include <utility>
#include <algorithm>

namespace tatami {

namespace sparse_utils {

template<typename Index_, class IndexStorage_, class PointerStorage_>
std::pair<size_t, size_t> extract_primary_dimension(Index_ i, const IndexStorage_& indices, const PointerStorage_& indptrs) {
    auto lower = sparse_utils::get_lower_limit(indptrs, i);
    return std::pair<size_t, size_t>(lower, sparse_utils::get_upper_limit(indices, indptrs, i) - lower);
}

template<typename Index_, class IndexStorage_, class PointerStorage_>
std::pair<size_t, size_t> extract_primary_dimension(
    Index_ i, 
    Index_ start, 
    Index_ length, 
    const IndexStorage_& indices, 
    const PointerStorage_& indptrs, 
    std::vector<std::pair<size_t, size_t> >& cached) 
{
    bool do_cache = !cached.empty();
    if (do_cache) {
        auto val = cached[i];
        if (val.first != -1) {
            return val;
        }
    }

    auto iIt = indices.begin() + sparse_utils::get_lower_limit(indptrs, i);
    auto eIt = indices.begin() + sparse_utils::get_upper_limit(indices, indptrs, i);

    if (iIt != eIt) {
        if (start > *iIt) {
            iIt = std::lower_bound(iIt, eIt, start);
        } 

        auto last = start + length;

        // Comparing the one-past-the-last requested index with the last observed index at 'eIt'.
        // If the former is less than the latter, then we need to do a binary search.
        // If greater than the latter, we restore 'eIt' to its original position.
        // If equal, then the decremented 'eIt' is already the one-past-the-last index, so we keep it as-is.
        --eIt;
        if (last < *eIt) {
            eIt = std::lower_bound(iIt, eIt, last);
        } else if (last > *eIt) {
            ++eIt;
        }
    }

    size_t outstart = iIt - indices.begin();
    size_t outlength = eIt - iIt;
    if (do_cache) {
        cached[i].first = outstart;
        cached[i].second = outlength;
    }

    return std::make_pair(outstart, outlength);
}

template<class ValueStorage_, typename Value_, typename Index_>
void transplant_primary_values(const ValueStorage_& values, const std::pair<size_t, size_t>& range, SparseRange<Value_, Index_>& output, Value_* out_values) {
    if (out_values) {
        if constexpr(has_data<Value_, ValueStorage_>::value) {
            output.value = values.data() + range.first;
        } else {
            auto vIt = values.begin() + range.first;
            std::copy(vIt, vIt + range.second, out_values);
            output.value = out_values;
        }
    } else {
        output.value = NULL;
    }
}

template<class IndexStorage_, typename Value_, typename Index_>
void transplant_primary_indices(const IndexStorage_& indices, const std::pair<size_t, size_t>& range, SparseRange<Value_, Index_>& output, Index_* out_indices) {
    if (out_indices) {
        if constexpr(has_data<Index_, IndexStorage_>::value) {
            output.index = indices.data() + range.first;
        } else {
            auto iIt = indices.begin() + range.first;
            std::copy(iIt, iIt + range.second, out_indices);
            output.index = out_indices;
        }
    } else {
        output.index = NULL;
    }
}

template<class ValueStorage_, class IndexStorage_, typename Value_, typename Index_>
void transplant_primary_expanded(const ValueStorage_& values, const IndexStorage_& indices, const std::pair<size_t, size_t>& range, Value_* out_values, Index_ start, Index_ length) {
    std::fill(out_values, out_values + length, static_cast<Value_>(0));
    auto vIt = values.begin() + range.first;
    auto iIt = indices.begin() + range.first;
    for (size_t x = 0; x < range.second; ++x, ++vIt, ++iIt) {
        out_values[*iIt - start] = *vIt;
    }
    return;
}

template<typename Index_, class IndexStorage_, class PointerStorage_, class Store_>
void primary_dimension(
    Index_ i, 
    const Index_* subset, 
    Index_ length, 
    const IndexStorage_& indices, 
    const PointerStorage_& indptrs, 
    std::vector<size_t>& cached, 
    Store_& store) 
{
    if (!length) {
        return;
    }

    auto iIt = indices.begin() + sparse_utils::get_lower_limit(indptrs, i);
    auto eIt = indices.begin() + sparse_utils::get_upper_limit(indices, indptrs, i);

    if (indices[0]) { // Only jumping ahead if the start is non-zero.
        bool do_cache = !cached.empty();
        if (do_cache) {
            if (cached[i] != -1) { // retrieving the jump from cache, if we came here before.
                iIt += cached[i];
            } else {
                auto iIt2 = std::lower_bound(iIt, eIt, subset[0]);
                cached[i] = iIt2 - iIt;
                iIt = iIt2;
            }
        } else {
            iIt = std::lower_bound(iIt, eIt, subset[0]);
        }
    } 

    if (iIt == eIt) {
        return;
    }

    Index_ counter = 0;
    while (counter < length) {
        auto current = subset[counter];

        while (iIt != eIt && current > *iIt) {
            ++iIt;
        }
        if (iIt == eIt) {
            break;
        }

        if (current == *iIt) {
            store.add(current, iIt - indices.begin());
        } else {
            store.skip(current);
        }
        ++counter;
    }

    return;
}

// Potential classes to use as Store_ in the primary_dimension() indexed extractor.
template<typename Value_, typename Index_, class ValueStorage_>
struct SimpleRawStore {
    SimpleRawStore(const ValueStorage_& iv, Value_* ov, Index_* oi) : in_values(iv), out_values(ov), out_indices(oi) {}

private:
    const ValueStorage_& in_values;
    Value_* out_values;
    Index_* out_indices;

public:
    Index_ n = 0;

    void add(Index_ i, size_t ptr) {
        ++n;
        if (out_indices) {
            *out_indices = i;
            ++out_indices;
        }
        if (out_values) {
            *out_values = in_values[ptr];
            ++out_values;
        }
        return;
    }

    void skip(Index_) {} 
};

template<typename Value_, typename Index_, class ValueStorage_>
struct SimpleExpandedStore {
    SimpleExpandedStore(const ValueStorage_& iv, Value_* ov) : in_values(iv), out_values(ov) {}

private:
    const ValueStorage_& in_values;
    Value_* out_values;

public:
    void add(Index_, size_t ptr) {
        *out_values = in_values[ptr];
        ++out_values;
        return;
    }

    void skip(Index_) {
        ++out_values;
    }
};

}

}

#endif
