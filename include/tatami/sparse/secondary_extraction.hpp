#ifndef TATAMI_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP
#define TATAMI_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP

#include <vector>
#include <type_traits>
#include <algorithm>

#include "utils.hpp"

namespace tatami {

namespace sparse_utils {

template<typename Index_, typename Pointer_ = size_t> 
struct SparseSecondaryExtractionCache {
private:
    // The cached position of the pointer at each primary element.
    // Specifically, 'indices[cached_indptrs[i]]' is the lower bound for 'last_request' in the primary element 'i'.
    std::vector<Pointer_> cached_indptrs; 

    // This vector contains the cached index being pointed to by 'cached_indptrs'.
    // We store this here as it is more cache-friendly than doing a look-up to 'indices' every time.
    //
    // More specifically, if 'last_increasing = true', we define 'cached_indices[i] := indices[cached_indptrs[i]]'.
    // If 'cached_indptrs[i]' is out of range, the cached index is instead set to 'max_index'.
    //
    // On the other hand, if 'last_increasing = false', we define 'cached_indices[i]' to be:
    // - 'last_request + 1', if 'indices[cached_indptrs[i]] == last_request'.
    // - 0, if 'cached_indptrs[i]' already refers to the start of the primary element 'i'.
    // - 'indices[cached_indptrs[i] - 1] + 1' otherwise.
    // This can be considered the "reverse lower bound", i.e., the largest element not greater than 'last_request'.
    std::vector<Index_> cached_indices;

    // Closest value in 'cached_indices' to the 'last_request', to see whether we can short-circuit the iteration.
    // If 'last_increasing = true', this is the minimum of values in 'cached_indices'.
    // If 'last_increasing = false', this is the maximum of values in 'cached_indices'.
    Index_ closest_cached_index = 0;

private:
    Index_ max_index;

    // What was the last requested index on the secondary dimension?
    Index_ last_request = 0;

    bool last_increasing = true;

public:
    SparseSecondaryExtractionCache() = default;

    SparseSecondaryExtractionCache(Index_ mi, Index_ length) : cached_indptrs(length), cached_indices(length), max_index(mi) {}

private:
    template<class ServeIndices_, class Store_>
    void search_above(Index_ secondary, Index_ index_primary, Index_ primary, const ServeIndices_& indices, Store_& store) {
        // Skipping if the curdex (corresponding to curptr) is already higher
        // than secondary. So, we only need to do more work if the request is
        // greater than the stored index. This also catches cases where we're
        // at the end of the dimension, as curdex is set to max_index.
        auto& curdex = cached_indices[index_primary];
        if (curdex > secondary) {
            return;
        }

        auto& curptr = cached_indptrs[index_primary];
        if (curdex == secondary) {
            store(primary, index_primary, cached_indptrs[index_primary]);
            return;
        }

        // Having a peek at the index of the next non-zero element; maybe we're
        // lucky enough that the requested index is below this, as would be the
        // case for consecutive or near-consecutive accesses.
        ++curptr;
        auto istart = indices.start(primary);
        auto iend = indices.end(primary);
        auto inext = istart + curptr;
        if (inext == iend) {
            curdex = max_index;
            return;
        }

        curdex = *inext;
        if (curdex > secondary) {
            return;
        }

        if (curdex == secondary) {
            store(primary, index_primary, curptr);
            return;
        }

        // Otherwise we need to search indices above the existing position. We
        // do a quick increment to cut down the search space a bit more. We
        // also define our own comparator to avoid problems with signed
        // comparisons, depending on the types of IndexIt_ and Index_.
        inext = std::lower_bound(inext + 1, iend, secondary, [](Index_ a, Index_b) -> bool { return a < b; });
        curptr = inext - istart;
        if (curptr == endptr) {
            curdex = max_index;
            return;
        }

        curdex = *inext;
        if (curdex > secondary) {
            return;
        }

        store(primary, index_primary, curptr);
        return;
    }

private:
    template<class ServeIndices_, class Store_>
    void search_below(Index_ secondary, Index_ index_primary, Index_ primary, const ServeIndices_& indices, Store_& store) {
        auto secondaryP1 = secondary + 1;
        auto& curdex = cached_indices[index_primary];
        if (curdex < secondaryP1) {
            return;
        }

        auto& curptr = cached_indptrs[index_primary];
        if (curdex == secondaryP1) {
            store(primary, index_primary, curptr);
            return;
        }

        // Can't decrement any further.
        if (curptr == 0) {
            curdex = 0;
            return;
        }

        --curptr;
        auto istart = indices.start(primary);
        auto inext = istart + curptr;
        curdex = *inext + 1;
        if (curdex < secondaryP1) {
            return;
        }

        if (curdex == secondaryP1) {
            store(primary, index_primary, curptr);
            return;
        }

        // Otherwise we need to search indices below the existing position.
        inext = std::lower_bound(istart, inext, secondary, [](Index_ a, Index_b) -> bool { return a < b; });
        curdex = *inext + 1;
        curptr = inext - istart;

        if (curdex == secondaryP1) {
            store(primary, index_primary, curptr);
            return;
        }

        if (curptr == 0) {
            curdex = 0;
            return;
        }

        --inext;
        curdex = *inext + 1;
        return;
    }

protected:
    template<class Index_, class PrimaryFunction_, class ServeIndices_, class Store_>
    void search_base(Index_ secondary, PrimaryFunction_ to_primary, const ServeIndices_& indices, Store_& store) {
        Index_ primary_length = cached_indices.size(); 
        if (primary_length == 0) {
            return;
        }

        if (secondary > last_request || (last_increasing && secondary == last_request)) {
            if (last_increasing) {
                if (secondary < closest_cached_index) {
                    last_request = secondary;
                    return; 
                }
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_above(secondary, p, to_primary(p), indices, store);
                }

            } else {
                // Need to reset the meaning of 'cached_indices'.
                last_increasing = true;
                for (Index_ p = 0; p < primary_length; ++p) {
                    auto primary = to_primary(p);
                    auto istart = indices.start(primary);
                    auto iend = indices.end(primary);
                    auto icurrent = istart + cached_indptrs[p];
                    cached_indices[p] = (icurrent == iend ? max_index : *icurrent);
                    search_above(secondary, p, primary, indices, store);
                }
            }

            closest_cached_index = *(std::min_element(cached_indices.begin(), cached_indices.end()));

        } else {
            if (!last_increasing) {
                if (secondary + 1 > closest_cached_index) {
                    last_request = secondary;
                    return;
                }
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_below(secondary, p, to_primary(p), indices, store);
                }

            } else {
                // Need to reset the meaning of 'cached_indices'.
                last_increasing = false;
                for (Index_ p = 0; p < primary_length; ++p) {
                    auto primary = to_primary(p);
                    auto istart = indices.start(primary);
                    auto iend = indices.end(primary);
                    auto icurrent = istart + cached_indptrs[p];
                    if (icurrent != iend && *icurrent == last_request) {
                        cached_indices[p] = last_request + 1;
                    } else {
                        cached_indices[p] = (istart == icurrent ? 0 : *(icurrent - 1) + 1);
                    }
                    search_below(secondary, p, primary, indices, store);
                }
            }

            closest_cached_index = *(std::max_element(cached_indices.begin(), cached_indices.end()));
        }

        last_request = secondary; 
        return;
    }
};

}

#endif
