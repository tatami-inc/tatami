#ifndef TATAMI_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP
#define TATAMI_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP

#include <vector>
#include <type_traits>
#include <algorithm>

#include "utils.hpp"

namespace tatami {

namespace sparse_utils {

template<typename Index_, class IndexServer_> 
struct SecondaryExtractionCache {
private:
    Index_ max_index;

    IndexServer_ indices;

    typedef typename IndexServer_::pointer_type Pointer_;

private:
    // The cached position of the pointer at each primary element.
    // Specifically, 'indices[cached_indptrs[i]]' is the lower bound for 'last_request' in the primary element 'i'.
    std::vector<Pointer_> cached_indptrs; 

    // This vector contains the cached index being pointed to by 'cached_indptrs'.
    // We store this here as it is more cache-friendly than doing a look-up to 'indices' every time.
    //
    // More specifically, if 'last_increasing = true', we define 'cached_indices[i] := indices[cached_indptrs[i]]'.
    // If 'cached_indptrs[i]' is out of range, the cached index is instead set to 'max_index'.
    // In short, 'cached_indices[i]' just contains the usual lower bound for 'last_request'.
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

    // What was the last requested index on the secondary dimension?
    Index_ last_request = 0;

    // Was the last requested index greater than its predecessor?
    bool last_increasing = true;

public:
    template<class PrimaryFunction_>
    SecondaryExtractionCache(IndexServer_ isrv, Index_ mi, Index_ length, PrimaryFunction_ to_primary) :
        indices(std::move(isrv)), max_index(mi), cached_indptrs(length), cached_indices(length) 
    {
        for (Index_ p = 0; p < length; ++p) {
            auto primary = to_primary(p);
            auto& curptr = cached_indptrs[p];
            curptr = indices.start_offset(primary);
            cached_indices[p] = (curptr == indices.end_offset(primary) ? max_index : *(indices.raw(primary) + curptr));
        }
    }

    auto size() const {
        return cached_indices.size();
    }

private:
    template<class Store_>
    void search_above(Index_ secondary, Index_ index_primary, Index_ primary, Store_& store) {
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
        auto endptr = indices.end_offset(primary);
        if (curptr == endptr) {
            curdex = max_index;
            return;
        }

        auto iraw = indices.raw(primary);
        auto inext = iraw + curptr;
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
        inext = std::lower_bound(inext + 1, iraw + endptr, secondary, [](Index_ a, Index_ b) -> bool { return a < b; });
        curptr = inext - iraw;
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
    template<class Store_>
    void search_below(Index_ secondary, Index_ index_primary, Index_ primary, Store_& store) {
        auto secondaryP1 = secondary + 1;
        auto& curdex = cached_indices[index_primary];
        if (curdex < secondaryP1) {
            return;
        }

        auto& curptr = cached_indptrs[index_primary];
        if (curdex == secondaryP1) {
            // Recall that 'curdex' is the reverse lower bound for 'last_request' (plus 1),
            // but 'curptr' points to the lower bound. This means that, if 'curptr' does not
            // already point to 'curdex', the former needs to be decremented so that 'store()' gets 
            // the right value. 'curptr' will only point to 'curdex' if 'curdex' is equal to
            // 'last_request + 1', as per the definition of 'cached_indices' above; thus,
            // the decrement should only occur if 'curdex != last_request + 1', which simplifies
            // to 'last_request != secondary' when we're inside this 'if' condition.
            curptr -= (last_request != secondary);
            store(primary, index_primary, curptr);
            return;
        }

        // Can't decrement any further.
        auto startptr = indices.start_offset(primary);
        if (curptr == startptr) {
            curdex = 0;
            return;
        }

        // Don't decrement 'curptr' yet; if 'curdex < secondary + 1', 'curptr' needs
        // to continue being the lower bound, while 'curdex' is the reverse lower bound.
        auto iraw = indices.raw(primary);
        auto inext = iraw + curptr - 1;
        curdex = *inext + 1;
        if (curdex < secondaryP1) {
            return;
        }

        if (curdex == secondaryP1) {
            // Decrement 'curptr' to match the definition of 'curdex', now that the 
            // lower bound and reverse lower bound are equal here.
            --curptr;
            store(primary, index_primary, curptr);
            return;
        }

        // Otherwise we need to search indices below the existing position.
        // Again, we use our own comparator to force casting and avoid
        // signed/unsigned comparisons.
        inext = std::lower_bound(iraw + startptr, inext, secondary, [](Index_ a, Index_ b) -> bool { return a < b; });
        curdex = *inext + 1;
        curptr = inext - iraw;

        if (curdex == secondaryP1) {
            // No need for decrement logic here, as both 'curdex' and 'curptr' are consistent right now.
            store(primary, index_primary, curptr);
            return;
        }

        if (curptr == startptr) {
            curdex = 0;
            return;
        }

        // Setting 'curdex' to the reverse lower bound again.
        --inext;
        curdex = *inext + 1;
        return;
    }

protected:
    template<class PrimaryFunction_, class Store_>
    void search_base(Index_ secondary, PrimaryFunction_ to_primary, Store_ store) {
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
                    search_above(secondary, p, to_primary(p), store);
                }

            } else {
                // Need to reset the meaning of 'cached_indices'.
                last_increasing = true;
                for (Index_ p = 0; p < primary_length; ++p) {
                    auto primary = to_primary(p);
                    auto curptr = cached_indptrs[p];
                    cached_indices[p] = (curptr == indices.end_offset(primary) ? max_index : *(indices.raw(primary) + curptr));
                    search_above(secondary, p, primary, store);
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
                    search_below(secondary, p, to_primary(p), store);
                }

            } else {
                // Need to reset the meaning of 'cached_indices'.
                last_increasing = false;
                for (Index_ p = 0; p < primary_length; ++p) {
                    auto primary = to_primary(p);
                    auto iraw = indices.raw(primary);
                    auto curptr = cached_indptrs[p];
                    if (curptr != indices.end_offset(primary) && *(iraw + curptr) == last_request) {
                        cached_indices[p] = last_request + 1;
                    } else if (curptr != indices.start_offset(primary)) {
                        cached_indices[p] = *(iraw + curptr - 1) + 1;
                    } else {
                        cached_indices[p] = 0;
                    }
                    search_below(secondary, p, primary, store);
                }
            }

            closest_cached_index = *(std::max_element(cached_indices.begin(), cached_indices.end()));
        }

        last_request = secondary; 
        return;
    }
};

// Subclasses for each selection typeblock_start.
template<typename Index_, class IndexServer_> 
struct FullSecondaryExtractionCache : public SecondaryExtractionCache<Index_, IndexServer_> {
    FullSecondaryExtractionCache(IndexServer_ isrv, Index_ mi, Index_ length) :
        SecondaryExtractionCache<Index_, IndexServer_>(std::move(isrv), mi, length, [](Index_ ip) -> Index_ { return ip; }) {}

    template<class Store_>
    void search(Index_ secondary, Store_ store) {
        this->search_base(secondary, [](Index_ ip) -> Index_ { return ip; }, std::move(store));
    }
};

template<typename Index_, class IndexServer_> 
struct BlockSecondaryExtractionCache : public SecondaryExtractionCache<Index_, IndexServer_> {
    BlockSecondaryExtractionCache(IndexServer_ isrv, Index_ mi, Index_ bs, Index_ bl) :
        SecondaryExtractionCache<Index_, IndexServer_>(std::move(isrv), mi, bl, [&](Index_ ip) -> Index_ { return bs + ip; }), block_start(bs) {}

    template<class Store_>
    void search(Index_ secondary, Store_ store) {
        this->search_base(secondary, [&](Index_ ip) -> Index_ { return ip + block_start; }, std::move(store));
    }

    Index_ block_start;
};

template<typename Index_, class IndexServer_> 
struct IndexSecondaryExtractionCache : public SecondaryExtractionCache<Index_, IndexServer_> {
    IndexSecondaryExtractionCache(IndexServer_ isrv, Index_ mi, std::vector<Index_> sub) :
        SecondaryExtractionCache<Index_, IndexServer_>(std::move(isrv), mi, sub.size(), [&](Index_ ip) -> Index_ { return sub[ip]; }), subset(std::move(sub)) {}

    template<class Store_>
    void search(Index_ secondary, Store_ store) {
        this->search_base(secondary, [&](Index_ ip) -> Index_ { return subset[ip]; }, std::move(store));
    }

    std::vector<Index_> subset;
};

}

}

#endif
