#ifndef TATAMI_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP
#define TATAMI_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP

#include "../base/Matrix.hpp"
#include "../utils/Index_to_container.hpp"

#include <vector>
#include <type_traits>
#include <algorithm>

namespace tatami {

namespace sparse_utils {

template<typename Index_, class IndexServer_> 
class SecondaryExtractionCache {
private:
    // Let 'my_indices[i] := my_indices_server.raw(i)' for some short-hand notation,
    // and assume that raw() produces some subscriptable object, e.g., a pointer.
    IndexServer_ my_indices_server;

    Index_ my_max_index;

    typedef typename IndexServer_::pointer_type Pointer;

    // The cached position of the pointer at each primary element.
    // Specifically, 'my_indices[i][my_cached_pointers[i]]' is the lower bound for 'my_last_request' in the primary element 'i'.
    std::vector<Pointer> my_cached_pointers; 

    // This vector contains the cached index being pointed to by 'my_cached_pointers'.
    // We store this here as it is more cache-friendly than doing a look-up to 'my_indices' every time.
    //
    // More specifically, if 'my_last_increasing = true', we define 'my_cached_indices[i] := my_indices[i][my_cached_pointers[i]]'.
    // If 'my_cached_pointers[i]' is out of range, the cached index is instead set to 'my_max_index'.
    // In short, 'my_cached_indices[i]' just contains the usual lower bound for 'my_last_request'.
    //
    // On the other hand, if 'my_last_increasing = false', we define 'my_cached_indices[i]' to be:
    // - 'my_last_request + 1', if 'my_indices[i][my_cached_pointers[i]] == my_last_request'.
    // - 0, if 'my_cached_pointers[i]' already refers to the start of the primary element 'i'.
    // - 'my_indices[i][my_cached_pointers[i] - 1] + 1' otherwise.
    // This can be considered the "reverse lower bound", i.e., the largest element not greater than 'my_last_request'.
    std::vector<Index_> my_cached_indices;

    // Closest value in 'my_cached_indices' to the 'my_last_request', to see whether we can short-circuit the iteration.
    // If 'my_last_increasing = true', this is the minimum of values in 'my_cached_indices'.
    // If 'my_last_increasing = false', this is the maximum of values in 'my_cached_indices'.
    Index_ my_closest_cached_index = 0;

    // What was the last requested index on the secondary dimension?
    Index_ my_last_request = 0;

    // Was the last requested index greater than its predecessor?
    bool my_last_increasing = true;

public:
    template<class PrimaryFunction_>
    SecondaryExtractionCache(IndexServer_ index_server, Index_ max_index, Index_ primary_length, PrimaryFunction_ to_primary) :
        my_indices_server(std::move(index_server)),
        my_max_index(max_index),
        my_cached_pointers(cast_Index_to_container_size<decltype(my_cached_pointers)>(primary_length)),
        my_cached_indices(cast_Index_to_container_size<decltype(my_cached_indices)>(primary_length))
    {
        for (Index_ p = 0; p < primary_length; ++p) {
            auto primary = to_primary(p);
            auto& curptr = my_cached_pointers[p];
            curptr = my_indices_server.start_offset(primary);
            my_cached_indices[p] = (curptr == my_indices_server.end_offset(primary) ? my_max_index : *(my_indices_server.raw(primary) + curptr));
        }
        if (primary_length) {
            my_closest_cached_index = *(std::min_element(my_cached_indices.begin(), my_cached_indices.end()));
        }
    }

    auto size() const {
        return my_cached_indices.size();
    }

private:
    template<class Store_>
    void search_above(Index_ secondary, Index_ index_primary, Index_ primary, Store_ store, bool& found) {
        // Skipping if the curdex (corresponding to curptr) is already higher
        // than secondary. So, we only need to do more work if the request is
        // greater than the stored index. This also catches cases where we're
        // at the end of the dimension, as curdex is set to my_max_index.
        auto& curdex = my_cached_indices[index_primary];
        if (curdex > secondary) {
            return;
        }

        auto& curptr = my_cached_pointers[index_primary];
        if (curdex == secondary) {
            store(primary, index_primary, my_cached_pointers[index_primary]);
            found = true;
            return;
        }

        // Having a peek at the index of the next non-zero element; maybe we're
        // lucky enough that the requested index is below this, as would be the
        // case for consecutive or near-consecutive accesses.
        ++curptr;
        auto endptr = my_indices_server.end_offset(primary);
        if (curptr == endptr) {
            curdex = my_max_index;
            return;
        }

        auto iraw = my_indices_server.raw(primary);
        auto inext = iraw + curptr;
        curdex = *inext;
        if (curdex > secondary) {
            return;
        }

        if (curdex == secondary) {
            store(primary, index_primary, curptr);
            found = true;
            return;
        }

        // Otherwise we need to search 'my_indices[primary]' above the existing
        // position. We do a quick increment to cut down the search space a bit
        // more. We also define our own comparator to avoid problems with
        // signed comparisons, depending on the types of IndexIt_ and Index_.
        inext = std::lower_bound(inext + 1, iraw + endptr, secondary, [](Index_ a, Index_ b) -> bool { return a < b; });
        curptr = inext - iraw;
        if (curptr == endptr) {
            curdex = my_max_index;
            return;
        }

        curdex = *inext;
        if (curdex > secondary) {
            return;
        }

        store(primary, index_primary, curptr);
        found = true;
        return;
    }

private:
    template<class Store_>
    void search_below(Index_ secondary, Index_ index_primary, Index_ primary, Store_ store, bool& found) {
        auto secondaryP1 = secondary + 1;
        auto& curdex = my_cached_indices[index_primary];
        if (curdex < secondaryP1) {
            return;
        }

        auto& curptr = my_cached_pointers[index_primary];
        if (curdex == secondaryP1) {
            // Recall that 'curdex' is the reverse lower bound for 'my_last_request' (plus 1),
            // but 'curptr' points to the lower bound. This means that, if 'curptr' does not
            // already point to 'curdex', the former needs to be decremented so that 'store()' gets 
            // the right value. 'curptr' will only point to 'curdex' if 'curdex' is equal to
            // 'my_last_request + 1', as per the definition of 'my_cached_indices' above; thus,
            // the decrement should only occur if 'curdex != my_last_request + 1', which simplifies
            // to 'my_last_request != secondary' when we're inside this 'if' condition.
            curptr -= (my_last_request != secondary);
            store(primary, index_primary, curptr);
            found = true;
            return;
        }

        // Can't decrement any further.
        auto startptr = my_indices_server.start_offset(primary);
        if (curptr == startptr) {
            curdex = 0;
            return;
        }

        // Don't decrement 'curptr' yet; if 'curdex < secondary + 1', 'curptr' needs
        // to continue being the lower bound, while 'curdex' is the reverse lower bound.
        auto iraw = my_indices_server.raw(primary);
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
            found = true;
            return;
        }

        // Otherwise we need to search 'my_indices[primary]' below the existing
        // position.  Again, we use our own comparator to force casting and
        // avoid signed/unsigned comparisons.
        inext = std::lower_bound(iraw + startptr, inext, secondary, [](Index_ a, Index_ b) -> bool { return a < b; });
        curdex = *inext + 1;
        curptr = inext - iraw;

        if (curdex == secondaryP1) {
            // No need for decrement logic here, as both 'curdex' and 'curptr' are consistent right now.
            store(primary, index_primary, curptr);
            found = true;
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

public:
    template<class PrimaryFunction_, class Store_>
    bool search(Index_ secondary, PrimaryFunction_ to_primary, Store_ store) {
        if (secondary > my_last_request || (my_last_increasing && secondary == my_last_request)) {
            bool found = false;

            if (my_last_increasing) {
                if (secondary < my_closest_cached_index) {
                    my_last_request = secondary;
                    return false; 
                }
                for (Index_ p = 0, plen = my_cached_indices.size(); p < plen; ++p) {
                    search_above(secondary, p, to_primary(p), store, found);
                }

            } else {
                // Need to reset the meaning of 'my_cached_indices'.
                my_last_increasing = true;
                for (Index_ p = 0, plen = my_cached_indices.size(); p < plen; ++p) {
                    auto primary = to_primary(p);
                    auto curptr = my_cached_pointers[p];
                    my_cached_indices[p] = (curptr == my_indices_server.end_offset(primary) ? my_max_index : *(my_indices_server.raw(primary) + curptr));
                    search_above(secondary, p, primary, store, found);
                }
            }

            if (found) {
                my_closest_cached_index = secondary;
            } else if (!my_cached_indices.empty()) {
                my_closest_cached_index = *(std::min_element(my_cached_indices.begin(), my_cached_indices.end()));
            }

        } else {
            bool found = false;

            if (!my_last_increasing) {
                if (secondary + 1 > my_closest_cached_index) {
                    my_last_request = secondary;
                    return false;
                }
                for (Index_ p = 0, plen = my_cached_indices.size(); p < plen; ++p) {
                    search_below(secondary, p, to_primary(p), store, found);
                }

            } else {
                // Need to reset the meaning of 'my_cached_indices'.
                my_last_increasing = false;
                for (Index_ p = 0, plen = my_cached_indices.size(); p < plen; ++p) {
                    auto primary = to_primary(p);
                    auto iraw = my_indices_server.raw(primary);
                    auto curptr = my_cached_pointers[p];

                    // Casting to Index_, as there's no guarantee that the
                    // stored index type can represent my_last_request.
                    if (curptr != my_indices_server.end_offset(primary) && static_cast<Index_>(*(iraw + curptr)) == my_last_request) {
                        my_cached_indices[p] = my_last_request + 1;
                    } else if (curptr != my_indices_server.start_offset(primary)) {
                        my_cached_indices[p] = *(iraw + curptr - 1) + 1;
                    } else {
                        my_cached_indices[p] = 0;
                    }
                    search_below(secondary, p, primary, store, found);
                }
            }

            if (found) {
                my_closest_cached_index = secondary + 1;
            } else if (!my_cached_indices.empty()) {
                my_closest_cached_index = *(std::max_element(my_cached_indices.begin(), my_cached_indices.end()));
            }
        }

        my_last_request = secondary; 
        return true;
    }
};

// Wrapper classes for each selection type.
template<typename Index_, class IndexServer_> 
class FullSecondaryExtractionCache {
public:
    FullSecondaryExtractionCache(IndexServer_ index_server, Index_ max_index, Index_ primary_length) :
        my_cache(std::move(index_server), max_index, primary_length, [](Index_ ip) -> Index_ { return ip; }) {}

    template<class Store_>
    bool search(Index_ secondary, Store_ store) {
        return my_cache.search(secondary, [](Index_ ip) -> Index_ { return ip; }, std::move(store));
    }

    auto size() const {
        return my_cache.size();
    }

private:
    SecondaryExtractionCache<Index_, IndexServer_> my_cache;
};

template<typename Index_, class IndexServer_> 
class BlockSecondaryExtractionCache {
public:
    BlockSecondaryExtractionCache(IndexServer_ index_server, Index_ max_index, Index_ block_start, Index_ block_length) :
        my_cache(std::move(index_server), max_index, block_length, Helper(block_start)), my_block_start(block_start) {}

    template<class Store_>
    bool search(Index_ secondary, Store_ store) {
        return my_cache.search(secondary, Helper(my_block_start), std::move(store));
    }

    auto size() const {
        return my_cache.size();
    }

private:
    SecondaryExtractionCache<Index_, IndexServer_> my_cache;
    Index_ my_block_start;

    // Just to avoid rewriting the lambda all the time.
    struct Helper {
        Helper(Index_ s) : shift(s) {}
        Index_ shift;
        Index_ operator()(Index_ ip) const {
            return ip + shift;
        }
    };
};

template<typename Index_, class IndexServer_> 
class IndexSecondaryExtractionCache {
public:
    IndexSecondaryExtractionCache(IndexServer_ index_server, Index_ max_index, VectorPtr<Index_> indices_ptr) :
        my_cache(std::move(index_server), max_index, indices_ptr->size(), Helper(*indices_ptr)),
        my_indices_ptr(std::move(indices_ptr)) 
    {}

    template<class Store_>
    bool search(Index_ secondary, Store_ store) {
        return my_cache.search(secondary, Helper(*my_indices_ptr), std::move(store));
    }

    auto size() const {
        return my_cache.size();
    }

private:
    SecondaryExtractionCache<Index_, IndexServer_> my_cache;
    VectorPtr<Index_> my_indices_ptr;

    // Just to avoid rewriting the lambda all the time.
    struct Helper {
        Helper(const std::vector<Index_>& s) : subset(s) {}
        const std::vector<Index_>& subset;
        Index_ operator()(Index_ ip) const {
            return subset[ip];
        }
    };
};

}

}

#endif
