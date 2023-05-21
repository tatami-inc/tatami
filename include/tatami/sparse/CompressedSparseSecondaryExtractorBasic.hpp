#ifndef TATAMI_COMPRESSED_SPARSE_SECONDARY_EXTRACTOR_BASIC_HPP
#define TATAMI_COMPRESSED_SPARSE_SECONDARY_EXTRACTOR_BASIC_HPP

#include <vector>

namespace tatami {

template<class Storage_>
using Stored = typename std::remove_reference<decltype(std::declval<Storage_>()[0])>::type;

template<typename Index_, typename StoredIndex_, typename CustomPointer_, class CustomPointerModifier_> 
struct CompressedSparseSecondaryExtractorBasic {
private:
    StoredIndex_ max_index;

    // The current position of the pointer at each primary element.
    std::vector<CustomPointer_> current_indptrs; 

    // Whether to move forward or back.
    bool lower_bound = true;

    // If 'lower_bound = true', this vector contains the current index being pointed to, i.e., 'current_indices[i] := indices[current_indptrs[i]]'.
    // We store this here as it is more cache-friendly than doing a look-up to 'indices'.
    // If 'current_indptrs[i]' is out of range, the current index is instead set to the maximum index (e.g., max rows for CSC matrices).
    //
    // If 'lower_bound = false', this vector instead contains:
    // - 'indices[current_indptrs[i] - 1]', if 'current_indptrs[i]' does not lie at the start of the primary dimension element.
    // - otherwise, 'decrement_fail', i.e., -1 or its unsigned equivalent.
    std::vector<StoredIndex_> current_indices;

    // Closest value in 'current_indices' to the 'last_request', to see whether we can short-circuit the iteration.
    StoredIndex_ closest_current_index;

    // What was the last requested index on the secondary dimension?
    StoredIndex_ last_request = 0;

public:
    CompressedSparseSecondaryExtractorBasic() = default;

    template<class IndexStorage_, class PointerStorage_>
    CompressedSparseSecondaryExtractorBasic(StoredIndex_ mi, const IndexStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
        max_index(mi), current_indices(length), current_indptrs(idp.begin() + start, idp.begin() + start + length)
    {
        /* Here, the general idea is to store a local copy of the actual
         * row indices (for CSC matrices; column indices, for CSR matrices)
         * so that we don't have to keep on doing cache-unfriendly look-ups
         * for the indices based on the pointers that we do have. This assumes
         * that the density is so low that updates to the local indices are
         * rare relative to the number of comparisons to those same indices.
         */
        auto idpIt = idp.begin() + start;
        for (Index_ i = 0; i < length; ++i, ++idpIt) {
            current_indices[i] = (*idpIt < *(idpIt + 1) ? idx[*idpIt] : max_index);
        }
        closest_current_index = (length ? *std::min_element(current_indices.begin(), current_indices.end()) : max_index);
        return;
    } 

    template<class IndexStorage_, class PointerStorage_>
    CompressedSparseSecondaryExtractorBasic(StoredIndex_ mi, const IndexStorage_& idx, const PointerStorage_& idp) :
        CompressedSparseSecondaryExtractorBasic(mi, idx, idp, static_cast<Index_>(0), static_cast<Index_>(idp.size() - 1)) {}

    template<class IndexStorage_, class PointerStorage_>
    CompressedSparseSecondaryExtractorBasic(StoredIndex_ mi, const IndexStorage_& idx, const PointerStorage_& idp, const Index_* subset, Index_ length) :
        max_index(mi), current_indices(length), current_indptrs(length)
    {
        for (Index_ i0 = 0; i0 < length; ++i0) {
            auto i = subset[i0];
            current_indptrs[i0] = idp[i];
            current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
        }
        closest_current_index = (length ? *std::min_element(current_indices.begin(), current_indices.end()) : max_index);
        return;
    }

public:
    template<bool reset_index_, class IndexStorage_, class PointerStorage_, class StoreFunction_, class SkipFunction_>
    void search_above_or_equal(
        StoredIndex_ secondary, 
        Index_ index_primary,
        Index_ primary,
        const IndexStorage_& indices,
        const PointerStorage_& indptrs,
        StoreFunction_& store,
        SkipFunction_& skip
    ) {
        auto& curdex = current_indices[index_primary];

        // Templated check to avoid having to incur the cost of hitting this,
        // given that index resets are irrelevant for the common case of
        // incremented iteration over consecutive secondary elements.
        if constexpr(reset_index_) {
            auto limit = indptrs[primary + 1];
            auto curptr = current_indptrs[index_primary];
            auto raw_ptr = CustomPointerModifier_::get(curptr);
            if (raw_ptr != limit) {
                curdex = indices[raw_ptr];
            } else {
                curdex = max_index;
            }
        }

        // Skipping if the curdex (corresponding to curptr) is already higher
        // than secondary.  So, we only need to do more work if the request is
        // greater than the stored index.  This also catches cases where we're
        // at the end of the dimension, as curdex is set to max_index.
        if (curdex > secondary) {
            skip(primary);
            return;
        }

        auto& curptr = current_indptrs[index_primary];
        if (curdex == secondary) {
            store(primary, curptr);
            return;
        }

        auto limit = indptrs[primary + 1];

        // Special case if the requested index is at the end of the matrix, in
        // which case we can just jump there directly rather than doing an
        // unnecessary binary search. Note that 'limit - 1' is always available
        // as this dimension element should be non-empty if secondary > curdex.
        if (secondary + 1 == max_index) {
            if (indices[limit - 1] == secondary) {
                CustomPointerModifier_::set(curptr, limit); // don't set directly to 'limit - 1' as this won't be the start of the run in semi-compressed mode.
                CustomPointerModifier_::decrement(curptr, indices, indptrs[primary]);
                curdex = secondary;
                store(primary, curptr);
            } else {
                CustomPointerModifier_::set(curptr, limit);
                curdex = max_index;
                skip(primary);
            }
            return;
        }

        // Having a peek at the index of the next non-zero element; maybe we're
        // lucky enough that the requested index is below this, as would be the
        // case for consecutive or near-consecutive accesses.
        CustomPointerModifier_::increment(curptr, indices, limit);
        auto raw_ptr = CustomPointerModifier_::get(curptr);
        if (raw_ptr == limit) {
            curdex = max_index;
            skip(primary);
            return;
        }

        curdex = indices[raw_ptr];
        if (curdex > secondary) {
            skip(primary);
            return;
        }

        if (curdex == secondary) {
            store(primary, curptr);
            return;
        }

        // Otherwise we need to search indices above the existing position.
        // We do a quick increment to cut down the search space; don't
        // need to pay the cost of using increment() here, as the lower
        // bound search is going to be faster than any increment.
        ++raw_ptr;
        Stored<PointerStorage_> next_ptr = std::lower_bound(indices.begin() + raw_ptr, indices.begin() + limit, secondary) - indices.begin();
        CustomPointerModifier_::set(curptr, next_ptr);

        if (next_ptr == limit) {
            curdex = max_index;
            skip(primary);
            return;
        }

        curdex = indices[next_ptr];
        if (curdex > secondary) {
            skip(primary);
            return;
        }

        store(primary, curptr);
        return;
    }

public:
    static constexpr StoredIndex_ decrement_fail = -1;

    template<bool check_index_, class IndexStorage_, class PointerStorage_, class StoreFunction_, class SkipFunction_>
    void search_below(
        StoredIndex_ secondary, 
        Index_ index_primary, 
        Index_ primary,
        const IndexStorage_& indices,
        const PointerStorage_& indptrs,
        StoreFunction_& store,
        SkipFunction_& skip
    ) {
        // In the context of this function, there's no need to check if
        // 'indices[curptr] == secondary'; we only enter this function if
        // 'last_request > secondary', and we already know that
        // 'indices[curptr] >= last_request > secondary'.

        // Checking if the next-lowest index is below the requested
        // 'secondary', at which point we can just quit without a
        // cache-unfriendly lookup to 'indices'. This is valid as we know that
        // we are not at the start of the dimension at this point. 
        auto& curdex = current_indices[index_primary];
        if constexpr(check_index_) {
            if (curdex < secondary || curdex == decrement_fail) { // need to check decrement_fail separately, in case StoredIndex_ is unsigned.
                skip(primary);
                return;
            }
        }

        curdex = decrement_fail;
        auto& curptr = current_indptrs[index_primary];

        // Can't decrement anymore, in which case we quit. 
        auto lower_limit = indptrs[primary];
        auto raw_ptr = CustomPointerModifier_::get(curptr);
        if (raw_ptr == lower_limit) {
            skip(primary);
            return;
        }

        // Special case if the requested index is at the end of the matrix, in
        // which case we can just jump there directly rather than doing an
        // unnecessary binary search.
        if (secondary == 0) {
            auto first_ptr = lower_limit;
            CustomPointerModifier_::set(curptr, first_ptr);

            if (indices[first_ptr] == secondary) {
                store(primary, curptr);
            } else {
                skip(primary);
            }
            return;
        }

        // Having a peek at the index of the next non-zero element and seeing
        // whether we can stop searching, as would be the case for consecutive
        // or near-consecutive accesses.
        --raw_ptr;
        auto candidate = indices[raw_ptr];
        if (candidate < secondary) {
            curdex = candidate;
            skip(primary);
            return;
        }

        if (candidate == secondary) {
            CustomPointerModifier_::decrement(curptr, indices, lower_limit);
            if (raw_ptr != lower_limit) {
                curdex = indices[raw_ptr - 1]; // cheap decrement to inspect the next-lowest element.
            }
            store(primary, curptr);
            return;
        }

        // Otherwise, searching indices below the (just-decremented) position.
        Stored<PointerStorage_> next_ptr = std::lower_bound(indices.begin() + lower_limit, indices.begin() + raw_ptr, secondary) - indices.begin();
        CustomPointerModifier_::set(curptr, next_ptr);

        if (next_ptr == indptrs[primary + 1]) {
            skip(primary);
            return;
        }

        if (indices[next_ptr] == secondary) {
            if (next_ptr != lower_limit) {
                curdex = indices[next_ptr - 1]; // cheap decrement to inspect the next-lowest element.
            }
            store(primary, curptr);
            return;
        }

        if (next_ptr != lower_limit) {
            curdex = indices[next_ptr - 1]; // cheap decrement to inspect the next-lowest element.
        }
        skip(primary);
        return;
    }

public:
    template<class IndexStorage_, class PointerStorage_, class PrimaryFunction_, class StoreFunction_, class SkipFunction_>
    bool search(
        StoredIndex_ secondary,
        Index_ primary_length,
        PrimaryFunction_ to_primary, 
        const IndexStorage_& indices,
        const PointerStorage_& indptrs,
        StoreFunction_ store,
        SkipFunction_ skip
    ) {
        if (secondary >= last_request) {
            if (lower_bound) {
                if (secondary < closest_current_index) {
                    return false;
                }
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_above_or_equal<false>(secondary, p, to_primary(p), indices, indptrs, store, skip);
                }

            } else {
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_above_or_equal<true>(secondary, p, to_primary(p), indices, indptrs, store, skip);
                }
                lower_bound = true;
            }

            // Doing a single min_element call is faster than successive min() calls
            // within search_above_or_equal() on GCC. Who knows why.
            if (primary_length) {
                closest_current_index = *std::min_element(current_indices.begin(), current_indices.end());
            }

        } else {
            if (!lower_bound) {
                if (closest_current_index == decrement_fail || secondary > closest_current_index) {
                    return false;
                }
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_below<true>(secondary, p, to_primary(p), indices, indptrs, store, skip);
                }

            } else {
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_below<false>(secondary, p, to_primary(p), indices, indptrs, store, skip);
                }
                lower_bound = false;
            }

            if constexpr(std::is_signed<StoredIndex_>::value) {
                if (primary_length) {
                    closest_current_index = *std::max_element(current_indices.begin(), current_indices.end());
                }
            } else {
                closest_current_index = decrement_fail;
                for (auto x : current_indices) {
                    if (x != decrement_fail && (x > closest_current_index || closest_current_index == decrement_fail)) {
                        closest_current_index = x;
                    }
                }
            }
        }

        last_request = secondary;
        return true;
    }
};

}

#endif
