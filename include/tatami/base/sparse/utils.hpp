#ifndef TATAMI_SPARSE_UTILS_HPP
#define TATAMI_SPARSE_UTILS_HPP

#include <vector>

namespace tatami {

namespace sparse {

template<class Storage_>
using Stored = typename std::remove_reference<decltype(std::declval<Storage_>()[0])>::type;

template<typename StoredIndex_, class CustomPointerModifier_> 
struct SecondaryExtractionWorkspaceBase {
private:
    StoredIndex_ max_index;

public:
    SecondaryExtractionWorkspaceBase(StoredIndex_ mi) : max_index(mi) {}

    SecondaryExtractionWorkspaceBase() = default;

public:
    template<typename Index_, class IndexStorage_, class PointerStorage_, class CustomPointer_>
    void search_above(StoredIndex_ secondary, Index_ primary, const IndexStorage_& indices, const PointerStorage_& indptrs, StoredIndex_& curdex, CustomPointer_& curptr) const {
        // Skipping if the curdex (corresponding to curptr) is already higher
        // than secondary.  So, we only need to do more work if the request is
        // greater than the stored index.  This also catches cases where we're
        // at the end of the dimension, as curdex is set to max_index.
        if (secondary <= curdex) {
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
            } else {
                CustomPointerModifier_::set(curptr, limit);
                curdex = max_index;
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
            return;
        }

        auto candidate = indices[raw_ptr];
        if (candidate >= secondary) {
            curdex = candidate;
            return;
        }

        // Otherwise we need to search indices above the existing position.
        ++raw_ptr;
        Stored<PointerStorage_> next_ptr = std::lower_bound(indices.begin() + raw_ptr, indices.begin() + limit, secondary) - indices.begin();
        CustomPointerModifier_::set(curptr, next_ptr);
        curdex = (next_ptr != limit ? indices[next_ptr] : max_index);
    }

    template<typename Index_, class IndexStorage_, class PointerStorage_, class CustomPointer_>
    void search_below(StoredIndex_ secondary, Index_ primary, const IndexStorage_& indices, const PointerStorage_& indptrs, StoredIndex_& curdex, CustomPointer_& curptr) const {
        if (secondary == curdex) {
            return;
        }

        auto lower_limit = indptrs[primary];
        if (CustomPointerModifier_::get(curptr) == lower_limit) {
            return;
        }

        // Special case if the requested index is at the end of the matrix, in
        // which case we can just jump there directly rather than doing an
        // unnecessary binary search.
        if (secondary == 0) {
            auto first_ptr = indptrs[primary];
            CustomPointerModifier_::set(curptr, first_ptr);
            curdex = indices[first_ptr];
            return;
        }

        // Having a peek at the index of the next non-zero element and
        // seeing whether we stop searching, as would be the case for
        // consecutive or near-consecutive accesses.
        CustomPointerModifier_::decrement(curptr, indices, lower_limit);
        auto raw_ptr = CustomPointerModifier_::get(curptr);
        auto candidate = indices[raw_ptr];
        if (candidate < secondary) {
            CustomPointerModifier_::increment(curptr, indices, indptrs[primary + 1]); // oops, went too far; undo the change.
            return;
        }

        curdex = candidate; 
        if (candidate == secondary) {
            return;
        }

        // Otherwise, searching indices below the existing position.
        Stored<PointerStorage_> next_ptr =  std::lower_bound(indices.begin() + lower_limit, indices.begin() + raw_ptr, secondary) - indices.begin();
        CustomPointerModifier_::set(curptr, next_ptr);
        curdex = (next_ptr != indptrs[primary + 1] ? indices[next_ptr] : max_index);
    }
};

template<typename Index_, typename StoredIndex_, class CustomPointer_, class CustomPointerModifier_>
struct SimpleSecondaryExtractionWorkspace {
public:
    SimpleSecondaryExtractionWorkspace() = default;

    template<class IndexStorage_, class PointerStorage_>
    SimpleSecondaryExtractionWorkspace(StoredIndex_ max_index, const IndexStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
        base(max_index), current_indices(length), current_indptrs(idp.begin() + start, idp.begin() + start + length)
    {
        /* Here, the general idea is to store a local copy of the actual
         * row indices (for CSC matrices; column indices, for CSR matrices)
         * so that we don't have to keep on doing cache-unfriendly look-ups
         * for the indices based on the pointers that we do have. This assumes
         * that the density is so low that updates to the local indices are
         * rare relative to the number of comparisons to those same indices.
         * Check out the `secondary_dimension()` function for how this is used.
         */
        auto idpIt = idp.begin() + start;
        for (Index_ i = 0; i < length; ++i, ++idpIt) {
            current_indices[i] = (*idpIt < *(idpIt + 1) ? idx[*idpIt] : max_index);
        }
        return;
    } 

    template<class IndexStorage_, class PointerStorage_>
    SimpleSecondaryExtractionWorkspace(StoredIndex_ max_index, const IndexStorage_& idx, const PointerStorage_& idp) :
        SimpleSecondaryExtractionWorkspace(max_index, idx, idp, static_cast<Index_>(0), static_cast<Index_>(idp.size() - 1)) {}

    template<class IndexStorage_, class PointerStorage_>
    SimpleSecondaryExtractionWorkspace(StoredIndex_ max_index, const IndexStorage_& idx, const PointerStorage_& idp, const Index_* subset, Index_ length) :
        base(max_index), current_indices(length), current_indptrs(length)
    {
        for (Index_ i0 = 0; i0 < length; ++i0) {
            auto i = subset[i0];
            current_indptrs[i0] = idp[i];
            current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
        }
        return;
    }

public:
    SecondaryExtractionWorkspaceBase<StoredIndex_, CustomPointerModifier_> base;

    // The current position of the pointer at each primary element.
    std::vector<CustomPointer_> current_indptrs; 

    // The current index being pointed to, i.e., current_indices[0] <= indices[current_ptrs[0]]. 
    // If current_ptrs[0] is out of range, the current index is instead set to the maximum index
    // (e.g., max rows for CSC matrices).
    std::vector<StoredIndex_> current_indices;

public:
    template<class IndexStorage_, class PointerStorage_>
    StoredIndex_ search_above(StoredIndex_ secondary, Index_ primary, Index_ index_primary, const IndexStorage_& indices, const PointerStorage_& indptrs) {
        auto& curdex = current_indices[index_primary];
        base.search_above(secondary, primary, indices, indptrs, curdex, current_indptrs[index_primary]);
        return curdex;
    }

    template<class IndexStorage_, class PointerStorage_>
    StoredIndex_ search_below(StoredIndex_ secondary, Index_ primary, Index_ index_primary, const IndexStorage_& indices, const PointerStorage_& indptrs) {
        auto& curdex = current_indices[index_primary];
        base.search_below(secondary, primary, indices, indptrs, curdex, current_indptrs[index_primary]);
        return curdex;
    }
};

}

}

#endif
