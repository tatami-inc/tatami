#ifndef TATAMI_FRAGMENTED_SPARSE_SECONDARY_EXTRACTOR_BASIC_HPP
#define TATAMI_FRAGMENTED_SPARSE_SECONDARY_EXTRACTOR_BASIC_HPP

#include <vector>

namespace tatami {

template<class Storage_>
using Stored = typename std::remove_reference<decltype(std::declval<Storage_>()[0])>::type;

// Clone of the CompressedSparseSecondaryExtractorBasic, but the indices are
// now in their own vectors rather than being in a single monolithic vector.
// This is not easily handled by templating of the Compressed class, hence the
// need for a copied class entirely.
template<typename Index_, typename StoredIndex_, typename CustomPointer_, class CustomPointerModifier_> 
struct FragmentedSparseSecondaryExtractorBasic {
private:
    StoredIndex_ max_index;

    std::vector<CustomPointer_> current_indptrs; 

    bool lower_bound = true;

    std::vector<StoredIndex_> current_indices;

    StoredIndex_ closest_current_index;

    StoredIndex_ last_request = 0;

public:
    FragmentedSparseSecondaryExtractorBasic() = default;

    FragmentedSparseSecondaryExtractorBasic(StoredIndex_ mi, Index_ length) : max_index(mi), current_indices(length), current_indptrs(length) {}

public:
    template<bool reset_index_, class IndexVectorStorage_, class StoreFunction_, class SkipFunction_>
    void search_above_or_equal(
        StoredIndex_ secondary, 
        Index_ index_primary,
        Index_ primary,
        const IndexVectorStorage_& indices,
        StoreFunction_& store,
        SkipFunction_& skip
    ) {
        auto& curdex = current_indices[index_primary];

        // Templated check to avoid having to incur the cost of hitting this,
        // given that index resets are irrelevant for the common case of
        // incremented iteration over consecutive secondary elements.
        if constexpr(reset_index_) {
            auto limit = indices[primary].size();
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

        auto limit = indices[primary].size();

        // Special case if the requested index is at the end of the matrix, in
        // which case we can just jump there directly rather than doing an
        // unnecessary binary search. Note that 'limit - 1' is always available
        // as this dimension element should be non-empty if secondary > curdex.
        if (secondary + 1 == max_index) {
            if (indices[limit - 1] == secondary) {
                CustomPointerModifier_::set(curptr, limit); // don't set directly to 'limit - 1' as this won't be the start of the run in semi-compressed mode.
                CustomPointerModifier_::decrement(curptr, indices, 0);
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

    template<bool check_index_, class IndexVectorStorage_, class StoreFunction_, class SkipFunction_>
    void search_below(
        StoredIndex_ secondary, 
        Index_ index_primary, 
        Index_ primary,
        const IndexVectorStorage_& indices,
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
        auto raw_ptr = CustomPointerModifier_::get(curptr);
        if (raw_ptr == 0) {
            skip(primary);
            return;
        }

        // Special case if the requested index is at the end of the matrix, in
        // which case we can just jump there directly rather than doing an
        // unnecessary binary search.
        if (secondary == 0) {
            CustomPointerModifier_::set(curptr, 0);

            if (indices[0] == secondary) {
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
            CustomPointerModifier_::decrement(curptr, indices, 0);
            if (raw_ptr != 0) {
                curdex = indices[raw_ptr - 1]; // cheap decrement to inspect the next-lowest element.
            }
            store(primary, curptr);
            return;
        }

        // Otherwise, searching indices below the (just-decremented) position.
        Stored<PointerStorage_> next_ptr = std::lower_bound(indices.begin(), indices.begin() + raw_ptr, secondary) - indices.begin();
        CustomPointerModifier_::set(curptr, next_ptr);
        if (next_ptr == raw_ptr) {
            skip(primary);
            return;
        }

        if (indices[next_ptr] == secondary) {
            if (next_ptr != 0) {
                curdex = indices[next_ptr - 1]; // cheap decrement to inspect the next-lowest element.
            }
            store(primary, curptr);
            return;
        }

        if (next_ptr != 0) {
            curdex = indices[next_ptr - 1]; // cheap decrement to inspect the next-lowest element.
        }
        skip(primary);
        return;
    }

public:
    template<class IndexVectorStorage_, class PrimaryFunction_, class StoreFunction_, class SkipFunction_>
    bool search(
        StoredIndex_ secondary,
        Index_ primary_length,
        PrimaryFunction_ to_primary, 
        const IndexVectorStorage_& indices,
        StoreFunction_ store,
        SkipFunction_ skip
    ) {
        if (secondary >= last_request) {
            if (lower_bound) {
                if (secondary < closest_current_index) {
                    return false;
                }
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_above_or_equal<false>(secondary, p, to_primary(p), indices, store, skip);
                }

            } else {
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_above_or_equal<true>(secondary, p, to_primary(p), indices, store, skip);
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
                    search_below<true>(secondary, p, to_primary(p), indices, store, skip);
                }

            } else {
                for (Index_ p = 0; p < primary_length; ++p) {
                    search_below<false>(secondary, p, to_primary(p), indices, store, skip);
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
