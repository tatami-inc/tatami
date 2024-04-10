#ifndef TATAMI_SPARSE_PRIMARY_EXTRACTION_HPP
#define TATAMI_SPARSE_PRIMARY_EXTRACTION_HPP

#include "../utils/has_data.hpp"

#include <algorithm>

namespace tatami {

namespace sparse_utils {

template<class Storage_, typename Pointer_, typename Data_>
const Data_* extract_primary_vector(const Storage_& input, Pointer_ offset, Pointer_ delta, Data_* buffer) {
    if constexpr(has_data<Data_, Storage_>::value) {
        return input.data() + offset;
    } else {
        auto it = input.begin() + offset;
        std::copy(it, it + delta, buffer);
        return buffer;
    }
}

template<class IndexIt_, typename Index_>
void refine_primary_limits(IndexIt_& indices_start, IndexIt_& indices_end, Index_ extent, Index_ smallest, Index_ largest_plus_one) {
    if (smallest) {
        // Using custom comparator to ensure that we cast to Index_ for signedness-safe comparisons.
        indices_start = std::lower_bound(indices_start, indices_end, smallest, [](Index_ a, Index_ b) -> bool { return a < b; });
    }

    if (largest_plus_one != extent) {
        indices_end = std::lower_bound(indices_start, indices_end, largest_plus_one, [](Index_ a, Index_ b) -> bool { return a < b; });
    }
}

template<class IndexIt_, typename Index_>
void refine_primary_block_limits(IndexIt_& indices_start, IndexIt_& indices_end, Index_ extent, Index_ block_start, Index_ block_length) {
    refine_primary_limits(indices_start, indices_end, extent, block_start, block_start + block_length);
}

template<typename Index_>
struct RetrievePrimarySubsetDense {
    RetrievePrimarySubsetDense(const std::vector<Index_>& subset, Index_ extent) : extent(extent) {
        if (!subset.empty()) {
            offset = subset.front();
            lastp1 = subset.back() + 1;
            size_t alloc = lastp1 - offset;
            present.resize(alloc);

            // Starting off at 1 to ensure that 0 is still a marker for
            // absence. It should be fine as subset.size() should fit inside
            // Index_ (otherwise nrow()/ncol() would give the wrong answer).
            Index_ counter = 1; 

            for (auto s : subset) {
                present[s - offset] = counter;
                ++counter;
            }
        }
    }

    template<class IndexIt_, class Store_>
    void populate(IndexIt_ indices_start, IndexIt_ indices_end, Store_ store) const {
        if (present.empty()) {
            return;
        }

        // Limiting the iteration to its boundaries based on the first and last subset index.
        auto original_start = indices_start;
        refine_primary_limits(indices_start, indices_end, extent, offset, lastp1);

        size_t counter = indices_start - original_start;
        for (; indices_start != indices_end; ++indices_start, ++counter) {
            auto ix = *indices_start;
            auto shift = present[ix - offset];
            if (shift) {
                store(shift - 1, counter);
            }
        }
    }

    Index_ extent;
    std::vector<Index_> present;
    Index_ offset = 0;
    Index_ lastp1 = 0;
};

template<typename Index_>
struct RetrievePrimarySubsetSparse {
    RetrievePrimarySubsetSparse(const std::vector<Index_>& subset, Index_ extent) : extent(extent) {
        if (!subset.empty()) {
            offset = subset.front();
            lastp1 = subset.back() + 1;
            size_t alloc = lastp1 - offset;
            present.resize(alloc);

            // Unlike the dense case, this is a simple present/absent signal,
            // as we don't need to map each structural non-zero back onto its 
            // corresponding location on a dense vector.
            for (auto s : subset) {
                present[s - offset] = 1;
            }
        }
    }

    template<class IndexIt_, class Store_>
    void populate(IndexIt_ indices_start, IndexIt_ indices_end, Store_ store) const {
        if (present.empty()) {
            return;
        }

        // Limiting the iteration to its boundaries based on the first and last subset index.
        auto original_start = indices_start;
        refine_primary_limits(indices_start, indices_end, extent, offset, lastp1);

        size_t counter = indices_start - original_start;
        for (; indices_start != indices_end; ++indices_start, ++counter) {
            auto ix = *indices_start;
            if (present[ix - offset]) {
                store(counter, ix);
            }
        }
    }

    Index_ extent;
    std::vector<unsigned char> present;
    Index_ offset = 0;
    Index_ lastp1 = 0;
};

}

}

#endif
