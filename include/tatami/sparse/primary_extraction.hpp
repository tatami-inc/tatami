#ifndef TATAMI_SPARSE_PRIMARY_EXTRACTION_HPP
#define TATAMI_SPARSE_PRIMARY_EXTRACTION_HPP

#include "../utils/has_data.hpp"
#include "../utils/Index_to_container.hpp"

#include <algorithm>

namespace tatami {

namespace sparse_utils {

template<class Storage_, typename Offset_, typename Data_>
const Data_* extract_primary_vector(const Storage_& input, const Offset_ offset, const Offset_ delta, Data_* const buffer) {
    if constexpr(has_data<Data_, Storage_>::value) {
        return input.data() + offset;
    } else {
        const auto it = input.begin() + offset;
        std::copy_n(it, delta, buffer);
        return buffer;
    }
}

template<class IndexIt_, typename Index_>
void refine_primary_limits(IndexIt_& indices_start, IndexIt_& indices_end, const Index_ extent, const Index_ smallest, const Index_ largest_plus_one) {
    if (smallest) {
        // Using custom comparator to ensure that we cast to Index_ for signedness-safe comparisons.
        indices_start = std::lower_bound(indices_start, indices_end, smallest, [](Index_ a, Index_ b) -> bool { return a < b; });
    }

    if (largest_plus_one != extent) {
        indices_end = std::lower_bound(indices_start, indices_end, largest_plus_one, [](Index_ a, Index_ b) -> bool { return a < b; });
    }
}

template<class IndexIt_, typename Index_>
void refine_primary_block_limits(IndexIt_& indices_start, IndexIt_& indices_end, const Index_ extent, const Index_ block_start, const Index_ block_length) {
    refine_primary_limits(indices_start, indices_end, extent, block_start, block_start + block_length);
}

template<typename Index_>
class RetrievePrimarySubsetDense {
public:
    RetrievePrimarySubsetDense(const std::vector<Index_>& subset, const Index_ extent) : my_extent(extent) {
        if (!subset.empty()) {
            my_offset = subset.front();
            my_lastp1 = subset.back() + 1; // +1 is safe as Index_ can hold the dimension extent, which is greater than all subset indices.
            resize_container_to_Index_size(my_present, my_lastp1 - my_offset);

            Index_ counter = 0; 
            for (const auto s : subset) {
                // Starting off at 1 to ensure that 0 is still a marker for
                // absence. It should be fine as subset.size() should fit inside
                // Index_ (otherwise nrow()/ncol() would give the wrong answer).
                ++counter;
                my_present[s - my_offset] = counter;
            }
        }
    }

private:
    Index_ my_extent;
    std::vector<Index_> my_present;
    Index_ my_offset = 0;
    Index_ my_lastp1 = 0;

public:
    template<class IndexIt_, class Store_>
    void populate(IndexIt_ indices_start, IndexIt_ indices_end, const Store_ store) const {
        if (my_present.empty()) {
            return;
        }

        // Limiting the iteration to its boundaries based on the first and last subset index.
        const auto original_start = indices_start;
        refine_primary_limits(indices_start, indices_end, my_extent, my_offset, my_lastp1);

        for (; indices_start != indices_end; ++indices_start) {
            const auto shift = my_present[*indices_start - my_offset];
            if (shift) {
                store(shift - 1, indices_start - original_start);
            }
        }
    }
};

template<typename Index_>
class RetrievePrimarySubsetSparse {
public:
    RetrievePrimarySubsetSparse(const std::vector<Index_>& subset, const Index_ extent) : my_extent(extent) {
        if (!subset.empty()) {
            my_offset = subset.front();
            my_lastp1 = subset.back() + 1; // +1 is safe as Index_ can hold the dimension extent, which is greater than all subset indices.
            resize_container_to_Index_size(my_present, my_lastp1 - my_offset);

            // Unlike the dense case, this is a simple present/absent signal,
            // as we don't need to map each structural non-zero back onto its 
            // corresponding location on a dense vector.
            for (const auto s : subset) {
                my_present[s - my_offset] = 1;
            }
        }
    }

private:
    Index_ my_extent;
    std::vector<unsigned char> my_present;
    Index_ my_offset = 0;
    Index_ my_lastp1 = 0;

public:
    template<class IndexIt_, class Store_>
    void populate(IndexIt_ indices_start, IndexIt_ indices_end, const Store_ store) const {
        if (my_present.empty()) {
            return;
        }

        // Limiting the iteration to its boundaries based on the first and last subset index.
        const auto original_start = indices_start;
        refine_primary_limits(indices_start, indices_end, my_extent, my_offset, my_lastp1);

        for (; indices_start != indices_end; ++indices_start) {
            const auto ix = *indices_start;
            if (my_present[ix - my_offset]) {
                store(indices_start - original_start, ix);
            }
        }
    }
};

}

}

#endif
