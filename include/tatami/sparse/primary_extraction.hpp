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
void refine_primary_block_limits(IndexIt_& indices_start, IndexIt_& indices_end, Index_ extent, Index_ block_start, Index_ block_length) {
    if (block_start) {
        // Using custom comparator to ensure that we cast to Index_ for signedness-safe comparisons.
        indices_start = std::lower_bound(indices_start, indices_end, block_start, [](Index_ a, Index_ b) -> bool { return a < b; });
    }

    auto block_end = block_start + block_length;
    if (block_end != extent) {
        indices_end = std::lower_bound(indices_start, indices_end, block_end, [](Index_ a, Index_ b) -> bool { return a < b; });
    }
}

template<class IndexIt_, typename Index_, class Store_>
void retrieve_primary_subset(IndexIt_ indices_start, IndexIt_ indices_end, const std::vector<Index_>& subset, Store_ store) {
    size_t nsub = subset.size();
    if (nsub == 0) {
        return;
    }

    size_t offset = 0;
    if (subset[0]) {
        // Using custom comparator to ensure that we cast to Index_ for signedness-safe comparisons.
        auto new_start = std::lower_bound(indices_start, indices_end, subset[0], [](Index_ a, Index_ b) -> bool { return a < b; });
        offset = new_start - indices_start;
        indices_start = new_start;
    }

    // Looping over the indices in the outer loop and 'subset' in the inner loop.
    // The indices should be sparser than subset, so we get a tighter inner loop.
    size_t s = 0;
    while (indices_start != indices_end) {
        auto curi = *indices_start;
        while (true) {
            if (s == nsub) {
                return;
            }
            if (curi == subset[s]) {
                store(s, offset, curi);
                ++s;
                break;
            }
            ++s;
        }
        ++indices_start;
        ++offset;
    }
}

}

}

#endif
