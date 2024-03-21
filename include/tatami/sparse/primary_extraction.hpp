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

template<typename Index_>
struct RetrievePrimarySubsetDense {
    RetrievePrimarySubsetDense(const std::vector<Index_>& subset) {
        if (!subset.empty()) {
            offset = subset.front();
            size_t alloc = subset.back() - offset + 1;
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
        size_t nmax = present.size();
        if (nmax == 0) {
            return;
        }

        size_t counter = 0;
        for (; indices_start != indices_end; ++indices_start, ++counter) {
            auto ix = *indices_start;
            size_t delta = ix - offset;
            if (delta < nmax) {
                auto shift = present[delta];
                if (shift) {
                    store(shift - 1, counter);
                }
            }
        }
    }

    std::vector<Index_> present;
    size_t offset = 0;
};

struct RetrievePrimarySubsetSparse {
    template<typename Index_>
    RetrievePrimarySubsetSparse(const std::vector<Index_>& subset) {
        if (!subset.empty()) {
            offset = subset.front();
            size_t alloc = subset.back() - offset + 1;
            present.resize(alloc);
            for (auto s : subset) {
                present[s - offset] = 1;
            }
        }
    }

    template<class IndexIt_, class Store_>
    void populate(IndexIt_ indices_start, IndexIt_ indices_end, Store_ store) const {
        size_t nmax = present.size();
        if (nmax == 0) {
            return;
        }

        size_t counter = 0;
        for (; indices_start != indices_end; ++indices_start, ++counter) {
            auto ix = *indices_start;
            size_t delta = ix - offset;
            if (delta < nmax && present[delta]) {
                store(counter, ix);
            }
        }
    }

    std::vector<unsigned char> present;
    size_t offset = 0;
};

}

}

#endif
