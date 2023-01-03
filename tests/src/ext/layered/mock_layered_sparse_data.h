#ifndef MOCK_LAYERED_SPARSE_DATA_H
#define MOCK_LAYERED_SPARSE_DATA_H

#include <random>
#include <vector>

template<bool ByRow>
void mock_layered_sparse_data(size_t NR, size_t NC, std::vector<size_t>& rows, std::vector<size_t>& cols, std::vector<int>& vals) {
    std::mt19937_64 rng(ByRow ? NR : NC);

    // This idea is to always have some values in all categories;
    // we wouldn't get this if we uniformly sampled.
    auto rand = [&]() -> int {
        int upper = 1;
        switch (rng() % 3) {
            case 0:
                upper = 256;
                break;
            case 1:
                upper = 65536;
                break;
            case 2:
                upper = 10000000; // that's big enough.
                break;
        }
        return rng() % upper;
    };

    if constexpr(!ByRow) {
        for (size_t c = 0; c < NC; ++c) {
            for (size_t r = 0; r < NR; ++r) {
                // i.e., around about two values per row. This gives
                // us something interesting w.r.t. multiple values at different limits
                // (e.g., about 1/9th of the rows will have both values at the lowest limit).
                if (rng() % NC < 2) { 
                    rows.push_back(r);
                    cols.push_back(c);
                    vals.push_back(rand());
                }
            }
        }
    } else {
        for (size_t r = 0; r < NR; ++r) {
            for (size_t c = 0; c < NC; ++c) {
                if (rng() % NC < 2) { 
                    rows.push_back(r);
                    cols.push_back(c);
                    vals.push_back(rand());
                }
            }
        }
    }

    return;
}

#endif
