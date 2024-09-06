#include <gtest/gtest.h>
#include "../custom_parallel.h"

#include <vector>
#include <algorithm>

#include "tatami/utils/parallelize.hpp"

#ifndef CUSTOM_PARALLEL_TEST
TEST(ParallelizeTest, BasicCheck) {
    std::vector<int> start(3, -1), length(3, -1);
    tatami::parallelize([&](int t, int s, int l) -> void {
        start[t] = s;
        length[t] = l;
    }, 100, 3);

    EXPECT_EQ(start.front(), 0);
    int last = length.front();
    for (int t = 1; t < 3; ++t) {
        EXPECT_EQ(last, start[t]);
        last += length[t];
    }
    EXPECT_EQ(last, 100);

    tatami::parallelize<false>([&](int t, int s, int l) -> void {
        start[t] = s;
        length[t] = l;
    }, 100, 3);
    EXPECT_EQ(start.front(), 0);
    EXPECT_EQ(length.front(), 100);
}

#else
TEST(ParallelizeTest, BasicCheck) {
    std::vector<std::vector<std::pair<int, int> > > parts(3);
    tatami::parallelize([&](int t, int s, int l) -> void {
        parts[t].emplace_back(s, l);
    }, 100, 3);

    std::vector<std::pair<int, int> > sorted_parts;
    for (const auto& p : parts) {
        sorted_parts.insert(sorted_parts.end(), p.begin(), p.end());
    }
    std::sort(sorted_parts.begin(), sorted_parts.end());

    EXPECT_EQ(sorted_parts.front().first, 0);
    int last = sorted_parts.front().second;
    for (size_t i = 1; i < sorted_parts.size(); ++i) {
        const auto& sp = sorted_parts[i];
        EXPECT_EQ(last, sp.first);
        last += sp.second;
    }
    EXPECT_EQ(last, 100);
}
#endif
