#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
#include <algorithm>

template<typename Function_, typename Index_>
void foo_parallel(Function_ fun, Index_ ntasks, int nthreads) {
    Index_ tasks_per_thread = ntasks / nthreads + (ntasks % nthreads > 0);
    Index_ start = 0;
    for (int t = 0; t < nthreads; ++t) {
        Index_ length = std::min(tasks_per_thread, ntasks - start);
        fun(t, start, length);
        start += length;
    }
}

// Put this before any tatami imports.
#define TATAMI_CUSTOM_PARALLEL foo_parallel
#endif

#include "tatami/utils/parallelize.hpp"

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
