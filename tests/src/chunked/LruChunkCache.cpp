#include <gtest/gtest.h>
#include "tatami/tatami.hpp"

TEST(LruChunkCache, Basic) {
    tatami::LruChunkCache<int, std::pair<int, int> > cache(3);

    int counter = 0;
    auto creator = []() -> std::pair<int, int> {
        return std::pair<int, int>(0, 0); 
    };
    auto populator = [&](int i, std::pair<int, int>& chunk) -> void {
        chunk.first = i;
        chunk.second = counter;
        ++counter;
        return;
    };

    auto out = cache.find_chunk(10, creator, populator); // new allocation
    EXPECT_EQ(out.first, 10);
    EXPECT_EQ(out.second, 0);

    out = cache.find_chunk(20, creator, populator); // new allocation
    EXPECT_EQ(out.first, 20);
    EXPECT_EQ(out.second, 1);

    out = cache.find_chunk(10, creator, populator); // retrieve from cache.
    EXPECT_EQ(out.first, 10);
    EXPECT_EQ(out.second, 0);

    out = cache.find_chunk(30, creator, populator); // new allocation, now we're full.
    EXPECT_EQ(out.first, 30);
    EXPECT_EQ(out.second, 2);

    out = cache.find_chunk(20, creator, populator); // retrieve from cache.
    EXPECT_EQ(out.first, 20);
    EXPECT_EQ(out.second, 1);

    out = cache.find_chunk(10, creator, populator); // retrieve from cache.
    EXPECT_EQ(out.first, 10);
    EXPECT_EQ(out.second, 0);

    out = cache.find_chunk(40, creator, populator); // evict the LRU chunk (i.e., 30) and fill with 40.
    EXPECT_EQ(out.first, 40);
    EXPECT_EQ(out.second, 3);

    out = cache.find_chunk(10, creator, populator); // retrieve from cache, as 10 is still present.
    EXPECT_EQ(out.first, 10);
    EXPECT_EQ(out.second, 0);

    out = cache.find_chunk(30, creator, populator); // evict the LRU chunk (i.e., 20) and fill, as it was previously evicted.
    EXPECT_EQ(out.first, 30);
    EXPECT_EQ(out.second, 4);

    out = cache.find_chunk(20, creator, populator); // evict the LRU chunk (i.e., 40) and fill, as it was evicted just previously.
    EXPECT_EQ(out.first, 20);
    EXPECT_EQ(out.second, 5);

    out = cache.find_chunk(10, creator, populator); // retrieve from cache.
    EXPECT_EQ(out.first, 10);
    EXPECT_EQ(out.second, 0);
}
