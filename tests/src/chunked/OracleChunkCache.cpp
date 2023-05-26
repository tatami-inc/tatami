#include <gtest/gtest.h>
#include "tatami/tatami.hpp"

class OracleChunkCacheTest : public ::testing::Test {
protected:
    struct TestChunk {
        unsigned char chunk_id;
        int cycle_number;
        bool alloc = false;
    };

    template<class Cache_>
    auto next(Cache_& cache, int& counter, int& nalloc) {
        return cache.next_chunk(
            [](int i) -> std::pair<unsigned char, int> {
                return std::make_pair<unsigned char, int>(i / 10, i % 10);
            },
            [](TestChunk& x, TestChunk& y) -> void {
                std::swap(x, y);
            },
            [](const TestChunk& chunk) -> bool {
                return chunk.alloc;
            },
            [&](TestChunk& chunk) -> void {
                ++nalloc;
                chunk.alloc = true;
            },
            [&](const std::vector<std::pair<unsigned char, int> >& chunks_in_need, std::vector<TestChunk>& chunk_data) -> void {
                for (const auto& x : chunks_in_need) {
                    auto& current = chunk_data[x.second];
                    current.chunk_id = x.first;
                    current.cycle_number = counter++;
                }
            }
        );
    }
};

TEST_F(OracleChunkCacheTest, Consecutive) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        33,
        44, // Cycle 2
        55,
        66,
        77, // Cycle 3
        88,
        99
    };

    tatami::OracleChunkCache<unsigned char, int, TestChunk> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 100, 3);
    int counter = 0;
    int nalloc = 0;

    for (int i = 1; i < 9; ++i) {
        auto out = next(cache, counter, nalloc); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(i));
        EXPECT_EQ(out.first->cycle_number, i - 1);
        EXPECT_TRUE(out.first->alloc);
        EXPECT_EQ(out.second, i);
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

TEST_F(OracleChunkCacheTest, AllPredictions) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        12, 
        31,
        23, // Cycle 2
        14, 
        28,
        45, 
        11, // Cycle 3
        36,
        32,
        42,
        24, // Cycle 4
        15
    };

    tatami::OracleChunkCache<unsigned char, int, TestChunk> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 100, 3);
    int counter = 0;
    int nalloc = 0;

    auto out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 22.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->cycle_number, 1);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 31.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->cycle_number, 2);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->cycle_number, 1);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 14.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 28.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->cycle_number, 1);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 8);

    out = next(cache, counter, nalloc); // fetching 45.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->cycle_number, 3);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 5);

    out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 36.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->cycle_number, 4);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 6);

    out = next(cache, counter, nalloc); // fetching 32.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->cycle_number, 4);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 42.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->cycle_number, 3);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 24.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->cycle_number, 5);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 15.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 5);

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

TEST_F(OracleChunkCacheTest, LimitedPredictions) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        12, 
        31, // Cycle 2 (forced by limited predictions)
        23, 
        14, 
        45, // Cycle 3
        51, 
        32, 
        34, // Cycle 4 (forced by limited predictions)
        15 
    };

    tatami::OracleChunkCache<unsigned char, int, TestChunk> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 3, 3);
    int counter = 0;
    int nalloc = 0;

    auto out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 22.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->cycle_number, 1);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 31.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->cycle_number, 2);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->cycle_number, 1);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 14.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 0);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 45.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->cycle_number, 3);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 5);

    out = next(cache, counter, nalloc); // fetching 51.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(5));
    EXPECT_EQ(out.first->cycle_number, 4);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 32.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->cycle_number, 2);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 34.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->cycle_number, 2);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 15.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->cycle_number, 5);
    EXPECT_TRUE(out.first->alloc);
    EXPECT_EQ(out.second, 5);

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

