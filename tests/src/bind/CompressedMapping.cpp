#include <gtest/gtest.h>

#include <vector>
#include <cstdint>
#include <random>

#include "tatami/bind/CompressedMapping.hpp"

class DelayedBindCompressedMappingTest : public ::testing::TestWithParam<std::tuple<int, int, int> > {};

TEST_P(DelayedBindCompressedMappingTest, Basic) {
    auto params = GetParam();
    auto min_size = std::get<0>(params);
    auto total_dim = std::get<1>(params); 
    auto variation = std::get<2>(params);

    // i.e., intervals range from [min_size, (variation + 1)*min_size].
    const int extra = variation * min_size + 1;

    std::mt19937_64 rng(min_size + total_dim * 13 + variation * 103);
    std::vector<int> cumulative(1);
    while (cumulative.back() != total_dim) {
        const int last = cumulative.back();
        const int chosen = min_size + (rng() % extra);
        auto actual = std::min(chosen, total_dim - last); 
        cumulative.push_back(last + actual);
    }

    // Merging the last interval with the second-last interval if it's too small.
    ASSERT_GE(cumulative.size(), 2);
    if (cumulative.back() - cumulative[cumulative.size() - 2] < min_size) {
        cumulative.pop_back();
        cumulative.back() = total_dim;
    }
    ASSERT_GE(cumulative.size(), 2);

    tatami::DelayedBind_internal::CompressedMapping<int> mapping(min_size, cumulative);
    for (int i = 0, end = cumulative.back(); i < end; ++i) {
        auto startp1 = cumulative.begin() + 1;
        EXPECT_EQ(mapping[i], std::upper_bound(startp1, cumulative.end(), i) - startp1);
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindCompressedMappingTest,
    ::testing::Combine(
        ::testing::Values(1, 2, 4, 5, 8, 11, 16, 20, 32, 51), // min size, alternating between powers of 2 and not.
        ::testing::Values(100, 256, 333), // total dimensions. We try even, powers of two, and odd.
        ::testing::Values(0, 2) // variation in intervals
    )
);

TEST(DelayedBind, CompressedMappingSmallType) {
    // Here we try to push the limits of what we can compact into the storage type.
    // Specifically, we try a bunch of bin sizes up to the maximum possible value for the cumulative size .

    // - Bin size is 64.
    // - Bin encodings are [0 << 6, 64 - 1], [0 << 6, 63 - 1]
    {
        std::int8_t min_size = 127;
        std::vector<std::int8_t> cumulative{ 0, 127 };
        tatami::DelayedBind_internal::CompressedMapping<std::int8_t> mapping(min_size, cumulative);
        for (int i = 0, end = cumulative.back(); i < end; ++i) {
            EXPECT_EQ(mapping[i], 0);
        }
    }

    // - Bin size is 32.
    // - Bin encodings are [0 << 5, 32 - 1], [0 << 5, 28 - 1], [1 << 5, 32 - 1], [1 << 5, 31 - 1]
    {
        std::int8_t min_size = 60;
        std::vector<std::int8_t> cumulative{ 0, 60, 127 };
        tatami::DelayedBind_internal::CompressedMapping<std::int8_t> mapping(min_size, cumulative);
        for (int i = 0, end = cumulative.back(); i < end; ++i) {
            auto startp1 = cumulative.begin() + 1;
            EXPECT_EQ(mapping[i], std::upper_bound(startp1, cumulative.end(), i) - startp1);
        }
    }

    // - Bin size is 16.
    // - Bin encodings are [0 << 4, 16 - 1], [0 << 4, 5 - 1],
    //                     [1 << 4, 12 - 1],
    //                     [2 << 4, 16 - 1], [2 << 4, 16 - 1], [2 << 4, 2 - 1],
    //                     [3 << 4, 7 - 1],
    //                     [4 << 4, 15 - 1]
    {
        std::int8_t min_size = 20;
        std::vector<std::int8_t> cumulative{ 0, 21, 44, 82, 103, 127 };
        tatami::DelayedBind_internal::CompressedMapping<std::int8_t> mapping(min_size, cumulative);
        for (int i = 0, end = cumulative.back(); i < end; ++i) {
            auto startp1 = cumulative.begin() + 1;
            EXPECT_EQ(mapping[i], std::upper_bound(startp1, cumulative.end(), i) - startp1);
        }
    }

    // - Bin size is 2.
    // - Bin encodings are [i << 1, 2 - 1] for all i in [0, 63), and then [63 << 1, 0].
    {
        std::int8_t min_size = 2;
        std::vector<std::int8_t> cumulative;
        for (int i = 0; i < 64; ++i) {
            cumulative.push_back(i * 2);
        }
        tatami::DelayedBind_internal::CompressedMapping<std::int8_t> mapping(min_size, cumulative);
        for (int i = 0, end = cumulative.back(); i < end; ++i) {
            EXPECT_EQ(mapping[i], i / 2);
        }
    }

    // - Bin size is 1.
    // - Bin encodings are [i << 0, 0] for all i in [0, 127).
    {
        std::int8_t min_size = 1;
        std::vector<std::int8_t> cumulative(128);
        std::iota(cumulative.begin(), cumulative.end(), 0);
        tatami::DelayedBind_internal::CompressedMapping<std::int8_t> mapping(min_size, cumulative);
        for (int i = 0, end = cumulative.back(); i < end; ++i) {
            EXPECT_EQ(mapping[i], i);
        }
    }
}

TEST(DelayedBind, CompressedMappingBigType) {
    // Checking that the binning approach still works when the numbers get pretty large.
    // Now the bin size is 2^19.
    {
        std::uint32_t min_size = 1000000;
        std::vector<std::uint32_t> cumulative{ 0, 1000000, 4000000, 6000000, 10000000 };
        tatami::DelayedBind_internal::CompressedMapping<std::uint32_t> mapping(min_size, cumulative);
        EXPECT_EQ(mapping[0], 0);
        EXPECT_EQ(mapping[999999], 0);
        EXPECT_EQ(mapping[1000000], 1);
        EXPECT_EQ(mapping[1000001], 1);
        EXPECT_EQ(mapping[2000000], 1);
        EXPECT_EQ(mapping[3999999], 1);
        EXPECT_EQ(mapping[4000000], 2);
        EXPECT_EQ(mapping[4000001], 2);
        EXPECT_EQ(mapping[5000000], 2);
        EXPECT_EQ(mapping[5999999], 2);
        EXPECT_EQ(mapping[6000000], 3);
        EXPECT_EQ(mapping[6000001], 3);
        EXPECT_EQ(mapping[8000000], 3);
        EXPECT_EQ(mapping[9999999], 3);
    }

    // Even bigger, bin size is now 2^28.
    {
        std::int32_t min_size = 500000000;
        std::vector<std::int32_t> cumulative{ 0, 1000000000, 1600000000, 2147483647 };
        tatami::DelayedBind_internal::CompressedMapping<std::int32_t> mapping(min_size, cumulative);
        EXPECT_EQ(mapping[0], 0);
        EXPECT_EQ(mapping[999999999], 0);
        EXPECT_EQ(mapping[1000000000], 1);
        EXPECT_EQ(mapping[1000000001], 1);
        EXPECT_EQ(mapping[1599999999], 1);
        EXPECT_EQ(mapping[1600000000], 2);
        EXPECT_EQ(mapping[2000000001], 2);
        EXPECT_EQ(mapping[2147483646], 2);
    }
}
