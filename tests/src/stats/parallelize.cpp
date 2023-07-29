#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/stats/utils.hpp"

TEST(ParallelizeTest, BasicCheck) {
    std::vector<int> start(3, -1), length(3, -1);

    tatami::parallelize([&](size_t t, int s, int l) -> void {
        start[t] = s;
        length[t] = l;
    }, 100, 3);

    EXPECT_EQ(start[0], 0);
    EXPECT_EQ(start[1], 34);
    EXPECT_EQ(start[2], 68);

    EXPECT_EQ(length[0], 34);
    EXPECT_EQ(length[1], 34);
    EXPECT_EQ(length[2], 32);
}

TEST(ParallelizeTest, TypeCheck1) {
    // Checking that the interval calculation is done correctly
    // when we're close to the overflow boundary.
    std::vector<unsigned char> start(2, -1), length(2, -1);
    tatami::parallelize([&](size_t t, unsigned char s, unsigned char l) -> void {
        start[t] = s;
        length[t] = l;
    }, static_cast<unsigned char>(255), 2);

    EXPECT_EQ(start[0], 0);
    EXPECT_EQ(start[1], 128);

    EXPECT_EQ(length[0], 128);
    EXPECT_EQ(length[1], 127);
}

TEST(ParallelizeTest, TypeCheck2) {
    // No overflow in the number of jobs.
    std::vector<unsigned char> start(1000, -1), length(1000, -1);
    tatami::parallelize([&](size_t t, unsigned char s, unsigned char l) -> void {
        start[t] = s;
        length[t] = l;
    }, static_cast<unsigned char>(2), 1000);

    EXPECT_EQ(start[0], 0);
    EXPECT_EQ(start[1], 1);

    EXPECT_EQ(length[0], 1);
    EXPECT_EQ(length[1], 1);

    start[0] = -1;
    length[0] = -1;
    start[1] = -1;
    length[1] = -1;
    EXPECT_EQ(start, std::vector<unsigned char>(1000, -1));
    EXPECT_EQ(length, std::vector<unsigned char>(1000, -1));
}

#ifndef CUSTOM_PARALLEL_TEST
TEST(ParallelizeTest, ErrorChecks) {
    EXPECT_ANY_THROW({
        try {
            tatami::parallelize([&](size_t, int, int) -> void {
                throw std::runtime_error("WHEE");
            }, 255, 2);
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("WHEE") != std::string::npos);
            throw;
        }
    });

    EXPECT_ANY_THROW({
        try {
            tatami::parallelize([&](size_t, int, int) -> void {
                throw 123;
            }, 255, 2);
        } catch (std::exception& e) {
            EXPECT_TRUE(std::string(e.what()).find("unknown error") != std::string::npos);
            throw;
        }
    });
}
#endif
