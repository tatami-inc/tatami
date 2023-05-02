#include <gtest/gtest.h>
#include "tatami/utils/Oracles.hpp"

TEST(OracleStream, BasicAccess) {
    auto test = std::make_unique<tatami::ConsecutiveOracle<int> >(0, 1000);
    tatami::OracleStream<int> streamer(std::move(test));

    int counter = 0;
    int prediction;
    while (streamer.next(prediction)) {
        EXPECT_EQ(counter, prediction);
        ++counter;
    }

    EXPECT_EQ(counter, 1000);
}

TEST(OracleStream, Replacement) {
    tatami::OracleStream<int> streamer;
    EXPECT_FALSE(streamer.active());

    auto test = std::make_unique<tatami::ConsecutiveOracle<int> >(0, 1000);
    streamer.set(std::move(test));
    EXPECT_TRUE(streamer.active());
    
    // A few hits.
    for (int counter = 0; counter < 20; ++counter) {
        int prediction;
        streamer.next(prediction);
        EXPECT_EQ(counter, prediction);
    }

    // Valid after replacement.
    auto test2 = std::make_unique<tatami::ConsecutiveOracle<int> >(1000, 2000);
    streamer.set(std::move(test2));

    int counter = 1000;
    int prediction;
    while (streamer.next(prediction)) {
        EXPECT_EQ(counter, prediction);
        ++counter;
    }

    EXPECT_EQ(counter, 2000);
}
