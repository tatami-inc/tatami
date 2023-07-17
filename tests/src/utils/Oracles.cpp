#include <gtest/gtest.h>
#include "tatami/utils/Oracles.hpp"

#include <random>
#include <vector>

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

    EXPECT_FALSE(streamer.next(prediction)); // any subsequent attempts to call it will always yield false.
}

TEST(OracleStream, FixedAccess) {
    std::mt19937_64 rng(42 * 42);
    std::vector<int> predictions(1234);
    for (auto& x : predictions) {
        x = rng() % 121;
    }

    auto test = std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size());
    tatami::OracleStream<int> streamer(std::move(test));

    auto pIt = predictions.begin();
    int prediction;
    while (streamer.next(prediction)) {
        EXPECT_EQ(*pIt, prediction);
        ++pIt;
    }

    EXPECT_TRUE(pIt == predictions.end());

    EXPECT_FALSE(streamer.next(prediction)); // any subsequent attempts to call it will always yield false.
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
    auto test2 = std::make_unique<tatami::ConsecutiveOracle<int> >(1000, 999);
    streamer.set(std::move(test2));

    int counter = 1000;
    int prediction;
    while (streamer.next(prediction)) {
        EXPECT_EQ(counter, prediction);
        ++counter;
    }

    EXPECT_EQ(counter, 1999);
}
