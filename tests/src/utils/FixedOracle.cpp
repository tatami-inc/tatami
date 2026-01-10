#include <gtest/gtest.h>
#include "tatami/utils/FixedOracle.hpp"

#include <random>
#include <vector>

TEST(FixedOracle, BasicAccess) {
    std::mt19937_64 rng(42 * 42);
    std::vector<int> predictions(1234);
    for (auto& x : predictions) {
        x = rng() % 121;
    }

    auto test = std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size());
    auto test_copy = std::make_unique<tatami::FixedVectorOracle<int> >(predictions);
    EXPECT_EQ(test->total(), predictions.size());
    EXPECT_EQ(test_copy->total(), predictions.size());

    for (std::size_t i = 0; i < predictions.size(); ++i) {
        EXPECT_EQ(test->get(i), predictions[i]);
        EXPECT_EQ(test_copy->get(i), predictions[i]);
    }
}

TEST(FixedOracle, OtherType) {
    // Check that it works with a different type for the underlying sequence.
    std::mt19937_64 rng(42 * 42);
    std::vector<unsigned> predictions(1234);
    for (auto& x : predictions) {
        x = rng() % 13;
    }

    auto test = std::make_unique<tatami::FixedViewOracle<int, const unsigned*> >(predictions.data(), predictions.size());
    auto test_copy = std::make_unique<tatami::FixedVectorOracle<int, std::vector<unsigned> > >(predictions);
    EXPECT_EQ(test->total(), predictions.size());
    EXPECT_EQ(test_copy->total(), predictions.size());

    for (std::size_t i = 0; i < predictions.size(); ++i) {
        EXPECT_EQ(test->get(i), predictions[i]);
        EXPECT_EQ(test_copy->get(i), predictions[i]);
    }
}
