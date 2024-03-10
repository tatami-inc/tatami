#include <gtest/gtest.h>
#include "tatami/utils/ConsecutiveOracle.hpp"

#include <random>
#include <vector>

class TestConsecutiveOracle : public ::testing::TestWithParam<std::tuple<int, int> > {};

TEST_P(TestConsecutiveOracle, BasicAccess) {
    auto param = GetParam();
    auto start = std::get<0>(param);
    auto len = std::get<1>(param);

    auto test = std::make_unique<tatami::ConsecutiveOracle<int> >(start, len);
    EXPECT_EQ(test->total(), len);
    for (int i = 0; i < len; ++i) {
        EXPECT_EQ(test->get(i), i + start);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ConsecutiveOracle,
    TestConsecutiveOracle,
    ::testing::Combine( 
        ::testing::Values(0, 10, 100),
        ::testing::Values(0, 100, 1000)
    )
);
