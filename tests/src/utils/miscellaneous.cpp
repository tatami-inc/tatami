#include <gtest/gtest.h>
#include "tatami/utils/has_data.hpp"

TEST(HasData, Basic) {
    {
        constexpr bool found = tatami::has_data<int, std::vector<int> >::value;
        EXPECT_TRUE(found);
    }

    {
        constexpr bool found = tatami::has_data<int, const std::vector<int> >::value;
        EXPECT_TRUE(found);
    }

    {
        constexpr bool found = tatami::has_data<double, std::vector<int> >::value;
        EXPECT_FALSE(found);
    }

    {
        constexpr bool found = tatami::has_data<double, std::pair<double, double> >::value;
        EXPECT_FALSE(found);
    }
}
