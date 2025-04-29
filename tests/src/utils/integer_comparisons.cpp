#include <gtest/gtest.h>
#include "tatami/utils/integer_comparisons.hpp"

TEST(IntegerComparisons, SafeNonNegativeEqual) {
    EXPECT_TRUE(tatami::safe_non_negative_equal(0, 0));
    EXPECT_FALSE(tatami::safe_non_negative_equal(-1, -1));
    EXPECT_FALSE(tatami::safe_non_negative_equal(0, -1));
    EXPECT_FALSE(tatami::safe_non_negative_equal(0, 1));
    EXPECT_TRUE(tatami::safe_non_negative_equal(1, 1));

    // Trying with a mix of signed and unsigned types.
    EXPECT_TRUE(tatami::safe_non_negative_equal(1u, 1));
    EXPECT_FALSE(tatami::safe_non_negative_equal(1u, -1));
    EXPECT_TRUE(tatami::safe_non_negative_equal(1, 1u));
    EXPECT_FALSE(tatami::safe_non_negative_equal(-1, 1u));
    EXPECT_TRUE(tatami::safe_non_negative_equal(1u, 1u));
}
