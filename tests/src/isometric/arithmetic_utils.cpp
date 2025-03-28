#include <gtest/gtest.h>

#include "tatami/isometric/unary/arithmetic_helpers.hpp"
#include "tatami_test/tatami_test.hpp"

#include <limits>
#include <vector>

template<tatami::ArithmeticOperation op_, bool right_>
void check_sparsity(double scalar) {
    bool res = tatami::delayed_arithmetic_actual_sparse<op_, right_>(scalar);
    bool ex = tatami::delayed_arithmetic<op_, right_, double, double>(0, scalar) == 0;
    EXPECT_EQ(res, ex);
}

TEST(DelayedUnaryIsometricArithmetic, NonIeee754Sparsity) {
    // Checking sparsity calculations for non-IEEE types.
    auto res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::ADD, true, int, int>(2);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::MULTIPLY, true, int, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::POWER, false, int, int>(2);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::DIVIDE, true, int, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::MODULO, true, int, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::DIVIDE, true, int, int>(0);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::MODULO, true, int, int>(0);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::DIVIDE, false, int, int>(2);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::MODULO, false, int, int>(2);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::INTEGER_DIVIDE, false, int, int>(2);
    EXPECT_FALSE(res);
}

TEST(DelayedUnaryIsometricArithmetic, NonIeee754Fill) {
    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::ADD, true, int, int>(2);
        EXPECT_EQ(res, 2);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::MULTIPLY, true, int, int>(2);
        EXPECT_EQ(res, 0);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::POWER, false, int, int>(2);
        EXPECT_EQ(res, 1);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, true, int, int>(2);
        EXPECT_EQ(res, 0);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::MODULO, true, int, int>(2);
        EXPECT_EQ(res, 0);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int, int>(2);
        EXPECT_EQ(res, 0);
    }

    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, true, int, int>(0); }, "division by zero");
    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, false, int, int>(2); }, "division by zero");

    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::MODULO, true, int, int>(0); }, "division by zero");
    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::MODULO, false, int, int>(2); }, "division by zero");

    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int, int>(0); }, "division by zero");
    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::INTEGER_DIVIDE, false, int, int>(2); }, "division by zero");
}

TEST(DelayedUnaryIsometricArithmetic, TrueIntegerDivide) {
    auto res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10, 5);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10, -5);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10, 5);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10, -5);
    EXPECT_EQ(res, 2);

    // Now with some remainders on the left.
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10, 3);
    EXPECT_EQ(res, 3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10, 3);
    EXPECT_EQ(res, -4);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10, -3);
    EXPECT_EQ(res, -4);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10, -3);
    EXPECT_EQ(res, 3);

    // And on the right.
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(5, 14);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(-5, 14);
    EXPECT_EQ(res, -3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(5, -14);
    EXPECT_EQ(res, -3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(-5, -14);
    EXPECT_EQ(res, 2);
}

TEST(DelayedUnaryIsometricArithmetic, DoublishIntegerDivide) {
    auto res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10.0, 5);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10.0, -5);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10.0, 5);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10.0, -5);
    EXPECT_EQ(res, 2);

    // Now with some fractions on the left.
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10.0, 3);
    EXPECT_EQ(res, 3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10.0, 3);
    EXPECT_EQ(res, -4);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(10.0, -3);
    EXPECT_EQ(res, -4);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(-10.0, -3);
    EXPECT_EQ(res, 3);

    // And on the right, for some testing.
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(5, 14.0);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(-5, 14.0);
    EXPECT_EQ(res, -3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(5, -14.0);
    EXPECT_EQ(res, -3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(-5, -14.0);
    EXPECT_EQ(res, 2);
}

TEST(DelayedUnaryIsometricArithmetic, TrueIntegerModulo) {
    auto res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, 5);
    EXPECT_EQ(res, 0);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, -5);
    EXPECT_EQ(res, 0);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, 5);
    EXPECT_EQ(res, 0);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, -5);
    EXPECT_EQ(res, 0);

    // Non-zero modulo on the left:
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, 3);
    EXPECT_EQ(res, 1);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, 3);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, -3);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, -3);
    EXPECT_EQ(res, -1);

    // And on the right:
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(5, 13);
    EXPECT_EQ(res, 3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(-5, 13);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(5, -13);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(-5, -13);
    EXPECT_EQ(res, -3);
}

TEST(DelayedUnaryIsometricArithmetic, DoublishIntegerModulo) {
    auto res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, 5.0);
    EXPECT_EQ(res, 0);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, -5.0);
    EXPECT_EQ(res, 0);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, 5.0);
    EXPECT_EQ(res, 0);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, -5.0);
    EXPECT_EQ(res, 0);

    // Non-zero modulo on the left:
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, 3.0);
    EXPECT_EQ(res, 1);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, 3.0);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(10, -3.0);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, true>(-10, -3.0);
    EXPECT_EQ(res, -1);

    // And on the right:
    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(5.0, 13.0);
    EXPECT_EQ(res, 3);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(-5.0, 13.0);
    EXPECT_EQ(res, -2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(5.0, -13.0);
    EXPECT_EQ(res, 2);

    res = tatami::delayed_arithmetic<tatami::ArithmeticOperation::MODULO, false>(-5.0, -13.0);
    EXPECT_EQ(res, -3);
}
