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
    auto res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::ADD, true, int>(2);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::MULTIPLY, true, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::POWER, false, int>(2);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::DIVIDE, true, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int>(2);
    EXPECT_TRUE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int>(0);
    EXPECT_FALSE(res);

    res = tatami::delayed_arithmetic_actual_sparse<tatami::ArithmeticOperation::INTEGER_DIVIDE, false, int>(2);
    EXPECT_FALSE(res);
}

TEST(DelayedUnaryIsometricArithmetic, NonIeee754Fill) {
    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::ADD, true, int>(2);
        EXPECT_EQ(res, 2);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::MULTIPLY, true, int>(2);
        EXPECT_EQ(res, 0);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::POWER, false, int>(2);
        EXPECT_EQ(res, 1);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, true, int>(2);
        EXPECT_EQ(res, 0);
    }

    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int>(2);
        EXPECT_EQ(res, 0);
    }

    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, true, int>(0); }, "division by zero");
    tatami_test::throws_error([&]() { tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, false, int>(2); }, "division by zero");
}

TEST(DelayedUnaryIsometricArithmetic, NonIeee754Ops) {
    {
        int scalar = 5;
        auto op = tatami::make_DelayedUnaryIsometricMultiplyScalar(scalar);
        EXPECT_TRUE(op.is_sparse());
    }

    {
        int scalar = 0;
        auto op = tatami::make_DelayedUnaryIsometricPowerScalar<false>(scalar);
        EXPECT_FALSE(op.is_sparse());
    }

    {
        auto op = tatami::make_DelayedUnaryIsometricDivideScalar<false, int>(5);
        tatami_test::throws_error([&]() { op.template fill<int, int>(true, 5); }, "division by zero");
    }
}

TEST(DelayedUnaryIsometricArithmetic, NonFiniteMultiply) {
    double scalar = std::numeric_limits<double>::infinity();
    auto op = tatami::make_DelayedUnaryIsometricMultiplyScalar(scalar);
    EXPECT_FALSE(op.is_sparse());
}
