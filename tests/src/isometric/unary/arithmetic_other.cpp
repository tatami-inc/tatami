#include <gtest/gtest.h>

#include "tatami/isometric/unary/arithmetic_helpers.hpp"

#include <limits>
#include <vector>

template<tatami::ArithmeticOperation op_, bool right_>
void check_sparsity(double scalar) {
    bool res = tatami::delayed_arithmetic_actual_sparse<op_, right_>(scalar);
    bool ex = tatami::delayed_arithmetic<op_, right_, double, double>(0, scalar) == 0;
    EXPECT_EQ(res, ex);
}

TEST(DelayedUnaryIsometricArithmetic, Sparsity) {
    // Assuming that doubles are IEEE-compatible, the various special values
    // are well-defined and can be produced by arithmetic. We use this to
    // empirically check the sparsity decisions.
    static_assert(std::numeric_limits<double>::is_iec559);

    std::vector<double> scalars {
        0.0,
        1.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    };

    for (auto s : scalars) {
        check_sparsity<tatami::ArithmeticOperation::ADD, true>(s);
        check_sparsity<tatami::ArithmeticOperation::ADD, false>(s);

        check_sparsity<tatami::ArithmeticOperation::SUBTRACT, true>(s);
        check_sparsity<tatami::ArithmeticOperation::SUBTRACT, false>(s);

        check_sparsity<tatami::ArithmeticOperation::MULTIPLY, true>(s);
        check_sparsity<tatami::ArithmeticOperation::MULTIPLY, false>(s);

        check_sparsity<tatami::ArithmeticOperation::DIVIDE, true>(s);
        check_sparsity<tatami::ArithmeticOperation::DIVIDE, false>(s);

        check_sparsity<tatami::ArithmeticOperation::POWER, true>(s);
        check_sparsity<tatami::ArithmeticOperation::POWER, false>(s);

        check_sparsity<tatami::ArithmeticOperation::MODULO, true>(s);
        check_sparsity<tatami::ArithmeticOperation::MODULO, false>(s);

        check_sparsity<tatami::ArithmeticOperation::INTEGER_DIVIDE, true>(s);
        check_sparsity<tatami::ArithmeticOperation::INTEGER_DIVIDE, false>(s);
    }
}

TEST(DelayedUnaryIsometricArithmetic, FillInt) {
    // Checking fill calculations for non-IEEE types.
    {
        auto res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::ADD, true, int>(2);
        EXPECT_EQ(res, 2);

        res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::MULTIPLY, true, int>(2);
        EXPECT_EQ(res, 0);

        res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::POWER, false, int>(2);
        EXPECT_EQ(res, 1);

        res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, true, int>(2);
        EXPECT_EQ(res, 0);

        res = tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::INTEGER_DIVIDE, true, int>(2);
        EXPECT_EQ(res, 0);
    }

    {
        std::string msg;
        try {
            tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, true, int>(0);
        } catch (std::exception& e) {
            msg = e.what();
        }
        EXPECT_TRUE(msg.find("division by zero") != std::string::npos);

        msg.clear();
        try {
            tatami::delayed_arithmetic_zero<tatami::ArithmeticOperation::DIVIDE, false, int>(2);
        } catch (std::exception& e) {
            msg = e.what();
        }
        EXPECT_TRUE(msg.find("division by zero") != std::string::npos);
    }
}
