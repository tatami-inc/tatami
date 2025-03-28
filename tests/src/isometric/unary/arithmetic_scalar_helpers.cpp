#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/arithmetic_helpers.hpp"
#include "tatami/isometric/unary/math_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricUtilsScalar { 
protected:
    inline static size_t nrow = 123, ncol = 89;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void assemble() {
        if (dense) {
            return;
        }

        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.seed = 511498129;
            return opt;
        }());

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.
    }
};

/*******************************
 ********* COMMUTATIVE *********
 *******************************/

class DelayedUnaryIsometricArithmeticCommutativeScalarTest : public ::testing::TestWithParam<double>, public DelayedUnaryIsometricUtilsScalar {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricArithmeticCommutativeScalarTest, Add) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAddScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_FALSE(dense_mod.is_sparse());

    if (val) {
        EXPECT_FALSE(sparse_mod.is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected.
    auto refvec = simulated;
    for (auto& r : refvec) {
        r += val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticCommutativeScalarTest, Multiply) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricMultiplyScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());
    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r *= val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticCommutativeScalarTest, NewType) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAddScalarHelper<float, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, op);
    tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, op);

    size_t nsize = simulated.size();
    auto frefvec = std::vector<float>(nsize);
    for (size_t i = 0; i < nsize; ++i) {
        frefvec[i] = simulated[i] + val;
    }

    tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);
    quick_test_all<float, int>(dense_fmod, fref);
    quick_test_all<float, int>(sparse_fmod, fref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticScalar,
    DelayedUnaryIsometricArithmeticCommutativeScalarTest,
    ::testing::Values(5, 0.1, -0.7, 0)
);

/***********************************
 ********* NON-COMMUTATIVE *********
 ***********************************/

class DelayedUnaryIsometricArithmeticNonCommutativeScalarTest : public ::testing::TestWithParam<std::tuple<double, bool> >, public DelayedUnaryIsometricUtilsScalar {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Subtract) {
    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (on_right) {
        op.reset(new tatami::DelayedUnaryIsometricSubtractScalarHelper<true, double, double, int, double>(val));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricSubtractScalarHelper<false, double, double, int, double>(val));
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod.is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod.is_sparse());
    }

    // Again, doing some light tests.
    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r -= val;
        } else {
            r = val - r;
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Divide) {
    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (on_right) {
        op.reset(new tatami::DelayedUnaryIsometricDivideScalarHelper<true, double, double, int, double>(val));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricDivideScalarHelper<false, double, double, int, double>(val));
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    if (on_right && val) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = careful_division(r, val);
        } else {
            r = careful_division(val, r);
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Power) {
    auto asimulated = simulated;
    for (auto& a : asimulated) {
        a = std::abs(a);
    }
    auto adense = std::make_shared<tatami::DenseMatrix<double, int, decltype(asimulated)> >(nrow, ncol, asimulated, true); // row major.
    auto asparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.

    auto my_param = GetParam();
    double val = std::abs(std::get<0>(my_param));
    bool on_right = std::get<1>(my_param);
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (on_right) {
        op.reset(new tatami::DelayedUnaryIsometricPowerScalarHelper<true, double, double, int, double>(val));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricPowerScalarHelper<false, double, double, int, double>(val));
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Again, doing some light tests.
    auto refvec = asimulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = std::pow(r, val);
        } else {
            r = std::pow(val, r);
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Modulo) {
    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (on_right) {
        op.reset(new tatami::DelayedUnaryIsometricModuloScalarHelper<true, double, double, int, double>(val));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricModuloScalarHelper<false, double, double, int, double>(val));
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = careful_modulo(r, val);
        } else {
            r = careful_modulo(val, r);
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, IntegerDivide) {
    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (on_right) {
        op.reset(new tatami::DelayedUnaryIsometricIntegerDivideScalarHelper<true, double, double, int, double>(val));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricIntegerDivideScalarHelper<false, double, double, int, double>(val));
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    auto refvec = simulated;
    for (auto& r : refvec) {
        // x == (x %% y) + y * (x %/% y)
        if (on_right) {
            r = std::floor(careful_division(r, val));
        } else {
            r = std::floor(careful_division(val, r));
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, NewType) {
    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<float, double, int> > op;
    if (on_right) {
        op.reset(new tatami::DelayedUnaryIsometricSubtractScalarHelper<true, float, double, int, double>(val));
    } else {
        op.reset(new tatami::DelayedUnaryIsometricSubtractScalarHelper<false, float, double, int, double>(val));
    }
    tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, op);
    tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, op);

    size_t nsize = simulated.size();
    auto frefvec = std::vector<float>(nsize);
    for (size_t i = 0; i < nsize; ++i) {
        if (on_right) {
            frefvec[i] = simulated[i] - val;
        } else {
            frefvec[i] = val - simulated[i];
        }
    }

    tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);
    quick_test_all<float, int>(dense_fmod, fref);
    quick_test_all<float, int>(sparse_fmod, fref);
}


INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricArithmeticScalar,
    DelayedUnaryIsometricArithmeticNonCommutativeScalarTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7, 0),
        ::testing::Values(true, false)
    )
);

/**********************************
 ********* SPECIAL VALUES *********
 **********************************/

TEST(DelayedUnaryIsometricArithmetic, NonIeee754Ops) {
    {
        tatami::DelayedUnaryIsometricMultiplyScalarHelper<int, int, int, int> op(5);
        EXPECT_TRUE(op.is_sparse());
        EXPECT_EQ(op.fill(true, 5), 0);
    }

    {
        tatami::DelayedUnaryIsometricPowerScalarHelper<true, int, int, int, int> op(5);
        EXPECT_TRUE(op.is_sparse());
        EXPECT_EQ(op.fill(true, 5), 0);
    }

    {
        tatami::DelayedUnaryIsometricPowerScalarHelper<true, int, int, int, int> op(0);
        EXPECT_FALSE(op.is_sparse());
        EXPECT_EQ(op.fill(true, 5), 1);
    }

    {
        tatami::DelayedUnaryIsometricDivideScalarHelper<false, int, int, int, int> op(5);
        EXPECT_FALSE(op.is_sparse());
        tatami_test::throws_error([&]() { op.fill(true, 5); }, "division by zero");
    }

    {
        tatami::DelayedUnaryIsometricIntegerDivideScalarHelper<false, int, int, int, int> op(5);
        EXPECT_FALSE(op.is_sparse());
        tatami_test::throws_error([&]() { op.fill(true, 5); }, "division by zero");
    }

    {
        tatami::DelayedUnaryIsometricModuloScalarHelper<false, int, int, int, int> op(5);
        EXPECT_FALSE(op.is_sparse());
        tatami_test::throws_error([&]() { op.fill(true, 5); }, "division by zero");
    }
}

TEST(DelayedUnaryIsometricArithmetic, NonFiniteMultiply) {
    double scalar = std::numeric_limits<double>::infinity();
    tatami::DelayedUnaryIsometricMultiplyScalarHelper<double, double, int, double> op(scalar);
    EXPECT_FALSE(op.is_sparse());
}

TEST(DelayedUnaryIsometricArithmeticScalar, BackCompatibility) {
    auto add = tatami::make_DelayedUnaryIsometricAddScalar(1);
    EXPECT_FALSE(add->is_sparse());
    auto sub = tatami::make_DelayedUnaryIsometricSubtractScalar<true>(1);
    EXPECT_FALSE(sub->is_sparse());
    auto mult = tatami::make_DelayedUnaryIsometricMultiplyScalar(1);
    EXPECT_TRUE(mult->is_sparse());
    auto div = tatami::make_DelayedUnaryIsometricDivideScalar<true>(1);
    EXPECT_TRUE(div->is_sparse());
    auto mod = tatami::make_DelayedUnaryIsometricModuloScalar<true>(1);
    EXPECT_TRUE(mod->is_sparse());
    auto pow = tatami::make_DelayedUnaryIsometricPowerScalar<true>(1);
    EXPECT_TRUE(pow->is_sparse());
    auto idiv = tatami::make_DelayedUnaryIsometricIntegerDivideScalar<true>(1);
    EXPECT_TRUE(idiv->is_sparse());
}
