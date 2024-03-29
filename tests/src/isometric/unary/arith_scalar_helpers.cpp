#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class ArithScalarUtils { 
protected:
    inline static size_t nrow = 123, ncol = 89;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void assemble() {
        if (dense) {
            return;
        }
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); // column major.
    }
};

/*******************************
 ********* COMMUTATIVE *********
 *******************************/

class ArithCommutativeScalarTest : public ::testing::TestWithParam<double>, public ArithScalarUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithCommutativeScalarTest, Addition) {
    double val = GetParam();
    auto op = tatami::make_DelayedAddScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod->sparse());
    } else {
        EXPECT_TRUE(sparse_mod->sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected.
    auto refvec = simulated;
    for (auto& r : refvec) {
        r += val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(ArithCommutativeScalarTest, Multiplication) {
    double val = GetParam();
    auto op = tatami::make_DelayedMultiplyScalarHelper(val);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r *= val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

INSTANTIATE_TEST_SUITE_P(
    ArithScalar,
    ArithCommutativeScalarTest,
    ::testing::Values(5, 0.1, -0.7, 0)
);

/***********************************
 ********* NON-COMMUTATIVE *********
 ***********************************/

class ArithNonCommutativeScalarTest : public ::testing::TestWithParam<std::tuple<double, bool> >, public ArithScalarUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(ArithNonCommutativeScalarTest, Subtraction) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedSubtractScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedSubtractScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod->sparse());
    } else {
        EXPECT_TRUE(sparse_mod->sparse());
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
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(ArithNonCommutativeScalarTest, Division) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedDivideScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedDivideScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = careful_division(r, val);
        } else {
            r = careful_division(val, r);
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ (val == 0));
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ (val == 0));
}

TEST_P(ArithNonCommutativeScalarTest, Power) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    tatami::DelayedAbsHelper op0;
    dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op0);
    sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op0);

    auto my_param = GetParam();
    double val = std::abs(std::get<0>(my_param));
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedPowerScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense_mod, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse_mod, op);
    } else {
        auto op = tatami::make_DelayedPowerScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense_mod, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse_mod, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    // Again, doing some light tests.
    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = std::pow(std::abs(r), val);
        } else {
            r = std::pow(val, std::abs(r));
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(ArithNonCommutativeScalarTest, Modulo) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedModuloScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedModuloScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = std::fmod(r, val);
        } else {
            r = std::fmod(val, r);
        }
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ !(on_right && val));
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ !(on_right && val));
}

TEST_P(ArithNonCommutativeScalarTest, IntegerDivision) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedIntegerDivideScalarHelper<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    } else {
        auto op = tatami::make_DelayedIntegerDivideScalarHelper<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);
    }

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
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
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ !(on_right && val));
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ !(on_right && val));
}

INSTANTIATE_TEST_SUITE_P(
    ArithScalar,
    ArithNonCommutativeScalarTest,
    ::testing::Combine(
        ::testing::Values(5, 0.1, -0.7, 0),
        ::testing::Values(true, false)
    )
);

/**********************************
 ********* SPECIAL VALUES *********
 **********************************/

TEST(ArithScalarTest, NonIeee754Multiply) {
    int scalar = 5;
    auto op = tatami::make_DelayedMultiplyScalarHelper(scalar);
    EXPECT_TRUE(op.is_sparse());
}

TEST(ArithScalarTest, NonFiniteMultiply) {
    double scalar = std::numeric_limits<double>::infinity();
    auto op = tatami::make_DelayedMultiplyScalarHelper(scalar);
    EXPECT_FALSE(op.is_sparse());
}

TEST(ArithScalarTest, NonIeee754Divide) {
    auto op = tatami::make_DelayedDivideScalarHelper<false, int>(5.0);
    tatami_test::throws_error([&]() { op.template fill<double>(5); }, "IEEE-754");
}
