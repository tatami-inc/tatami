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
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
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
    auto op = tatami::make_DelayedUnaryIsometricAddScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod->is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected.
    auto refvec = simulated;
    for (auto& r : refvec) {
        r += val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(DelayedUnaryIsometricArithmeticCommutativeScalarTest, Multiply) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricMultiplyScalar(val);
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r *= val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(DelayedUnaryIsometricArithmeticCommutativeScalarTest, NewType) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricAddScalar(val);
    auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
    auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

    size_t nsize = simulated.size();
    auto frefvec = std::vector<float>(nsize);
    for (size_t i = 0; i < nsize; ++i) {
        frefvec[i] = simulated[i] + val;
    }

    tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));
    quick_test_all(dense_fmod.get(), &fref);
    quick_test_all(sparse_fmod.get(), &fref);
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
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedUnaryIsometricSubtractScalar<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricSubtractScalar<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_FALSE(sparse_mod->is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Divide) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedUnaryIsometricDivideScalar<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricDivideScalar<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ (val == 0));
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ (val == 0));
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Power) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    tatami::DelayedUnaryIsometricAbs op0;
    dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op0);
    sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op0);

    auto my_param = GetParam();
    double val = std::abs(std::get<0>(my_param));
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedUnaryIsometricPowerScalar<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense_mod, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse_mod, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricPowerScalar<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense_mod, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse_mod, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, Modulo) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedUnaryIsometricModuloScalar<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricModuloScalar<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    auto refvec = simulated;
    for (auto& r : refvec) {
        if (on_right) {
            r = std::fmod(r, val);
        } else {
            r = std::fmod(val, r);
        }
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ !(on_right && val));
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ !(on_right && val));
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, IntegerDivide) {
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);
    if (on_right) {
        auto op = tatami::make_DelayedUnaryIsometricIntegerDivideScalar<true>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricIntegerDivideScalar<false>(val);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    if (on_right && val) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ !(on_right && val));
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ !(on_right && val));
}

TEST_P(DelayedUnaryIsometricArithmeticNonCommutativeScalarTest, NewType) {
    auto my_param = GetParam();
    double val = std::get<0>(my_param);
    bool on_right = std::get<1>(my_param);

    std::shared_ptr<tatami::Matrix<float, int> > dense_fmod, sparse_fmod;
    if (on_right) {
        auto op = tatami::make_DelayedUnaryIsometricSubtractScalar<true>(val);
        dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricSubtractScalar<false>(val);
        dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);
    }

    size_t nsize = simulated.size();
    auto frefvec = std::vector<float>(nsize);
    for (size_t i = 0; i < nsize; ++i) {
        if (on_right) {
            frefvec[i] = simulated[i] - val;
        } else {
            frefvec[i] = val - simulated[i];
        }
    }

    tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));
    quick_test_all(dense_fmod.get(), &fref);
    quick_test_all(sparse_fmod.get(), &fref);
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

TEST(DelayedUnaryIsometricArithmeticScalar, NonIeee754Multiply) {
    int scalar = 5;
    auto op = tatami::make_DelayedUnaryIsometricMultiplyScalar(scalar);
    EXPECT_TRUE(op.is_sparse());
}

TEST(DelayedUnaryIsometricArithmeticScalar, NonFiniteMultiply) {
    double scalar = std::numeric_limits<double>::infinity();
    auto op = tatami::make_DelayedUnaryIsometricMultiplyScalar(scalar);
    EXPECT_FALSE(op.is_sparse());
}

TEST(DelayedUnaryIsometricArithmeticScalar, NonIeee754Divide) {
    auto op = tatami::make_DelayedUnaryIsometricDivideScalar<false, int>(5.0);
    tatami_test::throws_error([&]() { op.template fill<double>(true, 5); }, "IEEE-754");
}
