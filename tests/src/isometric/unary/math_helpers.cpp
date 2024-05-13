#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class MathTest : public ::testing::Test {
protected:
    inline static size_t nrow = 82, ncol = 51;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    inline static std::shared_ptr<tatami::NumericMatrix> dense_unit, sparse_unit;
    inline static std::vector<double> simulated_unit;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, /* lower = */ -10, /* upper = */ 10);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.

        // Use a tighter range to get most values inside the domain of [-1, 1]
        // (but not all; we still leave a few outside for NaN testing purposes).
        simulated_unit = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, /* lower = */ -1.1, /* upper = */ 1.1);
        dense_unit = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated_unit));
        sparse_unit = tatami::convert_to_compressed_sparse<false, double, int>(dense_unit.get()); // column major.
    }
};

TEST_F(MathTest, Abs) {
    tatami::DelayedAbsHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::abs(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Sign) {
    tatami::DelayedSignHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (0 < r) - (r < 0);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Sqrt) {
    tatami::DelayedSqrtHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sqrt(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests; we assume that we have IEEE floats so sqrt(-1) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}

TEST_F(MathTest, Log) {
    // Trying with the natural base.
    {
        tatami::DelayedLogHelper op;
        auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        EXPECT_FALSE(dense_mod->sparse());
        EXPECT_FALSE(sparse_mod->sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log(r);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

        // Doing some light tests, assuming that log(-1) => NaN and log(0) => Inf.
        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }

    // Trying with another base.
    {
        tatami::DelayedLogHelper op(2.0);
        auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        EXPECT_FALSE(dense_mod->sparse());
        EXPECT_FALSE(sparse_mod->sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log(r) / std::log(2);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

        // Again, doing some light tests.
        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }
}

TEST_F(MathTest, Log1pBy) {
    // Trying with the natural base.
    {
        tatami::DelayedLog1pHelper op;
        auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        EXPECT_FALSE(dense_mod->sparse());
        EXPECT_TRUE(sparse_mod->sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log1p(r);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

        // Doing some light tests, assuming that log1p(-2) => NaN and log1p(-1) => Inf.
        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }

    // Trying with another base.
    {
        tatami::DelayedLog1pHelper op(2.0);
        auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

        EXPECT_FALSE(dense_mod->sparse());
        EXPECT_TRUE(sparse_mod->sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log1p(r) / std::log(2);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

        // Again, doing some light tests.
        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }
}

TEST_F(MathTest, Exp) {
    tatami::DelayedExpHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::exp(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Expm1) {
    tatami::DelayedExpm1Helper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::expm1(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Round) {
    tatami::DelayedRoundHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::round(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Ceiling) {
    tatami::DelayedCeilingHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::ceil(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Floor) {
    tatami::DelayedFloorHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::floor(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Trunc) {
    tatami::DelayedTruncHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::trunc(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Sin) {
    tatami::DelayedSinHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sin(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Cos) {
    tatami::DelayedCosHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::cos(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Tan) {
    tatami::DelayedTanHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tan(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Asin) {
    // Use a tighter range to get most values inside the domain of [-1, 1].
    tatami::DelayedAsinHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense_unit, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse_unit, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::asin(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests. We assume that asin(2) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}

TEST_F(MathTest, Acos) {
    tatami::DelayedAcosHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense_unit, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse_unit, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::acos(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests. We assume that acos(-2) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}

TEST_F(MathTest, Atan) {
    tatami::DelayedAtanHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::atan(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Sinh) {
    tatami::DelayedSinhHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sinh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Cosh) {
    tatami::DelayedCoshHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::cosh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Tanh) {
    tatami::DelayedTanhHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tanh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Asinh) {
    tatami::DelayedAsinhHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::asinh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests.
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(MathTest, Acosh) {
    tatami::DelayedAcoshHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::acosh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests. We assume that acosh(-1) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}

TEST_F(MathTest, Atanh) {
    tatami::DelayedAtanhHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense_unit, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse_unit, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::atanh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests. We assume that atanh(2) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}

TEST_F(MathTest, Gamma) {
    tatami::DelayedGammaHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tgamma(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests. We assume that gamma(-1) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}

TEST_F(MathTest, Lgamma) {
    tatami::DelayedLgammaHelper op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::lgamma(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    // Again, doing some light tests. We assume that lgamma(-1) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
}
