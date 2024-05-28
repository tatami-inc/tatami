#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/math_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricMathTest : public ::testing::Test {
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

TEST_F(DelayedUnaryIsometricMathTest, Abs) {
    tatami::DelayedUnaryIsometricAbs op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::abs(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sign) {
    {
        tatami::DelayedUnaryIsometricSign op;
        auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

        EXPECT_FALSE(dense_mod->is_sparse());
        EXPECT_TRUE(sparse_mod->is_sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            if (r < 0) {
                r = -1;
            } else if (r > 0) {
                r = 1;
            }
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

        quick_test_all(dense_mod.get(), &ref);
        quick_test_all(sparse_mod.get(), &ref);

        // Checking that it works for a different output type.
        {
            tatami::DelayedUnaryIsometricSign op;
            auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
            auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

            std::vector<float> frefvec(refvec.begin(), refvec.end());
            tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

            quick_test_all(dense_fmod.get(), &fref);
            quick_test_all(sparse_fmod.get(), &fref);
        }
    }

    // Throwing in some NaNs.
    {
        auto simulated_nan = simulated;
        simulated_nan[0] = std::numeric_limits<double>::quiet_NaN();
        auto dense_nan = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated_nan));
        auto sparse_nan = tatami::convert_to_compressed_sparse<false, double, int>(dense_nan.get()); // column major.

        tatami::DelayedUnaryIsometricSign op;
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation(dense_nan, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation(sparse_nan, op);

        auto refvec = simulated_nan;
        for (auto& r : refvec) {
            if (r < 0) {
                r = -1;
            } else if (r > 0) {
                r = 1;
            }
        }
        tatami::DenseRowMatrix<double, int> ref_nan(nrow, ncol, refvec);

        quick_test_all(dense_nan.get(), &ref_nan, /* has_nan = */ true);
        quick_test_all(sparse_nan.get(), &ref_nan, /* has_nan = */ true);

        // Checking that it works for a different output type that supports NaNs.
        {
            auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense_nan, op);
            auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse_nan, op);

            std::vector<float> frefvec(refvec.begin(), refvec.end());
            tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

            quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
            quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
        }

        // Checking that it works for a different output type that doesn't support NaNs.
        {
            auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<int>(dense_nan, op);
            auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<int>(sparse_nan, op);

            std::vector<int> frefvec(refvec.size());
            for (size_t i = 0; i < refvec.size(); ++i) {
                frefvec[i] = std::isnan(refvec[i]) ? 0 : refvec[i];
            }
            tatami::DenseRowMatrix<int, int> fref(nrow, ncol, std::move(frefvec));

            quick_test_all(dense_fmod.get(), &fref);
            quick_test_all(sparse_fmod.get(), &fref);
        }
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sqrt) {
    tatami::DelayedUnaryIsometricSqrt op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sqrt(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // Again, doing some light tests; we assume that we have IEEE floats so sqrt(<negative value>) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Log) {
    // Trying with the natural base.
    {
        tatami::DelayedUnaryIsometricLog op;
        auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

        EXPECT_FALSE(dense_mod->is_sparse());
        EXPECT_FALSE(sparse_mod->is_sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log(r);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

        // Doing some light tests, assuming that log(<negative value>) => NaN and log(0) => Inf.
        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }

    // Trying with another base.
    {
        tatami::DelayedUnaryIsometricLog op(2.0);
        auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

        EXPECT_FALSE(dense_mod->is_sparse());
        EXPECT_FALSE(sparse_mod->is_sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log(r) / std::log(2);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }

    // Checking that it works for a different output type.
    {
        tatami::DelayedUnaryIsometricLog op;
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(simulated.size());
        for (size_t i = 0; i < simulated.size(); ++i) {
            frefvec[i] = std::log(simulated[i]);
        }
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Log1p) {
    // Trying with the natural base.
    {
        tatami::DelayedUnaryIsometricLog1p op;
        auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

        EXPECT_FALSE(dense_mod->is_sparse());
        EXPECT_TRUE(sparse_mod->is_sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log1p(r);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

        // Doing some light tests, assuming that log1p(<less than 1>) => NaN and log1p(-1) => Inf.
        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }

    // Trying with another base.
    {
        tatami::DelayedUnaryIsometricLog1p op(2.0);
        auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

        EXPECT_FALSE(dense_mod->is_sparse());
        EXPECT_TRUE(sparse_mod->is_sparse());
        EXPECT_EQ(dense->nrow(), dense_mod->nrow());
        EXPECT_EQ(dense->ncol(), dense_mod->ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log1p(r) / std::log(2);
        }
        tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

        quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
        quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);
    }

    // Checking that it works for a different output type.
    {
        tatami::DelayedUnaryIsometricLog1p op;
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(simulated.size());
        for (size_t i = 0; i < simulated.size(); ++i) {
            frefvec[i] = std::log1p(simulated[i]);
        }
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Exp) {
    tatami::DelayedUnaryIsometricExp op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::exp(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Expm1) {
    tatami::DelayedUnaryIsometricExpm1 op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::expm1(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Round) {
    tatami::DelayedUnaryIsometricRound op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::round(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Ceiling) {
    tatami::DelayedUnaryIsometricCeiling op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::ceil(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Floor) {
    tatami::DelayedUnaryIsometricFloor op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::floor(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Trunc) {
    tatami::DelayedUnaryIsometricTrunc op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::trunc(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sin) {
    tatami::DelayedUnaryIsometricSin op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sin(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Cos) {
    tatami::DelayedUnaryIsometricCos op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::cos(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Tan) {
    tatami::DelayedUnaryIsometricTan op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tan(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Asin) {
    // Use a tighter range to get most values inside the domain of [-1, 1].
    tatami::DelayedUnaryIsometricAsin op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense_unit, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse_unit, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::asin(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // Assume assume that asin(<below -1 or above 1>) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense_unit, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse_unit, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Acos) {
    // Use a tighter range to get most values inside the domain of [-1, 1].
    tatami::DelayedUnaryIsometricAcos op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense_unit, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse_unit, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::acos(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // We assume that acos(<below -1 or above 1>) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense_unit, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse_unit, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Atan) {
    tatami::DelayedUnaryIsometricAtan op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::atan(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sinh) {
    tatami::DelayedUnaryIsometricSinh op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sinh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Cosh) {
    tatami::DelayedUnaryIsometricCosh op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::cosh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Tanh) {
    tatami::DelayedUnaryIsometricTanh op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tanh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Asinh) {
    tatami::DelayedUnaryIsometricAsinh op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::asinh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Acosh) {
    tatami::DelayedUnaryIsometricAcosh op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::acosh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // We assume that acosh(<less than 1>) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Atanh) {
    tatami::DelayedUnaryIsometricAtanh op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense_unit, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse_unit, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::atanh(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // Again, doing some light tests. We assume that atanh(<below -1 or above 1>) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense_unit, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse_unit, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Gamma) {
    tatami::DelayedUnaryIsometricGamma op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tgamma(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    // We assume that gamma(<less than or equal to zero>) => NaN.
    quick_test_all(dense_mod.get(), &ref, /* has_nan = */ true);
    quick_test_all(sparse_mod.get(), &ref, /* has_nan = */ true);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref, /* has_nan = */ true);
        quick_test_all(sparse_fmod.get(), &fref, /* has_nan = */ true);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Lgamma) {
    tatami::DelayedUnaryIsometricLgamma op;
    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::lgamma(r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, refvec);

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);

    // Checking that it works for a different output type.
    {
        auto dense_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(dense, op);
        auto sparse_fmod = tatami::make_DelayedUnaryIsometricOperation<float>(sparse, op);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseRowMatrix<float, int> fref(nrow, ncol, std::move(frefvec));

        quick_test_all(dense_fmod.get(), &fref);
        quick_test_all(sparse_fmod.get(), &fref);
    }
}
