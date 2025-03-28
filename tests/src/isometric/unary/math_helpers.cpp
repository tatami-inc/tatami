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
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = 4239867123;
            return opt;
        }());
        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.

        // Use a tighter range to get most values inside the domain of [-1, 1]
        // (but not all; we still leave a few outside for NaN testing purposes).
        simulated_unit = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.lower = -1.1;
            opt.upper = 1.1;
            opt.seed = 9234723;
            return opt;
        }());
        dense_unit.reset(new tatami::DenseMatrix<double, int, decltype(simulated_unit)>(nrow, ncol, simulated_unit, true)); // row major.
        sparse_unit = tatami::convert_to_compressed_sparse<double, int>(*dense_unit, false, {}); // column major.
    }
};

TEST_F(DelayedUnaryIsometricMathTest, Abs) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAbsHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::abs(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAbsHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sign) {
    {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricSignHelper<double, double, int> >();
        tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
        tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

        EXPECT_FALSE(dense_mod.is_sparse());
        EXPECT_TRUE(sparse_mod.is_sparse());
        EXPECT_EQ(nrow, dense_mod.nrow());
        EXPECT_EQ(ncol, dense_mod.ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            if (r < 0) {
                r = -1;
            } else if (r > 0) {
                r = 1;
            }
        }
        tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

        quick_test_all<double, int>(dense_mod, ref);
        quick_test_all<double, int>(sparse_mod, ref);

        // Checking that it works for a different output type.
        {
            auto fop = std::make_shared<tatami::DelayedUnaryIsometricSign<float, double, int> >();
            tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
            tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

            std::vector<float> frefvec(refvec.begin(), refvec.end());
            tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

            quick_test_all<float, int>(dense_fmod, fref);
            quick_test_all<float, int>(sparse_fmod, fref);
        }
    }

    // Throwing in some NaNs.
    {
        auto simulated_nan = simulated;
        simulated_nan[0] = std::numeric_limits<double>::quiet_NaN();
        auto dense_nan = std::make_shared<tatami::DenseMatrix<double, int, decltype(simulated_nan)> >(nrow, ncol, simulated_nan, true); // row major.
        auto sparse_nan = tatami::convert_to_compressed_sparse<double, int>(*dense_nan, false, {}); // column major.

        auto op = std::make_shared<tatami::DelayedUnaryIsometricSignHelper<double, double, int> >();
        tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense_nan, op);
        tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse_nan, op);

        auto refvec = simulated_nan;
        for (auto& r : refvec) {
            if (r < 0) {
                r = -1;
            } else if (r > 0) {
                r = 1;
            }
        }
        tatami::DenseMatrix<double, int, decltype(refvec)> ref_nan(nrow, ncol, refvec, true);

        quick_test_all<double, int>(dense_mod, ref_nan);
        quick_test_all<double, int>(sparse_mod, ref_nan);

        // Checking that it works for a different output type that supports NaNs.
        {
            auto fop = std::make_shared<tatami::DelayedUnaryIsometricSignHelper<float, double, int> >();
            tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense_nan, fop);
            tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse_nan, fop);

            std::vector<float> frefvec(refvec.begin(), refvec.end());
            tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

            quick_test_all<float, int>(dense_fmod, fref);
            quick_test_all<float, int>(sparse_fmod, fref);
        }

        // Checking that it works for a different output type that doesn't support NaNs.
        {
            auto iop = std::make_shared<tatami::DelayedUnaryIsometricSignHelper<int, double, int> >();
            tatami::DelayedUnaryIsometricOperation<int, double, int> dense_imod(dense_nan, iop);
            tatami::DelayedUnaryIsometricOperation<int, double, int> sparse_imod(sparse_nan, iop);

            std::vector<int> irefvec(refvec.size());
            for (size_t i = 0; i < refvec.size(); ++i) {
                irefvec[i] = std::isnan(refvec[i]) ? 0 : refvec[i];
            }
            tatami::DenseMatrix<int, int, decltype(irefvec)> iref(nrow, ncol, std::move(irefvec), true);

            quick_test_all<int, int>(dense_imod, iref);
            quick_test_all<int, int>(sparse_imod, iref);
        }
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sqrt) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricSqrtHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sqrt(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // Again, doing some light tests; we assume that we have IEEE floats so sqrt(<negative value>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricSqrtHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Log) {
    // Trying with the natural base.
    {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricLogHelper<double, double, int, double> >();
        tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
        tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

        EXPECT_FALSE(dense_mod.is_sparse());
        EXPECT_FALSE(sparse_mod.is_sparse());
        EXPECT_EQ(nrow, dense_mod.nrow());
        EXPECT_EQ(ncol, dense_mod.ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log(r);
        }
        tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

        // Doing some light tests, assuming that log(<negative value>) => NaN and log(0) => Inf.
        quick_test_all<double, int>(dense_mod, ref); 
        quick_test_all<double, int>(sparse_mod, ref);
    }

    // Trying with another base.
    {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricLogHelper<double, double, int, double> >(2);
        tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
        tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

        EXPECT_FALSE(dense_mod.is_sparse());
        EXPECT_FALSE(sparse_mod.is_sparse());
        EXPECT_EQ(nrow, dense_mod.nrow());
        EXPECT_EQ(ncol, dense_mod.ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log(r) / std::log(2);
        }
        tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

        quick_test_all<double, int>(dense_mod, ref);
        quick_test_all<double, int>(sparse_mod, ref);
    }

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricLogHelper<float, double, int, double> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(simulated.size());
        for (size_t i = 0; i < simulated.size(); ++i) {
            frefvec[i] = std::log(simulated[i]);
        }
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Log1p) {
    // Trying with the natural base.
    {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricLog1pHelper<double, double, int, double> >();
        tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
        tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

        EXPECT_FALSE(dense_mod.is_sparse());
        EXPECT_TRUE(sparse_mod.is_sparse());
        EXPECT_EQ(nrow, dense_mod.nrow());
        EXPECT_EQ(ncol, dense_mod.ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log1p(r);
        }
        tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

        // Doing some light tests, assuming that log1p(<less than 1>) => NaN and log1p(-1) => Inf.
        quick_test_all<double, int>(dense_mod, ref);
        quick_test_all<double, int>(sparse_mod, ref);
    }

    // Trying with another base.
    {
        auto op = std::make_shared<tatami::DelayedUnaryIsometricLog1pHelper<double, double, int, double> >(2.0);
        tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
        tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

        EXPECT_FALSE(dense_mod.is_sparse());
        EXPECT_TRUE(sparse_mod.is_sparse());
        EXPECT_EQ(nrow, dense_mod.nrow());
        EXPECT_EQ(ncol, dense_mod.ncol());

        auto refvec = simulated;
        for (auto& r : refvec) {
            r = std::log1p(r) / std::log(2);
        }
        tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

        quick_test_all<double, int>(dense_mod, ref);
        quick_test_all<double, int>(sparse_mod, ref);
    }

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricLog1pHelper<float, double, int, double> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(simulated.size());
        for (size_t i = 0; i < simulated.size(); ++i) {
            frefvec[i] = std::log1p(simulated[i]);
        }
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Exp) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricExpHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::exp(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricExpHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Expm1) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricExpm1Helper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::expm1(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricExpm1Helper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Round) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricRoundHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::round(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricRoundHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Ceiling) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricCeilingHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::ceil(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricCeilingHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Floor) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricFloorHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::floor(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricFloorHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Trunc) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricTruncHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::trunc(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricTruncHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sin) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricSinHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sin(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricSinHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Cos) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricCosHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::cos(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricCosHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Tan) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricTanHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tan(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricTanHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Asin) {
    // Use a tighter range to get most values inside the domain of [-1, 1].
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAsinHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense_unit, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse_unit, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::asin(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // Assume assume that asin(<below -1 or above 1>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAsinHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense_unit, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse_unit, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Acos) {
    // Use a tighter range to get most values inside the domain of [-1, 1].
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAcosHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense_unit, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse_unit, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated_unit;
    for (auto& r : refvec) {
        r = std::acos(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // We assume that acos(<below -1 or above 1>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAcosHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense_unit, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse_unit, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Atan) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAtanHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::atan(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAtanHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Sinh) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricSinhHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::sinh(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricSinhHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Cosh) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricCoshHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::cosh(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricCoshHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Tanh) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricTanhHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tanh(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricTanhHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Asinh) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAsinhHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::asinh(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAsinhHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Acosh) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAcoshHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::acosh(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // We assume that acosh(<less than 1>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAcoshHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Atanh) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricAtanhHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::atanh(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // Again, doing some light tests. We assume that atanh(<below -1 or above 1>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricAtanhHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Gamma) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricGammaHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::tgamma(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // We assume that gamma(<less than or equal to zero>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricGammaHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}

TEST_F(DelayedUnaryIsometricMathTest, Lgamma) {
    auto op = std::make_shared<tatami::DelayedUnaryIsometricLgammaHelper<double, double, int> >();
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(nrow, dense_mod.nrow());
    EXPECT_EQ(ncol, dense_mod.ncol());

    auto refvec = simulated;
    for (auto& r : refvec) {
        r = std::lgamma(r);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, refvec, true);

    // We assume that gamma(<less than or equal to zero>) => NaN.
    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // Checking that it works for a different output type.
    {
        auto fop = std::make_shared<tatami::DelayedUnaryIsometricLgammaHelper<float, double, int> >();
        tatami::DelayedUnaryIsometricOperation<float, double, int> dense_fmod(dense, fop);
        tatami::DelayedUnaryIsometricOperation<float, double, int> sparse_fmod(sparse, fop);

        std::vector<float> frefvec(refvec.begin(), refvec.end());
        tatami::DenseMatrix<float, int, decltype(frefvec)> fref(nrow, ncol, std::move(frefvec), true);

        quick_test_all<float, int>(dense_fmod, fref);
        quick_test_all<float, int>(sparse_fmod, fref);
    }
}
