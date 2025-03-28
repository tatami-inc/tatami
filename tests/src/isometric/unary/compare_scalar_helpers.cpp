#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/compare_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricCompareScalarTest : public ::testing::TestWithParam<int> {
protected:
    inline static size_t nrow = 123, ncol = 89;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.1;
            opt.lower = -2;
            opt.upper = 2;
            opt.seed = 1867423;
            return opt;
        }());

        for (auto& x : simulated) {
            if (x) {
                // Rounding for easier tests of exact equality.
                x = std::round(x);
                if (x == 0) {
                    x = 1;
                }
            }
        }

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.
    }
};

TEST_P(DelayedUnaryIsometricCompareScalarTest, Equal) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricEqualScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (val) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r == val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, GreaterThan) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricGreaterThanScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (val >= 0) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r > val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, LessThan) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricLessThanScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (val <= 0) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r < val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, GreaterThanOrEqual) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricGreaterThanOrEqualScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (val > 0) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r >= val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, LessThanOrEqual) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricLessThanOrEqualScalarHelper<double, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);

    if (val < 0) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r <= val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, NotEqual) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricNotEqualScalarHelper<double, double, int, double> >(val);
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

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r != val;
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, NewType) {
    double val = GetParam();
    auto op = std::make_shared<tatami::DelayedUnaryIsometricNotEqualScalarHelper<uint8_t, double, int, double> >(val);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, op);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, op);

    EXPECT_FALSE(dense_umod.is_sparse());
    if (val) {
        EXPECT_FALSE(sparse_umod.is_sparse());
    } else {
        EXPECT_TRUE(sparse_umod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) { 
        urefvec[i] = simulated[i] != val;
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricCompareScalar,
    DelayedUnaryIsometricCompareScalarTest,
    ::testing::Values(0, -1, 1)
);

/********************************************
 ********************************************/

class DelayedUnaryIsometricSpecialCompareTest : public ::testing::TestWithParam<bool> {
protected:
    inline static size_t nrow = 123, ncol = 89;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    inline static std::vector<double> simulated;

    static void SetUpTestSuite() {
        std::mt19937_64 rng(nrow * ncol);
        std::uniform_real_distribution dist;

        simulated.resize(nrow * ncol);
        for (size_t i = 0; i < simulated.size(); ++i) {
            double val = dist(rng);
            if (val < 0.05) {
                simulated[i] = std::numeric_limits<double>::quiet_NaN();
            } else if (val < 0.1) {
                simulated[i] = std::numeric_limits<double>::infinity();
            } else if (val < 0.15) {
                simulated[i] = -std::numeric_limits<double>::infinity();
            } else if (val < 0.2) {
                simulated[i] = val;
            }
        }

        dense.reset(new tatami::DenseMatrix<double, int, decltype(simulated)>(nrow, ncol, simulated, true)); // row major.
        sparse = tatami::convert_to_compressed_sparse<double, int>(*dense, false, {}); // column major.
    }
};

TEST_P(DelayedUnaryIsometricSpecialCompareTest, NaN) {
    bool pass = GetParam();
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (pass) {
        op.reset(new tatami::DelayedUnaryIsometricIsnanHelper<true, double, double, int>);
    } else {
        op.reset(new tatami::DelayedUnaryIsometricIsnanHelper<false, double, double, int>);
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);
    if (pass) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        if (pass) {
            r = std::isnan(r);
        } else {
            r = !std::isnan(r);
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialCompareTest, Infinity) {
    bool pass = GetParam();
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (pass) {
        op.reset(new tatami::DelayedUnaryIsometricIsinfHelper<true, double, double, int>);
    } else {
        op.reset(new tatami::DelayedUnaryIsometricIsinfHelper<false, double, double, int>);
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);
    if (pass) {
        EXPECT_TRUE(sparse_mod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        if (pass) {
            r = std::isinf(r);
        } else {
            r = !std::isinf(r);
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialCompareTest, Finite) {
    bool pass = GetParam();
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > op;
    if (pass) {
        op.reset(new tatami::DelayedUnaryIsometricIsfiniteHelper<true, double, double, int>);
    } else {
        op.reset(new tatami::DelayedUnaryIsometricIsfiniteHelper<false, double, double, int>);
    }
    tatami::DelayedUnaryIsometricOperation<double, double, int> dense_mod(dense, op);
    tatami::DelayedUnaryIsometricOperation<double, double, int> sparse_mod(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_EQ(dense_mod.nrow(), nrow);
    EXPECT_EQ(dense_mod.ncol(), ncol);
    if (pass) {
        EXPECT_FALSE(sparse_mod.is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        if (pass) {
            r = std::isfinite(r);
        } else {
            r = !std::isfinite(r);
        }
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialCompareTest, NewType) {
    bool pass = GetParam();
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<uint8_t, double, int> > op;
    if (pass) {
        op.reset(new tatami::DelayedUnaryIsometricIsnanHelper<true, uint8_t, double, int>);
    } else {
        op.reset(new tatami::DelayedUnaryIsometricIsnanHelper<false, uint8_t, double, int>);
    }
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> dense_umod(dense, op);
    tatami::DelayedUnaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse, op);

    EXPECT_FALSE(dense_umod.is_sparse());
    if (pass) {
        EXPECT_TRUE(sparse_umod.is_sparse());
    } else {
        EXPECT_FALSE(sparse_umod.is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) {
        if (pass) {
            urefvec[i] = std::isnan(simulated[i]);
        } else {
            urefvec[i] = !std::isnan(simulated[i]);
        }
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricSpecialCompare,
    DelayedUnaryIsometricSpecialCompareTest,
    ::testing::Values(true, false)
);
