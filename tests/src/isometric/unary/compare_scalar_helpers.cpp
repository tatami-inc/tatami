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

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
    }
};

TEST_P(DelayedUnaryIsometricCompareScalarTest, Equal) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricEqualScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);

    if (val) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r == val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, GreaterThan) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricGreaterThanScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);

    if (val >= 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r > val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, LessThan) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricLessThanScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);

    if (val <= 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r < val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, GreaterThanOrEqual) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricGreaterThanOrEqualScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);

    if (val > 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r >= val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, LessThanOrEqual) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricLessThanOrEqualScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);

    if (val < 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r <= val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, NotEqual) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricNotEqualScalar(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);

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
        r = r != val;
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricCompareScalarTest, NewType) {
    double val = GetParam();
    auto op = tatami::make_DelayedUnaryIsometricNotEqualScalar(val);

    auto dense_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, op);
    auto sparse_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, op);

    EXPECT_FALSE(dense_umod->is_sparse());
    if (val) {
        EXPECT_FALSE(sparse_umod->is_sparse());
    } else {
        EXPECT_TRUE(sparse_umod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<uint8_t> urefvec(simulated.size());
    for (size_t i = 0; i < simulated.size(); ++i) { 
        urefvec[i] = simulated[i] != val;
    }
    tatami::DenseRowMatrix<uint8_t, int> uref(nrow, ncol, std::move(urefvec));

    quick_test_all<uint8_t, int>(*dense_umod, uref);
    quick_test_all<uint8_t, int>(*sparse_umod, uref);
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

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); // column major.
    }
};

TEST_P(DelayedUnaryIsometricSpecialCompareTest, NaN) {
    bool pass = GetParam();
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricIsnan();
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricIsnan<false>();
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);
    if (pass) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialCompareTest, Infinity) {
    bool pass = GetParam();
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricIsinf();
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricIsinf<false>();
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);
    if (pass) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialCompareTest, Finite) {
    bool pass = GetParam();
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricIsfinite();
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricIsfinite<false>();
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense_mod->nrow(), nrow);
    EXPECT_EQ(dense_mod->ncol(), ncol);
    if (pass) {
        EXPECT_FALSE(sparse_mod->is_sparse());
    } else {
        EXPECT_TRUE(sparse_mod->is_sparse());
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
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialCompareTest, NewType) {
    bool pass = GetParam();
    std::shared_ptr<tatami::Matrix<uint8_t, int> > dense_umod, sparse_umod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricIsnan();
        dense_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, op);
        sparse_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricIsnan<false>();
        dense_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(dense, op);
        sparse_umod = tatami::make_DelayedUnaryIsometricOperation<uint8_t>(sparse, op);
    }

    EXPECT_FALSE(dense_umod->is_sparse());
    if (pass) {
        EXPECT_TRUE(sparse_umod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_umod->is_sparse());
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
    tatami::DenseRowMatrix<uint8_t, int> uref(nrow, ncol, std::move(urefvec));

    quick_test_all<uint8_t, int>(*dense_umod, uref);
    quick_test_all<uint8_t, int>(*sparse_umod, uref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricSpecialCompare,
    DelayedUnaryIsometricSpecialCompareTest,
    ::testing::Values(true, false)
);
