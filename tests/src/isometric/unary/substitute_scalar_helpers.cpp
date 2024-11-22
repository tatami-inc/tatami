#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOperation.hpp"
#include "tatami/isometric/unary/substitute_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedUnaryIsometricSubstituteScalarTest : public ::testing::TestWithParam<std::tuple<double, double> > {
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
            opt.seed = 128736123;
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

TEST_P(DelayedUnaryIsometricSubstituteScalarTest, Equal) {
    auto params = GetParam();
    auto comp = std::get<0>(params);
    auto sub = std::get<1>(params);
    auto op = tatami::make_DelayedUnaryIsometricSubstituteEqualScalar(comp, sub);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (comp != 0 || sub == 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (r == comp ? sub : r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteScalarTest, GreaterThan) {
    auto params = GetParam();
    auto comp = std::get<0>(params);
    auto sub = std::get<1>(params);
    auto op = tatami::make_DelayedUnaryIsometricSubstituteGreaterThanScalar(comp, sub);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (comp >= 0 || sub == 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (r > comp ? sub : r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteScalarTest, LessThan) {
    auto params = GetParam();
    auto comp = std::get<0>(params);
    auto sub = std::get<1>(params);
    auto op = tatami::make_DelayedUnaryIsometricSubstituteLessThanScalar(comp, sub);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (comp <= 0 || sub == 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (r < comp ? sub : r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteScalarTest, GreaterThanOrEqual) {
    auto params = GetParam();
    auto comp = std::get<0>(params);
    auto sub = std::get<1>(params);
    auto op = tatami::make_DelayedUnaryIsometricSubstituteGreaterThanOrEqualScalar(comp, sub);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (comp > 0 || sub == 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (r >= comp ? sub : r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteScalarTest, LessThanOrEqual) {
    auto params = GetParam();
    auto comp = std::get<0>(params);
    auto sub = std::get<1>(params);
    auto op = tatami::make_DelayedUnaryIsometricSubstituteLessThanOrEqualScalar(comp, sub);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (comp < 0 || sub == 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (r <= comp ? sub : r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSubstituteScalarTest, NotEqual) {
    auto params = GetParam();
    auto comp = std::get<0>(params);
    auto sub = std::get<1>(params);
    auto op = tatami::make_DelayedUnaryIsometricSubstituteNotEqualScalar(comp, sub);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (comp == 0 || sub == 0) {
        EXPECT_TRUE(sparse_mod->is_sparse());
    } else {
        EXPECT_FALSE(sparse_mod->is_sparse());
    }

    // Toughest tests are handled by 'arith_vector.hpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = (r != comp ? sub : r);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricSubstituteScalar,
    DelayedUnaryIsometricSubstituteScalarTest,
    ::testing::Combine(
        ::testing::Values(0.0, -1.0, 1.0),
        ::testing::Values(0.0, 2.0)
    )
);

/********************************************
 ********************************************/

class DelayedUnaryIsometricSpecialSubstituteTest : public ::testing::TestWithParam<bool> {
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

TEST_P(DelayedUnaryIsometricSpecialSubstituteTest, NaN) {
    bool pass = GetParam();
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricSubstituteIsnan(69.0);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricSubstituteIsnan<false>(69.0);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);
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
        if ((pass && std::isnan(r)) || (!pass && !std::isnan(r))) {
            r = 69;
        }
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialSubstituteTest, Infinity) {
    bool pass = GetParam();
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricSubstituteIsinf(69.0);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricSubstituteIsinf<false>(69.0);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);
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
        if ((pass && std::isinf(r)) || (!pass && !std::isinf(r))) {
            r = 69;
        }
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_P(DelayedUnaryIsometricSpecialSubstituteTest, Finite) {
    bool pass = GetParam();
    std::shared_ptr<tatami::NumericMatrix> dense_mod, sparse_mod;

    if (pass) {
        auto op = tatami::make_DelayedUnaryIsometricSubstituteIsfinite(69.0);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    } else {
        auto op = tatami::make_DelayedUnaryIsometricSubstituteIsfinite<false>(69.0);
        dense_mod = tatami::make_DelayedUnaryIsometricOperation(dense, op);
        sparse_mod = tatami::make_DelayedUnaryIsometricOperation(sparse, op);
    }

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);
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
        if ((pass && std::isfinite(r)) || (!pass && !std::isfinite(r))) {
            r = 69;
        }
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnaryIsometricSpecialSubstitute,
    DelayedUnaryIsometricSpecialSubstituteTest,
    ::testing::Values(true, false)
);

