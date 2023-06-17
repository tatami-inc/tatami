#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class CompareScalarTest : public ::testing::TestWithParam<int> { 
protected:
    size_t nrow = 123, ncol = 89;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse;
    std::vector<double> simulated;
protected:
    void SetUp() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1, -2, 2);
        for (auto& x : simulated) {
            if (x) {
                // Rounding for easier tests of exact equality.
                x = std::round(x);
                if (x == 0) {
                    x = 1;
                }
            }
        }

        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column major.
        return;
    }
};

TEST_P(CompareScalarTest, Equal) {
    double val = GetParam();
    auto op = tatami::make_DelayedEqualScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r == val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(CompareScalarTest, GreaterThan) {
    double val = GetParam();
    auto op = tatami::make_DelayedGreaterThanScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val >= 0) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r > val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(CompareScalarTest, LessThan) {
    double val = GetParam();
    auto op = tatami::make_DelayedLessThanScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val <= 0) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r < val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(CompareScalarTest, GreaterThanOrEqual) {
    double val = GetParam();
    auto op = tatami::make_DelayedGreaterThanOrEqualScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val > 0) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r >= val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(CompareScalarTest, LessThanOrEqual) {
    double val = GetParam();
    auto op = tatami::make_DelayedLessThanOrEqualScalarHelper(val);

    auto dense_mod = tatami::make_DelayedUnaryIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedUnaryIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_EQ(dense->nrow(), nrow);
    EXPECT_EQ(dense->ncol(), ncol);

    if (val < 0) {
        EXPECT_TRUE(sparse_mod->sparse());
    } else {
        EXPECT_FALSE(sparse_mod->sparse());
    }

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r <= val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_P(CompareScalarTest, NotEqual) {
    double val = GetParam();
    auto op = tatami::make_DelayedNotEqualScalarHelper(val);

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

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    auto refvec = simulated;
    for (auto& r : refvec) {
        r = r != val;
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

INSTANTIATE_TEST_CASE_P(
    CompareScalar,
    CompareScalarTest,
    ::testing::Values(0, -1, 1)
);
