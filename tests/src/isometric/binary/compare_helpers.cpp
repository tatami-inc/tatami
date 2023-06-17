#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class BinaryCompareTest : public ::testing::Test {
protected:
    size_t nrow = 151, ncol = 71;
    std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    std::vector<double> simulated_left, simulated_right;
protected:
    void SetUp() {
        simulated_left = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.2, /* lower = */ 1, /* upper = */ 4, /* seed */ 12345);
        for (auto& x : simulated_left) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_left = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated_left));
        sparse_left = tatami::convert_to_sparse<false>(dense_left.get()); // column major.

        simulated_right = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.2, /* lower = */ 1, /* upper = */ 4, /* seed */ 67890);
        for (auto& x : simulated_right) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_right = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated_right));
        sparse_right = tatami::convert_to_sparse<false>(dense_right.get()); // column major.
        return;
    }
};

TEST_F(BinaryCompareTest, Equal) {
    auto op = tatami::make_DelayedBinaryEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] == simulated_right[i]);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, GreaterThan) {
    auto op = tatami::make_DelayedBinaryGreaterThanHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] > simulated_right[i]);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, LessThan) {
    auto op = tatami::make_DelayedBinaryLessThanHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] < simulated_right[i]);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, GreaterThanOrEqual) {
    auto op = tatami::make_DelayedBinaryGreaterThanOrEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] >= simulated_right[i]);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, LessThanOrEqual) {
    auto op = tatami::make_DelayedBinaryLessThanOrEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] <= simulated_right[i]);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, NotEqual) {
    auto op = tatami::make_DelayedBinaryNotEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    // Toughest tests are handled by the Vector case; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the scalar operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] != simulated_right[i]);
    }
    tatami::DenseRowMatrix<double> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}
