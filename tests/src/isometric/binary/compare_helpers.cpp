#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class BinaryCompareTest : public ::testing::Test {
protected:
    inline static size_t nrow = 151, ncol = 71;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    inline static std::vector<double> simulated_left, simulated_right;

    static void SetUpTestSuite() {
        simulated_left = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.2, /* lower = */ 1, /* upper = */ 4, /* seed */ 12345);
        for (auto& x : simulated_left) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_left = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated_left));
        sparse_left = tatami::convert_to_compressed_sparse<false, double, int>(dense_left.get()); // column major.

        simulated_right = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.2, /* lower = */ 1, /* upper = */ 4, /* seed */ 67890);
        for (auto& x : simulated_right) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_right = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated_right));
        sparse_right = tatami::convert_to_compressed_sparse<false, double, int>(dense_right.get()); // column major.
        return;
    }
};

TEST_F(BinaryCompareTest, Equal) {
    auto op = tatami::make_DelayedBinaryEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] == simulated_right[i]);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, GreaterThan) {
    auto op = tatami::make_DelayedBinaryGreaterThanHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] > simulated_right[i]);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, LessThan) {
    auto op = tatami::make_DelayedBinaryLessThanHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] < simulated_right[i]);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, GreaterThanOrEqual) {
    auto op = tatami::make_DelayedBinaryGreaterThanOrEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] >= simulated_right[i]);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, LessThanOrEqual) {
    auto op = tatami::make_DelayedBinaryLessThanOrEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] <= simulated_right[i]);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}

TEST_F(BinaryCompareTest, NotEqual) {
    auto op = tatami::make_DelayedBinaryNotEqualHelper();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOp(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOp(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] != simulated_right[i]);
    }
    tatami::DenseRowMatrix<double, int> ref(nrow, ncol, std::move(refvec));

    quick_test_all(dense_mod.get(), &ref);
    quick_test_all(sparse_mod.get(), &ref);
}
