#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "tatami/isometric/binary/compare_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedBinaryIsometricCompareTest : public ::testing::Test {
protected:
    inline static size_t nrow = 151, ncol = 71;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    inline static std::vector<double> simulated_left, simulated_right;

    static void SetUpTestSuite() {
        simulated_left = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.18;
            opt.lower = 1;
            opt.upper = 4;
            opt.seed = 13579;
            return opt;
        }());
        for (auto& x : simulated_left) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_left.reset(new tatami::DenseMatrix<double, int, decltype(simulated_left)>(nrow, ncol, simulated_left, true)); // row major
        sparse_left = tatami::convert_to_compressed_sparse<double, int>(*dense_left, false, {}); // column major.

        simulated_right = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.17;
            opt.lower = 1;
            opt.upper = 4;
            opt.seed = 24680;
            return opt;
        }());
        for (auto& x : simulated_right) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_right.reset(new tatami::DenseRowMatrix<double, int, decltype(simulated_right)>(nrow, ncol, simulated_right)); // row major
        sparse_right = tatami::convert_to_compressed_sparse<double, int>(*dense_right, false, {}); // column major.
        return;
    }
};

TEST_F(DelayedBinaryIsometricCompareTest, Equal) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricEqualHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] == simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);

    // For coverage purposes.
    EXPECT_FALSE(op->zero_depends_on_row());
    EXPECT_FALSE(op->non_zero_depends_on_row());
    EXPECT_FALSE(op->zero_depends_on_column());
    EXPECT_FALSE(op->non_zero_depends_on_column());
}

TEST_F(DelayedBinaryIsometricCompareTest, GreaterThan) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricGreaterThanHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] > simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, LessThan) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricLessThanHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] < simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, GreaterThanOrEqual) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricGreaterThanOrEqualHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] >= simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, LessThanOrEqual) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricLessThanOrEqualHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] <= simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, NotEqual) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricNotEqualHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = (simulated_left[i] != simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, NewType) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricEqualHelper<uint8_t, double, int> >();
    tatami::DelayedBinaryIsometricOperation<uint8_t, double, int> dense_umod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse_left, sparse_right, op);

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<uint8_t> urefvec(nrow * ncol);
    for (size_t i = 0; i < urefvec.size(); ++i) {
        urefvec[i] = simulated_left[i] == simulated_right[i];
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

TEST(DelayedBinaryIsometricCompare, BackCompatibility) {
    auto eq = tatami::make_DelayedBinaryIsometricEqual();
    EXPECT_FALSE(eq->is_sparse());
    auto gt = tatami::make_DelayedBinaryIsometricGreaterThan();
    EXPECT_TRUE(gt->is_sparse());
    auto lt = tatami::make_DelayedBinaryIsometricLessThan();
    EXPECT_TRUE(lt->is_sparse());
    auto gte = tatami::make_DelayedBinaryIsometricGreaterThanOrEqual();
    EXPECT_FALSE(gte->is_sparse());
    auto lte = tatami::make_DelayedBinaryIsometricLessThanOrEqual();
    EXPECT_FALSE(lte->is_sparse());
    auto neq = tatami::make_DelayedBinaryIsometricNotEqual();
    EXPECT_TRUE(neq->is_sparse());
}
