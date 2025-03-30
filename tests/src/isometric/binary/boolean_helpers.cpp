#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <tuple>
#include <cstdint>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOperation.hpp"
#include "tatami/isometric/binary/boolean_helpers.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

class DelayedBinaryIsometricBooleanTest : public ::testing::Test {
protected:
    inline static size_t nrow = 123, ncol = 155;
    inline static std::shared_ptr<tatami::NumericMatrix> dense_left, sparse_left, dense_right, sparse_right;
    inline static std::vector<double> simulated_left, simulated_right;

    static void SetUpTestSuite() {
        simulated_left = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = 9876;
            return opt;
        }());
        dense_left.reset(new tatami::DenseMatrix<double, int, decltype(simulated_left)>(nrow, ncol, simulated_left, true)); // row major.
        sparse_left = tatami::convert_to_compressed_sparse<double, int>(*dense_left, false, {}); // column major.

        simulated_right = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.25;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = 54321;
            return opt;
        }());
        dense_right.reset(new tatami::DenseMatrix<double, int, decltype(simulated_right)>(nrow, ncol, simulated_right, true)); // row major.
        sparse_right = tatami::convert_to_compressed_sparse<double, int>(*dense_right, false, {}); // column major.
        return;
    }
};

TEST_F(DelayedBinaryIsometricBooleanTest, EQUAL) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricBooleanEqualHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = static_cast<bool>(simulated_left[i]) == static_cast<bool>(simulated_right[i]);
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

TEST_F(DelayedBinaryIsometricBooleanTest, AND) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricBooleanAndHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = static_cast<bool>(simulated_left[i]) && static_cast<bool>(simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricBooleanTest, OR) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricBooleanOrHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = static_cast<bool>(simulated_left[i]) || static_cast<bool>(simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricBooleanTest, XOR) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricBooleanXorHelper<double, double, int> >();
    tatami::DelayedBinaryIsometricOperation<double, double, int> dense_mod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<double, double, int> sparse_mod(sparse_left, sparse_right, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<double> refvec(nrow * ncol);
    for (size_t i = 0; i < refvec.size(); ++i) {
        refvec[i] = static_cast<bool>(simulated_left[i]) != static_cast<bool>(simulated_right[i]);
    }
    tatami::DenseMatrix<double, int, decltype(refvec)> ref(nrow, ncol, std::move(refvec), true);

    quick_test_all<double, int>(dense_mod, ref);
    quick_test_all<double, int>(sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricBooleanTest, NewType) {
    auto op = std::make_shared<tatami::DelayedBinaryIsometricBooleanEqualHelper<uint8_t, double, int> >();
    tatami::DelayedBinaryIsometricOperation<uint8_t, double, int> dense_umod(dense_left, dense_right, op);
    tatami::DelayedBinaryIsometricOperation<uint8_t, double, int> sparse_umod(sparse_left, sparse_right, op);

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<uint8_t> urefvec(nrow * ncol);
    for (size_t i = 0; i < urefvec.size(); ++i) {
        urefvec[i] = static_cast<bool>(simulated_left[i]) == static_cast<bool>(simulated_right[i]);
    }
    tatami::DenseMatrix<uint8_t, int, decltype(urefvec)> uref(nrow, ncol, std::move(urefvec), true);

    quick_test_all<uint8_t, int>(dense_umod, uref);
    quick_test_all<uint8_t, int>(sparse_umod, uref);
}

TEST(DelayedBinaryIsometricBoolean, BackCompatibility) {
    auto eq = tatami::make_DelayedBinaryIsometricBooleanEqual();
    EXPECT_FALSE(eq->is_sparse());
    auto andop = tatami::make_DelayedBinaryIsometricBooleanAnd();
    EXPECT_TRUE(andop->is_sparse());
    auto orop = tatami::make_DelayedBinaryIsometricBooleanOr();
    EXPECT_TRUE(orop->is_sparse());
    auto xorop = tatami::make_DelayedBinaryIsometricBooleanXor();
    EXPECT_TRUE(xorop->is_sparse());
}
