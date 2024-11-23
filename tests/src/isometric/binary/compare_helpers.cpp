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
        dense_left = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated_left));
        sparse_left = tatami::convert_to_compressed_sparse<false, double, int>(dense_left.get()); // column major.

        simulated_right = tatami_test::simulate_vector<double>(nrow * ncol, []{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.17;
            opt.lower = 1;
            opt.upper = 4;
            opt.seed = 24680;
            return opt;
        }());
        for (auto& x : simulated_right) { x = std::round(x); } // Rounding for easier tests of exact equality.
        dense_right = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated_right));
        sparse_right = tatami::convert_to_compressed_sparse<false, double, int>(dense_right.get()); // column major.
        return;
    }
};

TEST_F(DelayedBinaryIsometricCompareTest, Equal) {
    auto op = tatami::make_DelayedBinaryIsometricEqual();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOperation(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOperation(sparse_left, sparse_right, op);

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

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, GreaterThan) {
    auto op = tatami::make_DelayedBinaryIsometricGreaterThan();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOperation(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOperation(sparse_left, sparse_right, op);

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

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, LessThan) {
    auto op = tatami::make_DelayedBinaryIsometricLessThan();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOperation(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOperation(sparse_left, sparse_right, op);

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

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, GreaterThanOrEqual) {
    auto op = tatami::make_DelayedBinaryIsometricGreaterThanOrEqual();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOperation(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOperation(sparse_left, sparse_right, op);

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

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, LessThanOrEqual) {
    auto op = tatami::make_DelayedBinaryIsometricLessThanOrEqual();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOperation(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOperation(sparse_left, sparse_right, op);

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

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, NotEqual) {
    auto op = tatami::make_DelayedBinaryIsometricNotEqual();
    auto dense_mod = tatami::make_DelayedBinaryIsometricOperation(dense_left, dense_right, op);
    auto sparse_mod = tatami::make_DelayedBinaryIsometricOperation(sparse_left, sparse_right, op);

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

    quick_test_all<double, int>(*dense_mod, ref);
    quick_test_all<double, int>(*sparse_mod, ref);
}

TEST_F(DelayedBinaryIsometricCompareTest, NewType) {
    auto op = tatami::make_DelayedBinaryIsometricEqual();
    auto dense_umod = tatami::make_DelayedBinaryIsometricOperation<uint8_t>(dense_left, dense_right, op);
    auto sparse_umod = tatami::make_DelayedBinaryIsometricOperation<uint8_t>(sparse_left, sparse_right, op);

    // Toughest tests are handled by 'arith_helpers.cpp'; they would
    // be kind of redundant here, so we'll just do something simple
    // to check that the operation behaves as expected. 
    std::vector<uint8_t> urefvec(nrow * ncol);
    for (size_t i = 0; i < urefvec.size(); ++i) {
        urefvec[i] = simulated_left[i] == simulated_right[i];
    }
    tatami::DenseRowMatrix<uint8_t, int> uref(nrow, ncol, std::move(urefvec));

    quick_test_all<uint8_t, int>(*dense_umod, uref);
    quick_test_all<uint8_t, int>(*sparse_umod, uref);
}
