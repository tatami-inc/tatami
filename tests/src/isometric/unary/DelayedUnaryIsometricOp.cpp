#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/unary/DelayedUnaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

TEST(DelayedUnaryIsometricOp, ConstOverload) {
    int nrow = 23, ncol = 42;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));

    auto vec = std::vector<double>(nrow);
    auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
    auto mat = tatami::make_DelayedUnaryIsometricOp(dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), dense->nrow());
    EXPECT_EQ(mat->ncol(), dense->ncol());
}

TEST(DelayedUnaryIsometricOp, Mock) {
    int nrow = 23, ncol = 42;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(simulated)));
    auto sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); 

    auto udense = tatami::make_DelayedUnaryIsometricOp(dense, tatami::DelayedUnaryBasicMockHelper());
    EXPECT_FALSE(udense->sparse());
    EXPECT_EQ(udense->sparse_proportion(), 0);

    auto usparse = tatami::make_DelayedUnaryIsometricOp(sparse, tatami::DelayedUnaryAdvancedMockHelper());
    EXPECT_TRUE(usparse->sparse());
    EXPECT_EQ(usparse->sparse_proportion(), 1);

    auto ref = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::vector<double>(nrow * ncol)));

    // Spamming a whole stack of tests.
    tatami_test::TestAccessParameters params;
    tatami_test::test_full_access(params, udense.get(), ref.get());
    tatami_test::test_block_access(params, udense.get(), ref.get(), 10, 30);
    tatami_test::test_indexed_access(params, udense.get(), ref.get(), 5, 5);

    tatami_test::test_full_access(params, usparse.get(), ref.get());
    tatami_test::test_block_access(params, usparse.get(), ref.get(), 10, 30);
    tatami_test::test_indexed_access(params, usparse.get(), ref.get(), 5, 5);

    params.use_row = false; 
    tatami_test::test_full_access(params, udense.get(), ref.get());
    tatami_test::test_block_access(params, udense.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, udense.get(), ref.get(), 2, 3);

    tatami_test::test_full_access(params, usparse.get(), ref.get());
    tatami_test::test_block_access(params, usparse.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, usparse.get(), ref.get(), 2, 3);
}