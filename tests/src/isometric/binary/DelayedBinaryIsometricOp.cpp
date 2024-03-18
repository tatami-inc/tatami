#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/isometric/binary/DelayedBinaryIsometricOp.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"
#include "../utils.h"

TEST(DelayedBinaryIsometricOp, ConstOverload) {
    int nrow = 23, ncol = 42;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(simulated)));

    auto op = tatami::make_DelayedBinaryAddHelper();
    auto mat = tatami::make_DelayedBinaryIsometricOp(dense, dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), nrow);
    EXPECT_EQ(mat->ncol(), ncol);
}

TEST(DelayedBinaryIsometricOp, Mock) {
    int nrow = 23, ncol = 42;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(simulated)));
    auto sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); 

    auto bvdense = tatami::make_DelayedBinaryIsometricOp(dense, dense, tatami::DelayedBinaryMockVariableDenseHelper());
    EXPECT_FALSE(bvdense->sparse());

    auto bcdense = tatami::make_DelayedBinaryIsometricOp(sparse, sparse, tatami::DelayedBinaryMockConstantDenseHelper());
    EXPECT_FALSE(bcdense->sparse());

    auto bsparse = tatami::make_DelayedBinaryIsometricOp(sparse, sparse, tatami::DelayedBinaryMockSparseHelper());
    EXPECT_TRUE(bsparse->sparse());

    auto ref = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::vector<double>(nrow * ncol)));

    // Spamming a whole stack of tests.
    tatami_test::TestAccessParameters params;
    tatami_test::test_full_access(params, bvdense.get(), ref.get());
    tatami_test::test_block_access(params, bvdense.get(), ref.get(), 10, 30);
    tatami_test::test_indexed_access(params, bvdense.get(), ref.get(), 5, 5);

    tatami_test::test_full_access(params, bcdense.get(), ref.get());
    tatami_test::test_block_access(params, bcdense.get(), ref.get(), 10, 30);
    tatami_test::test_indexed_access(params, bcdense.get(), ref.get(), 5, 5);


    tatami_test::test_full_access(params, bsparse.get(), ref.get());
    tatami_test::test_block_access(params, bsparse.get(), ref.get(), 10, 30);
    tatami_test::test_indexed_access(params, bsparse.get(), ref.get(), 5, 5);

    params.use_row = false; 
    tatami_test::test_full_access(params, bvdense.get(), ref.get());
    tatami_test::test_block_access(params, bvdense.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, bvdense.get(), ref.get(), 2, 3);

    tatami_test::test_full_access(params, bcdense.get(), ref.get());
    tatami_test::test_block_access(params, bcdense.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, bcdense.get(), ref.get(), 2, 3);

    tatami_test::test_full_access(params, bsparse.get(), ref.get());
    tatami_test::test_block_access(params, bsparse.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, bsparse.get(), ref.get(), 2, 3);
}
