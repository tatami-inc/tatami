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

TEST(DelayedBinaryIsometricOp, Misshappen) {
    std::vector<double> src(200);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(10, 20, src));
    auto dense2 = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(20, 10, src));
    tatami_test::throws_error([&]() {
        tatami::make_DelayedBinaryIsometricOp(dense, dense2, tatami::DelayedBinaryBasicMockHelper());
    }, "should be the same");
}

class DelayedBinaryIsometricOpMockTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static int nrow = 23, ncol = 42;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<const tatami::NumericMatrix> dense, sparse, bdense, bsparse, ref;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(simulated)));
        sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); 

        bdense = tatami::make_DelayedBinaryIsometricOp(dense, dense, tatami::DelayedBinaryBasicMockHelper());
        bsparse = tatami::make_DelayedBinaryIsometricOp(sparse, sparse, tatami::DelayedBinaryAdvancedMockHelper());
        ref.reset(new tatami::DenseRowMatrix<double>(nrow, ncol, std::vector<double>(nrow * ncol)));
    }
};

TEST_P(DelayedBinaryIsometricOpMockTest, Basic) {
    EXPECT_FALSE(bdense->sparse());
    EXPECT_EQ(bdense->sparse_proportion(), 0);
    EXPECT_TRUE(bsparse->sparse());
    EXPECT_EQ(bsparse->sparse_proportion(), 1);

    // Spamming a whole stack of tests.
    tatami_test::TestAccessParameters params;
    auto tparam = GetParam();
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);

    tatami_test::test_full_access(params, bdense.get(), ref.get());
    tatami_test::test_block_access(params, bdense.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, bdense.get(), ref.get(), 3, 5);

    tatami_test::test_full_access(params, bsparse.get(), ref.get());
    tatami_test::test_block_access(params, bsparse.get(), ref.get(), 5, 20);
    tatami_test::test_indexed_access(params, bsparse.get(), ref.get(), 2, 4);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBinary,
    DelayedBinaryIsometricOpMockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row access
        ::testing::Values(true, false)  // oracle usage
    )
);
