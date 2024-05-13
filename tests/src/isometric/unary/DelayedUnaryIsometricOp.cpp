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
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulated));

    auto vec = std::vector<double>(nrow);
    auto op = tatami::make_DelayedAddVectorHelper<0>(vec);
    auto mat = tatami::make_DelayedUnaryIsometricOp(dense, std::move(op));

    // cursory checks.
    EXPECT_EQ(mat->nrow(), dense->nrow());
    EXPECT_EQ(mat->ncol(), dense->ncol());
}

class DelayedUnaryIsometricOpMockTest : public ::testing::TestWithParam<std::tuple<bool, bool> > {
protected:
    inline static int nrow = 57, ncol = 37;
    inline static std::vector<double> simulated;
    inline static std::shared_ptr<const tatami::NumericMatrix> dense, sparse, udense, usparse, ref;

    static void SetUpTestSuite() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.1);
        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));
        sparse = tatami::convert_to_compressed_sparse<false, double, int>(dense.get()); 

        udense = tatami::make_DelayedUnaryIsometricOp(dense, tatami::DelayedUnaryBasicMockHelper());
        usparse = tatami::make_DelayedUnaryIsometricOp(sparse, tatami::DelayedUnaryAdvancedMockHelper());
        ref.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::vector<double>(nrow * ncol)));
    }
};

TEST_P(DelayedUnaryIsometricOpMockTest, Mock) {
    EXPECT_FALSE(udense->sparse());
    EXPECT_EQ(udense->sparse_proportion(), 0);

    EXPECT_TRUE(usparse->sparse());
    EXPECT_EQ(usparse->sparse_proportion(), 1);

    // Spamming a whole stack of tests.
    tatami_test::TestAccessParameters params;
    auto tparam = GetParam();
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);

    tatami_test::test_full_access(params, udense.get(), ref.get());
    tatami_test::test_block_access(params, udense.get(), ref.get(), 5, 30);
    tatami_test::test_indexed_access(params, udense.get(), ref.get(), 3, 5);

    tatami_test::test_full_access(params, usparse.get(), ref.get());
    tatami_test::test_block_access(params, usparse.get(), ref.get(), 5, 30);
    tatami_test::test_indexed_access(params, usparse.get(), ref.get(), 2, 4);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedUnary,
    DelayedUnaryIsometricOpMockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row access
        ::testing::Values(true, false)  // oracle usage
    )
);
