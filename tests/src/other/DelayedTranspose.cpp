#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedTranspose.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

template<class PARAM>
class TransposeTest: public ::testing::TestWithParam<PARAM> {
protected:
    size_t nrow = 199, ncol = 201;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse, tdense, tsparse, ref;
    std::vector<double> simulated;
protected:
    void SetUp() {
        simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
        sparse = tatami::convert_to_compressed_sparse<false>(dense.get()); // column-major.
        tdense = tatami::make_DelayedTranspose(dense);
        tsparse = tatami::make_DelayedTranspose(sparse);

        std::vector<double> refvec(nrow * ncol);
        for (size_t r = 0; r < nrow; ++r) {
            for (size_t c = 0; c < ncol; ++c) {
                refvec[c * nrow + r] = simulated[r * ncol + c];
            }
        }
        ref.reset(new tatami::DenseRowMatrix<double>(ncol, nrow, refvec));
    }
};

using TransposeFullTest = TransposeTest<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t> >;

TEST_P(TransposeFullTest, Row) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    EXPECT_EQ(tdense->ncol(), nrow);
    EXPECT_EQ(tdense->nrow(), ncol);
    EXPECT_FALSE(tdense->sparse());
    EXPECT_EQ(tdense->sparse_proportion(), 0);
    EXPECT_FALSE(tdense->prefer_rows());
    EXPECT_EQ(tdense->prefer_rows_proportion(), 0);

    EXPECT_EQ(tsparse->ncol(), nrow);
    EXPECT_EQ(tsparse->nrow(), ncol);
    EXPECT_TRUE(tsparse->sparse());
    EXPECT_EQ(tsparse->sparse_proportion(), 1);
    EXPECT_TRUE(tsparse->prefer_rows());
    EXPECT_EQ(tsparse->prefer_rows_proportion(), 1);

    tatami_test::test_full_access(params, tdense.get(), ref.get());
    tatami_test::test_full_access(params, tsparse.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    TransposeTest,
    TransposeFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

using TransposeBlockTest = TransposeTest<std::tuple<bool, bool, tatami_test::TestAccessOrder, int, std::pair<double, double> > >;

TEST_P(TransposeBlockTest, Sliced) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (params.use_row ? tdense->ncol() : tdense->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, tdense.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(params, tsparse.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    TransposeTest,
    TransposeBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 4), // jumps (to test workspace memory)
        ::testing::Values(
            std::make_pair(0, 0.44),
            std::make_pair(0.21, 0.89), 
            std::make_pair(0.33, 1)
        )
    )
);

using TransposeIndexTest = TransposeTest<std::tuple<bool, bool, tatami_test::TestAccessOrder, int, std::pair<double, double> > >;

TEST_P(TransposeIndexTest, Column) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters params;
    params.use_row = std::get<0>(tparam);
    params.use_oracle = std::get<1>(tparam);
    params.order = std::get<2>(tparam);
    params.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (params.use_row ? tdense->ncol() : tdense->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(params, tdense.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, tsparse.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    TransposeTest,
    TransposeIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::make_pair(0, 0.05),
            std::make_pair(0.2, 0.1), 
            std::make_pair(0.7, 0.03)
        )
    )
);

TEST(TransposeTest, ConstOverload) {
    int nrow = 10, ncol = 15;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05);
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, simulated));
    auto tdense = tatami::make_DelayedTranspose(dense);

    // Cursory checks.
    EXPECT_EQ(dense->nrow(), tdense->ncol());
    EXPECT_EQ(dense->ncol(), tdense->nrow());
}
