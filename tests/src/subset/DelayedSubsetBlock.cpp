#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/subset/DelayedSubsetBlock.hpp"
//#include "tatami/subset/make_DelayedSubset.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

template<class PARAM> 
class SubsetBlockTest : public ::testing::TestWithParam<PARAM> {
protected:
    int NR = 192, NC = 132;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse, ref, dense_block, sparse_block;
    std::vector<double> simulated;

    int block_length;
protected:
    void SetUp() {
        simulated = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.2);
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, simulated));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.
        return;
    }

    void extra_assemble(const PARAM& param) {
        double full =  (std::get<0>(param) ? NR : NC);
        int first = full * std::get<1>(param).first;
        int last = full * std::get<1>(param).second;
        block_length = last - first;

        if (std::get<0>(param)) {
            std::vector<double> sub(simulated.data() + first * NC, simulated.data() + last * NC);
            ref.reset(new tatami::DenseRowMatrix<double>(block_length, NC, std::move(sub)));
            dense_block = tatami::make_DelayedSubsetBlock<0>(dense, first, block_length);
            sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, first, block_length);
        } else {
            std::vector<double> sub;
            sub.reserve(NR * block_length);
            for (int r = 0; r < NR; ++r) {
                auto row = simulated.data() + r * NC;
                sub.insert(sub.end(), row + first, row + last);
            }
            ref.reset(new tatami::DenseRowMatrix<double>(NR, block_length, std::move(sub)));
            dense_block = tatami::make_DelayedSubsetBlock<1>(dense, first, block_length);
            sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, first, block_length);
        }
    }
};

/*****************************
 *****************************/

using SubsetBlockFullAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, bool, bool, tatami_test::TestAccessOrder, int> >;

TEST_P(SubsetBlockFullAccessTest, Full) {
    auto tparam = GetParam();
    extra_assemble(tparam);

    EXPECT_EQ(ref->nrow(), dense_block->nrow());
    EXPECT_EQ(ref->ncol(), dense_block->ncol());
    if (std::get<0>(tparam)) {
        EXPECT_EQ(block_length, dense_block->nrow());
        EXPECT_EQ(dense->ncol(), dense_block->ncol());
    } else {
        EXPECT_EQ(dense->nrow(), dense_block->nrow());
        EXPECT_EQ(block_length, dense_block->ncol());
    }

    EXPECT_FALSE(dense_block->sparse());
    EXPECT_EQ(dense_block->sparse_proportion(), 0);
    EXPECT_TRUE(sparse_block->sparse());
    EXPECT_EQ(sparse_block->sparse_proportion(), 1);

    EXPECT_TRUE(dense_block->prefer_rows());
    EXPECT_EQ(dense_block->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_block->prefer_rows());
    EXPECT_EQ(sparse_block->prefer_rows_proportion(), 0);

    tatami_test::TestAccessParameters param;
    param.use_row = std::get<2>(tparam);
    param.use_oracle = std::get<3>(tparam);
    param.order = std::get<4>(tparam);
    param.jump = std::get<5>(tparam);

    tatami_test::test_full_access(param, dense_block.get(), ref.get());
    tatami_test::test_full_access(param, sparse_block.get(), ref.get());
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values( // the block dimensions.
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(true, false), // row or column access.
        ::testing::Values(true, false), // whether to use an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3) // jump, to check the workspace memory
    )
);

/*****************************
 *****************************/

using SubsetBlockSlicedAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, bool, bool, tatami_test::TestAccessOrder, int, std::pair<double, double> > >;

TEST_P(SubsetBlockSlicedAccessTest, Sliced) {
    auto tparam = GetParam();
    extra_assemble(tparam);

    tatami_test::TestAccessParameters param;
    param.use_row = std::get<2>(tparam);
    param.use_oracle = std::get<3>(tparam);
    param.order = std::get<4>(tparam);
    param.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (param.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(param, dense_block.get(), ref.get(), FIRST, LAST);
    tatami_test::test_block_access(param, sparse_block.get(), ref.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(true, false), // row or column access.
        ::testing::Values(true, false), // whether to use an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::make_pair(0.0, 0.45), 
            std::make_pair(0.33, 0.66),
            std::make_pair(0.56, 1.0)
        )
    )
);

/*****************************
 *****************************/

using SubsetBlockIndexedAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, bool, bool, tatami_test::TestAccessOrder, int, std::pair<double, double> > >;

TEST_P(SubsetBlockIndexedAccessTest, Indexed) {
    auto tparam = GetParam();
    extra_assemble(tparam);

    tatami_test::TestAccessParameters param;
    param.use_row = std::get<2>(tparam);
    param.use_oracle = std::get<3>(tparam);
    param.order = std::get<4>(tparam);
    param.jump = std::get<5>(tparam);

    auto interval_info = std::get<6>(tparam);
    auto len = (param.use_row ? ref->ncol() : ref->nrow());
    size_t FIRST = interval_info.first * len, STEP = interval_info.second * len;

    tatami_test::test_indexed_access(param, dense_block.get(), ref.get(), FIRST, STEP);
    tatami_test::test_indexed_access(param, sparse_block.get(), ref.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    DelayedSubsetBlock,
    SubsetBlockIndexedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1.0)
        ),
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::make_pair(0.0, 0.05), 
            std::make_pair(0.33, 0.06),
            std::make_pair(0.56, 0.02)
        )
    )
);

/****************************************************
 ****************************************************/

TEST(DelayedSubsetBlock, ConstOverload) {
    int NR = 9, NC = 7;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
    auto sub = tatami::make_DelayedSubsetBlock<1>(dense, static_cast<int>(5), static_cast<int>(3));
    EXPECT_EQ(sub->ncol(), 3);
    EXPECT_EQ(sub->nrow(), NR);
}

//TEST(DelayedSubsetBlock, CorrectMaker) {
//    int NR = 90, NC = 50;
//    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_dense_vector<double>(NR * NC)));
//
//    // Checking that the make function dispatches correctly to the block subset class.
//    {
//        std::vector<int> indices { 5, 6, 7, 8, 9, 10 };
//        auto sub = tatami::make_DelayedSubset<1>(dense, indices);
//        EXPECT_EQ(sub->ncol(), 6);
//        EXPECT_EQ(sub->nrow(), NR);
//
//        auto ref = tatami::make_DelayedSubsetBlock<1>(dense, static_cast<int>(5), static_cast<int>(6));
//        tatami_test::test_simple_row_access(sub.get(), ref.get(), true, 1);
//        tatami_test::test_simple_column_access(sub.get(), ref.get(), true, 1);
//    }
//
//    // Checking that it behaves correctly with an empty index vector.
//    {
//        std::vector<int> indices;
//        auto sub = tatami::make_DelayedSubset<0>(dense, indices);
//        EXPECT_EQ(sub->ncol(), NC);
//        EXPECT_EQ(sub->nrow(), 0);
//    }
//}
