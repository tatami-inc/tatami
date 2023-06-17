#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/subset/DelayedSubsetBlock.hpp"
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
            dense_block = tatami::make_DelayedSubsetBlock<0>(dense, first, last);
            sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, first, last);
        } else {
            std::vector<double> sub;
            sub.reserve(NR * block_length);
            for (int r = 0; r < NR; ++r) {
                auto row = simulated.data() + r * NC;
                sub.insert(sub.end(), row + first, row + last);
            }
            ref.reset(new tatami::DenseRowMatrix<double>(NR, block_length, std::move(sub)));
            dense_block = tatami::make_DelayedSubsetBlock<1>(dense, first, last);
            sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, first, last);
        }
    }
};

/*****************************
 *****************************/

using SubsetBlockFullAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, bool, int> >;

TEST_P(SubsetBlockFullAccessTest, Row) {
    auto param = GetParam();
    extra_assemble(param);

    EXPECT_EQ(ref->nrow(), dense_block->nrow());
    EXPECT_EQ(ref->ncol(), dense_block->ncol());
    if (std::get<0>(param)) {
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

    bool FORWARD = std::get<2>(param);
    bool JUMP = std::get<3>(param);
    tatami_test::test_simple_row_access(dense_block.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_block.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(SubsetBlockFullAccessTest, Column) {
    auto param = GetParam();
    extra_assemble(param);
    bool FORWARD = std::get<2>(param);
    bool JUMP = std::get<3>(param);

    tatami_test::test_simple_column_access(dense_block.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_block.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubsetBlock,
    SubsetBlockFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(true, false), // iterate forwards or backwards.
        ::testing::Values(1, 3) // jump, to check the workspace memory
    )
);

/*****************************
 *****************************/

using SubsetBlockSlicedAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, bool, int, std::vector<double> > >;

TEST_P(SubsetBlockSlicedAccessTest, Row) {
    auto param = GetParam();
    extra_assemble(param);

    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    int FIRST = interval_info[0] * dense_block->ncol(), LAST = interval_info[1] * dense_block->ncol();

    tatami_test::test_sliced_row_access(dense_block.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(sparse_block.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(SubsetBlockSlicedAccessTest, Column) {
    auto param = GetParam();
    extra_assemble(param);

    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    int FIRST = interval_info[0] * dense_block->nrow(), LAST = interval_info[1] * dense_block->nrow();

    tatami_test::test_sliced_column_access(dense_block.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(sparse_block.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubsetBlock,
    SubsetBlockSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(true, false), // iterate forwards or backwards.
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<double>({ 0, 0.45 }), 
            std::vector<double>({ 0.33, 0.66 }),
            std::vector<double>({ 0.56, 1 })
        )        
    )
);

/*****************************
 *****************************/

using SubsetBlockIndexedAccessTest = SubsetBlockTest<std::tuple<bool, std::pair<double, double>, bool, int, std::vector<double> > >;

TEST_P(SubsetBlockIndexedAccessTest, Row) {
    auto param = GetParam();
    extra_assemble(param);

    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    int FIRST = interval_info[0] * dense_block->ncol(), STEP = interval_info[1] * dense_block->ncol();

    tatami_test::test_indexed_row_access(dense_block.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(sparse_block.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(SubsetBlockIndexedAccessTest, Column) {
    auto param = GetParam();
    extra_assemble(param);

    bool FORWARD = std::get<2>(param);
    int JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    int FIRST = interval_info[0] * dense_block->nrow(), STEP = interval_info[1] * dense_block->nrow();

    tatami_test::test_indexed_column_access(dense_block.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(sparse_block.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubsetBlock,
    SubsetBlockIndexedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(true, false), // iterate forwards or backwards.
        ::testing::Values(1, 3), // jump, to check the workspace memory
        ::testing::Values(
            std::vector<double>({ 0, 0.05 }), 
            std::vector<double>({ 0.33, 0.06 }),
            std::vector<double>({ 0.56, 0.02 })
        )        
    )
);

/*****************************
 *****************************/

class SubsetBlockOracleTest : public ::testing::TestWithParam<std::tuple<bool, std::pair<double, double>, bool> > {
protected:
    int NR = 155, NC = 102;
    std::shared_ptr<tatami::NumericMatrix> dense_block, sparse_block, wrapped_dense_block, wrapped_sparse_block;

protected:
    template<class PARAM>
    void assemble(const PARAM& param) {
        auto simulated = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.2);
        auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, simulated));
        auto sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.

        double full =  (std::get<0>(param) ? NR : NC);
        int first = full * std::get<1>(param).first;
        int last = full * std::get<1>(param).second;
        auto block_length = last - first;

        if (std::get<0>(param)) {
            dense_block = tatami::make_DelayedSubsetBlock<0>(dense, first, last);
            sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, first, last);
            wrapped_dense_block = tatami::make_DelayedSubsetBlock<0>(tatami_test::make_CrankyMatrix(dense), first, last);
            wrapped_sparse_block = tatami::make_DelayedSubsetBlock<0>(tatami_test::make_CrankyMatrix(sparse), first, last);
        } else {
            dense_block = tatami::make_DelayedSubsetBlock<1>(dense, first, last);
            sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, first, last);
            wrapped_dense_block = tatami::make_DelayedSubsetBlock<1>(tatami_test::make_CrankyMatrix(dense), first, last);
            wrapped_sparse_block = tatami::make_DelayedSubsetBlock<1>(tatami_test::make_CrankyMatrix(sparse), first, last);
        }
    }
};

TEST_P(SubsetBlockOracleTest, Validate) {
    auto param = GetParam();
    assemble(param);
    auto random = std::get<2>(param);

    EXPECT_FALSE(dense_block->uses_oracle(true));
    EXPECT_TRUE(wrapped_dense_block->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_dense_block.get(), dense_block.get(), random);
    tatami_test::test_oracle_column_access(wrapped_sparse_block.get(), sparse_block.get(), random);

    tatami_test::test_oracle_row_access(wrapped_dense_block.get(), dense_block.get(), random);
    tatami_test::test_oracle_row_access(wrapped_sparse_block.get(), sparse_block.get(), random);
}

INSTANTIATE_TEST_CASE_P(
    DelayedSubsetBlock,
    SubsetBlockOracleTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row or column subsetting, respectively.
        ::testing::Values(
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.8),
            std::make_pair(0.4, 1)
        ),
        ::testing::Values(true, false)  // use random or consecutive oracle.
    )
);

/****************************************************
 ****************************************************/

TEST(DelayedSubsetBlock, ConstOverload) {
    int NR = 9, NC = 7;
    auto dense = std::shared_ptr<const tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
    auto sub = tatami::make_DelayedSubsetBlock<1>(dense, static_cast<int>(5), static_cast<int>(8));
    EXPECT_EQ(sub->ncol(), 3);
    EXPECT_EQ(sub->nrow(), NR);
}

