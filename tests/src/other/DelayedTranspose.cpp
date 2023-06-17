#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/other/DelayedTranspose.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

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
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.
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

using TransposeFullTest = TransposeTest<std::tuple<bool, size_t> >;

TEST_P(TransposeFullTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

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

    tatami_test::test_simple_row_access(tdense.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(tsparse.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(TransposeFullTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_column_access(tdense.get(), ref.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(tsparse.get(), ref.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    TransposeTest,
    TransposeFullTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

using TransposeBlockTest = TransposeTest<std::tuple<bool, int, std::vector<double> > >;

TEST_P(TransposeBlockTest, Row) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * tdense->ncol(), LAST = interval_info[1] * tdense->ncol();

    tatami_test::test_sliced_row_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(TransposeBlockTest, Column) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * tdense->nrow(), LAST = interval_info[1] * tdense->nrow();

    tatami_test::test_sliced_column_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    TransposeTest,
    TransposeBlockTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4), // jumps (to test workspace memory)
        ::testing::Values(
            std::vector<double>({ 0, 0.44 }),
            std::vector<double>({ 0.21, 0.89 }), 
            std::vector<double>({ 0.33, 1 })
        )
    )
);

using TransposeIndexTest = TransposeTest<std::tuple<bool, int, std::vector<double> > >;

TEST_P(TransposeIndexTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    tatami_test::test_indexed_column_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(TransposeIndexTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    tatami_test::test_indexed_row_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    TransposeTest,
    TransposeIndexTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.05 }),
            std::vector<double>({ 0.2, 0.1 }), 
            std::vector<double>({ 0.7, 0.03 })
        )
    )
);

class TransposeOracleTest : public ::testing::TestWithParam<bool> {
protected:
    size_t nrow = 199, ncol = 201;
    std::shared_ptr<tatami::NumericMatrix> tdense, tsparse, wrapped_dense, wrapped_sparse;

    void extra_assemble() {
        auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05);
        auto dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(nrow, ncol, std::move(simulated)));
        auto sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.
        tdense = tatami::make_DelayedTranspose(dense);
        tsparse = tatami::make_DelayedTranspose(sparse);
        wrapped_dense = tatami::make_DelayedTranspose(tatami_test::make_CrankyMatrix(dense));
        wrapped_sparse = tatami::make_DelayedTranspose(tatami_test::make_CrankyMatrix(sparse));
    }
};

TEST_P(TransposeOracleTest, Validate) {
    auto random = GetParam();
    extra_assemble();
    EXPECT_FALSE(tdense->uses_oracle(true));
    EXPECT_TRUE(wrapped_dense->uses_oracle(true));

    tatami_test::test_oracle_column_access(wrapped_dense.get(), tdense.get(), random);
    tatami_test::test_oracle_column_access(wrapped_sparse.get(), tsparse.get(), random);

    tatami_test::test_oracle_row_access(wrapped_dense.get(), tdense.get(), random);
    tatami_test::test_oracle_row_access(wrapped_sparse.get(), tsparse.get(), random);
}

INSTANTIATE_TEST_CASE_P(
    TransposeTest,
    TransposeOracleTest,
    ::testing::Values(true, false)  // use random or consecutive oracle.
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
