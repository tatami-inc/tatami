#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"
#include "../_tests/simulate_vector.h"

template<class PARAM>
class TransposeTest: public ::testing::TestWithParam<PARAM> {
protected:
    size_t nrow = 199, ncol = 201;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse, tdense, tsparse, ref;
    std::vector<double> simulated;
protected:
    void SetUp() {
        simulated = simulate_sparse_vector<double>(nrow * ncol, 0.05);
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
    EXPECT_FALSE(tdense->prefer_rows());
    EXPECT_FALSE(tdense->sparse());

    EXPECT_EQ(tsparse->ncol(), nrow);
    EXPECT_EQ(tsparse->nrow(), ncol);
    EXPECT_TRUE(tsparse->prefer_rows());
    EXPECT_TRUE(tsparse->sparse());

    test_simple_row_access(tdense.get(), ref.get(), FORWARD, JUMP);
    test_simple_row_access(tsparse.get(), ref.get(), FORWARD, JUMP);
}

TEST_P(TransposeFullTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    test_simple_column_access(tdense.get(), ref.get(), FORWARD, JUMP);
    test_simple_column_access(tsparse.get(), ref.get(), FORWARD, JUMP);
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

    test_sliced_row_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_row_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(TransposeBlockTest, Column) {
    auto param = GetParam();
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * tdense->nrow(), LAST = interval_info[1] * tdense->nrow();

    test_sliced_column_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_column_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
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

    test_indexed_column_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_column_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(TransposeIndexTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    test_indexed_row_access(tdense.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_row_access(tsparse.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
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
