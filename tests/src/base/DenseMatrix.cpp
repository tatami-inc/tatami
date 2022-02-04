#include <gtest/gtest.h>

#include <vector>
#include <deque>
#include <numeric>

#include "tatami/base/DenseMatrix.hpp"
#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/simulate_vector.h"

TEST(DenseMatrix, Basic) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double> mat(10, 20, contents);
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    for (size_t i = 0; i < mat.ncol(); ++i) {
        auto start = contents.begin() + i * mat.nrow();
        std::vector<double> expected(start, start + mat.nrow());
        EXPECT_EQ(mat.column(i), expected);
    }

    for (size_t i = 0; i < mat.nrow(); ++i) {
        std::vector<double> expected(mat.ncol());
        for (size_t j = 0; j < mat.ncol(); ++j) {
            expected[j] = contents[j * mat.nrow() + i];
        }
        EXPECT_EQ(mat.row(i), expected);
    }
}

TEST(DenseMatrix, OddsAndEnds) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    // Checks run properly.
    contents.clear();
    EXPECT_ANY_THROW({
        tatami::DenseColumnMatrix<double> mat(10, 20, contents);
    });
    EXPECT_ANY_THROW({
        tatami::DenseColumnMatrix<double> mat(10, 20, std::move(contents));
    });

    std::deque<double> more_contents(200);
    std::iota(more_contents.begin(), more_contents.end(), 1);
    tatami::DenseColumnMatrix<double, int, std::deque<double> > mat2(10, 20, more_contents);
    EXPECT_EQ(more_contents.size(), 200);
}

/*************************************
 *************************************/

class DenseTestMethods {
protected:
    size_t nrow = 200, ncol = 100;
    std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column;

    void assemble() {
        dense_row.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulate_dense_vector<double>(nrow * ncol, 0.05)));

        std::vector<double> buffer(nrow * ncol);
        for (size_t i = 0; i < ncol; ++i) {
            dense_row->column_copy(i, buffer.data() + i * nrow);
        }
        dense_column.reset(new tatami::DenseColumnMatrix<double, int>(nrow, ncol, std::move(buffer)));

        return;
    }
};

class DenseUtilsTest : public ::testing::Test, public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_F(DenseUtilsTest, Basic) {
    size_t NC = dense_column->ncol(), NR = dense_column->nrow();
    EXPECT_EQ(NC, ncol);
    EXPECT_EQ(NR, nrow);

    EXPECT_FALSE(dense_column->prefer_rows());
    EXPECT_TRUE(dense_row->prefer_rows());

    EXPECT_FALSE(dense_column->sparse());
    EXPECT_FALSE(dense_row->sparse());
}

/*************************************
 *************************************/

class DenseFullAccessTest : public ::testing::TestWithParam<std::tuple<bool, int> >, public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(DenseFullAccessTest, Basic) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    test_simple_column_access(dense_row.get(), dense_column.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    DenseMatrix,
    DenseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

class DenseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<size_t> > >, public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(DenseSlicedAccessTest, Basic) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    test_sliced_column_access(dense_column.get(), dense_row.get(), FORWARD, JUMP, FIRST, LEN, SHIFT);
}

INSTANTIATE_TEST_CASE_P(
    DenseMatrix,
    DenseSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);
