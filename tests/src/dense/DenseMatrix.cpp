#include <gtest/gtest.h>

#include <vector>
#include <deque>
#include <numeric>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami_test/tatami_test.hpp"

TEST(DenseMatrix, Basic) {
    std::vector<double> contents(200);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double> mat(10, 20, contents);

    EXPECT_FALSE(mat.sparse());
    EXPECT_EQ(mat.sparse_proportion(), 0);
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.prefer_rows_proportion(), 0);

    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    {
        auto wrk = mat.dense_column();
        for (size_t i = 0; i < mat.ncol(); ++i) {
            auto start = contents.begin() + i * mat.nrow();
            std::vector<double> expected(start, start + mat.nrow());
            EXPECT_EQ(wrk->fetch(i), expected);
        }
    }

    {
        auto wrk = mat.dense_row();
        for (size_t i = 0; i < mat.nrow(); ++i) {
            std::vector<double> expected(mat.ncol());
            for (size_t j = 0; j < mat.ncol(); ++j) {
                expected[j] = contents[j * mat.nrow() + i];
            }
            EXPECT_EQ(wrk->fetch(i), expected);
        }
    }

    EXPECT_FALSE(mat.uses_oracle(true));
    {
        auto wrk = mat.dense_row();
        wrk->set_oracle(nullptr); // no-op.
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
        auto simulated = tatami_test::simulate_dense_vector<double>(nrow * ncol, 0.05);

        auto transposed = std::vector<double>(nrow * ncol);
        for (size_t c = 0; c < ncol; ++c) {
            for (size_t r = 0; r < nrow; ++r) {
                transposed[c * nrow + r] = simulated[r * ncol + c];
            }
        }

        dense_row.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));
        dense_column.reset(new tatami::DenseColumnMatrix<double, int>(nrow, ncol, std::move(transposed)));
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

TEST_P(DenseFullAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_column_access(dense_row.get(), dense_column.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(dense_column.get(), dense_row.get(), FORWARD, JUMP);
}

TEST_P(DenseFullAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_row_access(dense_row.get(), dense_column.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(dense_column.get(), dense_row.get(), FORWARD, JUMP);
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

class DenseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(DenseSlicedAccessTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    tatami_test::test_sliced_column_access(dense_column.get(), dense_row.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(dense_row.get(), dense_column.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(DenseSlicedAccessTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    tatami_test::test_sliced_row_access(dense_column.get(), dense_row.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(dense_row.get(), dense_column.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    DenseMatrix,
    DenseSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.45 }),
            std::vector<double>({ 0.2, 0.8 }), 
            std::vector<double>({ 0.7, 1 })
        )
    )
);

/*************************************
 *************************************/

class DenseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public DenseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(DenseIndexedAccessTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    tatami_test::test_indexed_column_access(dense_column.get(), dense_row.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(dense_row.get(), dense_column.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(DenseIndexedAccessTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    tatami_test::test_indexed_row_access(dense_column.get(), dense_row.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(dense_row.get(), dense_column.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    DenseMatrix,
    DenseIndexedAccessTest,
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

/*************************************
 *************************************/

TEST(DenseMatrix, TypeOverflow) {
    // Check for correct pointer arithmetic when indices are supplied in a smaller integer type;
    // the class should convert them into size_t before doing any work.
    std::vector<double> contents(20000);
    double counter = -105;
    for (auto& i : contents) { i = counter++; }

    tatami::DenseColumnMatrix<double> ref(100, 200, contents);
    tatami::DenseColumnMatrix<double, unsigned char> limited(100, 200, contents);

    EXPECT_EQ(limited.nrow(), 100);
    EXPECT_EQ(limited.ncol(), 200);

    {
        auto rwrk = ref.dense_column();
        auto lwrk = limited.dense_column();
        for (int i = 0; i < ref.ncol(); ++i) {
            EXPECT_EQ(rwrk->fetch(i), lwrk->fetch(i));
        }
    }

    {
        auto rwrk = ref.dense_row();
        auto lwrk = limited.dense_row();
        for (int i = 0; i < ref.nrow(); ++i) {
            EXPECT_EQ(rwrk->fetch(i), lwrk->fetch(i));
        }
    }

    {
        auto rwrk = ref.dense_column(59, 189);
        auto lwrk = limited.dense_column(59, 189);
        for (int i = 0; i < ref.ncol(); ++i) {
            EXPECT_EQ(rwrk->fetch(i), lwrk->fetch(i));
        }
    }

    {
        auto rwrk = ref.dense_row(59, 89);
        auto lwrk = limited.dense_row(59, 89);
        for (int i = 0; i < ref.nrow(); ++i) {
            EXPECT_EQ(rwrk->fetch(i), lwrk->fetch(i));
        }
    }

    {
        std::vector<int> indices{ 11, 33, 55, 99, 111, 122, 155, 177, 199 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());

        auto rwrk = ref.dense_column(std::move(indices));
        auto lwrk = limited.dense_column(std::move(uindices));
        for (int i = 0; i < ref.ncol(); ++i) {
            EXPECT_EQ(rwrk->fetch(i), lwrk->fetch(i));
        }
    }

    {
        std::vector<int> indices{ 10, 20, 40, 60, 80, 99 };
        std::vector<unsigned char> uindices(indices.begin(), indices.end());

        auto rwrk = ref.dense_row(std::move(indices));
        auto lwrk = limited.dense_row(std::move(uindices));
        for (int i = 0; i < ref.nrow(); ++i) {
            EXPECT_EQ(rwrk->fetch(i), lwrk->fetch(i));
        }
    }
}
