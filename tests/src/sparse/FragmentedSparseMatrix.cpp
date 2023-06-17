#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/FragmentedSparseMatrix.hpp"
#include "tatami/utils/ArrayView.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(FragmentedSparseMatrix, ConstructionEmpty) {
    std::vector<std::vector<double> > values(20);
    std::vector<std::vector<int> > indices(20);

    tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 20, values, indices);
    EXPECT_TRUE(mat.sparse());
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    // Comparing access for an empty matrix.
    tatami::DenseColumnMatrix<double, int> dense(10, 20, std::vector<double>(200));
    tatami_test::test_simple_column_access(&mat, &dense, true, 1);
    tatami_test::test_simple_row_access(&mat, &dense, true, 1);

    // Same for row-major.
    values.resize(10);
    indices.resize(10);
    tatami::FragmentedSparseRowMatrix<double, int> rmat(10, 20, values, indices);
    EXPECT_TRUE(rmat.sparse());
    EXPECT_EQ(rmat.nrow(), 10);
    EXPECT_EQ(rmat.ncol(), 20);
    tatami_test::test_simple_column_access(&rmat, &dense, true, 1);
    tatami_test::test_simple_row_access(&rmat, &dense, true, 1);
}

TEST(FragmentedSparseMatrix, ConstructionFail) {
    std::vector<std::vector<double> > values(20);
    std::vector<std::vector<int> > indices(19);
    tatami_test::throws_error([&]() { tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 20, values, indices); }, "same length");

    indices.resize(20);
    tatami_test::throws_error([&]() { tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 10, values, indices); }, "should be equal to number of columns");
    tatami_test::throws_error([&]() { tatami::FragmentedSparseRowMatrix<double, int> mat(10, 10, values, indices); }, "should be equal to number of rows");

    indices.front().resize(10);
    tatami_test::throws_error([&]() { tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 20, values, indices); }, "same length");

    values.front().resize(10);
    tatami_test::throws_error([&]() { tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 20, values, indices); }, "strictly increasing");
    tatami_test::throws_error([&]() { tatami::FragmentedSparseRowMatrix<double, int> mat(20, 10, values, indices); }, "strictly increasing");

    values.front().resize(1);
    indices.front().resize(1);
    indices.front()[0] = -1;
    tatami_test::throws_error([&]() { tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 20, values, indices); }, "non-negative");

    indices.front()[0] = 10001;
    tatami_test::throws_error([&]() { tatami::FragmentedSparseColumnMatrix<double, int> mat(10, 20, values, indices); }, "non-negative");
}

/*************************************
 *************************************/

class FragmentedSparseTestMethods {
protected:
    size_t nrow = 200, ncol = 100;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse_row, sparse_column;

    void assemble() {
        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05)));

        {
            std::vector<std::vector<double> > values(nrow);
            std::vector<std::vector<int> > indices(nrow);
            std::vector<double> buffer(ncol);

            auto wrk = dense->dense_row();
            for (size_t r = 0; r < nrow; ++r) {
                auto content = wrk->fetch(r, buffer.data());
                for (size_t c = 0; c < ncol; ++c) {
                    if (content[c]) {
                        values[r].push_back(content[c]);
                        indices[r].push_back(c);
                    }
                }
            }

            sparse_row.reset(new tatami::FragmentedSparseRowMatrix<double, int>(nrow, ncol, std::move(values), std::move(indices)));
        }

        {
            std::vector<std::vector<double> > values(ncol);
            std::vector<std::vector<int> > indices(ncol);
            std::vector<double> buffer(nrow);

            auto wrk = dense->dense_column();
            for (size_t c = 0; c < ncol; ++c) {
                auto content = wrk->fetch(c, buffer.data());
                for (size_t r = 0; r < nrow; ++r) {
                    if (content[r]) {
                        values[c].push_back(content[r]);
                        indices[c].push_back(r);
                    }
                }
            }

            sparse_column.reset(new tatami::FragmentedSparseColumnMatrix<double, int>(nrow, ncol, std::move(values), std::move(indices)));
        }

        return;
    }
};

class FragmentedSparseUtilsTest : public ::testing::Test, public FragmentedSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_F(FragmentedSparseUtilsTest, Basic) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();
    EXPECT_EQ(NC, ncol);
    EXPECT_EQ(NR, nrow);
    EXPECT_EQ(sparse_row->ncol(), ncol);
    EXPECT_EQ(sparse_row->nrow(), nrow);

    EXPECT_FALSE(dense->sparse());
    EXPECT_TRUE(sparse_row->sparse());
    EXPECT_EQ(sparse_row->sparse_proportion(), 1);
    EXPECT_TRUE(sparse_column->sparse());
    EXPECT_EQ(sparse_column->sparse_proportion(), 1);

    EXPECT_TRUE(sparse_row->prefer_rows());
    EXPECT_EQ(sparse_row->prefer_rows_proportion(), 1);
    EXPECT_FALSE(sparse_column->prefer_rows());
    EXPECT_EQ(sparse_column->prefer_rows_proportion(), 0);

    EXPECT_FALSE(sparse_row->uses_oracle(true));
    {
        auto wrk = sparse_row->sparse_column();
        wrk->set_oracle(nullptr); // no-op.
    }
}

/*************************************
 *************************************/

class FragmentedSparseFullAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t> >, public FragmentedSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(FragmentedSparseFullAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

TEST_P(FragmentedSparseFullAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    FragmentedSparseMatrix,
    FragmentedSparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

class FragmentedSparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public FragmentedSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(FragmentedSparseSlicedAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    tatami_test::test_sliced_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(FragmentedSparseSlicedAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    tatami_test::test_sliced_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    FragmentedSparseMatrix,
    FragmentedSparseSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.51 }),
            std::vector<double>({ 0.25, 0.9 }), 
            std::vector<double>({ 0.63, 1 })
        )
    )
);

/*************************************
 *************************************/

class FragmentedSparseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public FragmentedSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(FragmentedSparseIndexedAccessTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    tatami_test::test_indexed_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(FragmentedSparseIndexedAccessTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    tatami_test::test_indexed_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    FragmentedSparseMatrix,
    FragmentedSparseIndexedAccessTest,
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

TEST(FragmentedSparseMatrix, ArrayView) {
    std::vector<tatami::ArrayView<double> > values(20, tatami::ArrayView<double>(NULL, 0));
    std::vector<tatami::ArrayView<int> > indices(20, tatami::ArrayView<int>(NULL, 0));

    tatami::FragmentedSparseColumnMatrix<double, int, decltype(values), decltype(indices)> mat(10, 20, values, indices);
    EXPECT_TRUE(mat.sparse());
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
}
