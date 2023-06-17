#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(CompressedSparseMatrix, ConstructionEmpty) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr);
    EXPECT_TRUE(mat.sparse());
    EXPECT_FALSE(mat.prefer_rows());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    // Comparing access for an empty matrix.
    tatami::DenseColumnMatrix<double, int> dense(10, 20, std::vector<double>(200));
    tatami_test::test_simple_column_access(&mat, &dense, true, 1);
    tatami_test::test_simple_row_access(&mat, &dense, true, 1);

    // Same for row-major.
    indptr.resize(11);
    tatami::CompressedSparseRowMatrix<double, int> rmat(10, 20, values, indices, indptr);
    EXPECT_TRUE(rmat.sparse());
    EXPECT_EQ(rmat.nrow(), 10);
    EXPECT_EQ(rmat.ncol(), 20);
    tatami_test::test_simple_column_access(&rmat, &dense, true, 1);
    tatami_test::test_simple_row_access(&rmat, &dense, true, 1);
}

TEST(CompressedSparseMatrix, ConstructionFail) {
    std::vector<double> values { 0 };
    std::vector<int> indices;
    std::vector<size_t> indptr(21);
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "same length");

    indices.push_back(0);
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 10, values, indices, indptr); }, "should be equal to 'ncol");
    tatami_test::throws_error([&]() { tatami::CompressedSparseRowMatrix<double, int> mat(10, 10, values, indices, indptr); }, "should be equal to 'nrow");

    indptr[0] = 1;
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "should be zero");

    indptr[0] = 0;
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "should be equal to length");

    indptr.back() = 1;
    indptr[1] = 1;
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "non-decreasing");

    indptr[1] = 0;
    indices[0] = -1;
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "non-negative");

    indices[0] = 10001;
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "non-negative");

    indices[0] = 3;
    indices.push_back(2);
    values.push_back(2);
    indptr.back() = indices.size();
    tatami_test::throws_error([&]() { tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr); }, "strictly increasing");
    tatami_test::throws_error([&]() { tatami::CompressedSparseRowMatrix<double, int> mat(20, 10, values, indices, indptr); }, "strictly increasing");
}

/*************************************
 *************************************/

class SparseTestMethods {
protected:
    size_t nrow = 200, ncol = 100;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse_row, sparse_column;

    void assemble() {
        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05)));
        sparse_row = tatami::convert_to_sparse<true>(dense.get());
        sparse_column = tatami::convert_to_sparse<false>(dense.get());
        return;
    }
};

class SparseUtilsTest : public ::testing::Test, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_F(SparseUtilsTest, Basic) {
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

class SparseFullAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t> >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseFullAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    tatami_test::test_simple_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

TEST_P(SparseFullAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

class SparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseSlicedAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    tatami_test::test_sliced_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(SparseSlicedAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    tatami_test::test_sliced_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseSlicedAccessTest,
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

class SparseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseIndexedAccessTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    tatami_test::test_indexed_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(SparseIndexedAccessTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    tatami_test::test_indexed_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseIndexedAccessTest,
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

TEST(CompressedSparseMatrix, SecondarySkip) {
    int nrow = 201, ncol = 12;
    auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05);

    // Injecting some all-zero rows that should be skipped.
    {
        auto start = simulated.data();
        std::fill(start + ncol * 12,  start + ncol * 50,  0);
        std::fill(start + ncol * 70,  start + ncol * 100, 0);
        std::fill(start + ncol * 130, start + ncol * 140, 0);
        std::fill(start + ncol * 160, start + ncol * 190, 0);
    }

    tatami::DenseRowMatrix<double, int> dense(nrow, ncol, std::move(simulated));

    // Confirm that we injected the all-zero rows correctly!
    {
        int all_zero = 0;
        auto wrk = dense.dense_row();
        for (int r = 0; r < nrow; ++r) {
            auto extracted = wrk->fetch(r);
            int non_zero = false;
            for (auto x : extracted) {
                non_zero += x != 0;
            }
            all_zero += (non_zero == 0);
        }
        EXPECT_TRUE(all_zero > 100); 
    }

    // Make a column-sparse compressed matrix, so that we can check
    // that secondary extraction correctly skips the all-zero rows.
    auto sparse_column = tatami::convert_to_sparse<false>(&dense);
    tatami_test::test_simple_row_access(sparse_column.get(), &dense, true, 1); // forward
    tatami_test::test_simple_row_access(sparse_column.get(), &dense, false, 1); // backward
}
