#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/SemiCompressedSparseMatrix.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(SemiCompressedSparseMatrix, ConstructionEmpty) {
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    tatami::SemiCompressedSparseColumnMatrix<double, int> mat(10, 20, indices, indptr);
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
    tatami::SemiCompressedSparseRowMatrix<double, int> rmat(10, 20, indices, indptr);
    EXPECT_TRUE(rmat.sparse());
    EXPECT_EQ(rmat.nrow(), 10);
    EXPECT_EQ(rmat.ncol(), 20);
    tatami_test::test_simple_column_access(&rmat, &dense, true, 1);
    tatami_test::test_simple_row_access(&rmat, &dense, true, 1);
}

/*************************************
 *************************************/

class SemiSparseTestMethods {
protected:
    size_t nrow = 243, ncol = 199;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse_row, sparse_column;

    void assemble() {
        auto simulated = tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05, 1, 4);
        for (auto& s : simulated) {
            s = static_cast<int>(s);
        }
        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, std::move(simulated)));

        // Assembling my lovely indexed matrix.
        {
            std::vector<int> expanded_by_col;
            std::vector<size_t> hits(ncol + 1);
            auto dcwork = dense->dense_column();
            for (size_t c = 0; c < ncol; ++c) {
                auto observed = dcwork->fetch(c);
                for (size_t r = 0; r < nrow; ++r) {
                    if (observed[r]) {
                        expanded_by_col.insert(expanded_by_col.end(), observed[r], r);
                        hits[c + 1] += observed[r];
                    }
                }
                hits[c + 1] += hits[c];
            }
            sparse_column.reset(new tatami::SemiCompressedSparseColumnMatrix<double, int>(nrow, ncol, std::move(expanded_by_col), std::move(hits)));
        }

        // And again, in the other direction.
        {
            std::vector<int> expanded_by_row;
            std::vector<size_t> hits(nrow + 1);
            auto drwork = dense->dense_row();
            for (size_t r = 0; r < nrow; ++r) {
                auto observed = drwork->fetch(r);
                for (size_t c = 0; c < ncol; ++c) {
                    if (observed[c]) {
                        expanded_by_row.insert(expanded_by_row.end(), observed[c], c);
                        hits[r + 1] += observed[c];
                    }
                }
                hits[r + 1] += hits[r];
            }
            sparse_row.reset(new tatami::SemiCompressedSparseRowMatrix<double, int>(nrow, ncol, std::move(expanded_by_row), std::move(hits)));
        }

        return;
    }
};

class SemiSparseUtilsTest : public ::testing::Test, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_F(SemiSparseUtilsTest, Basic) {
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

class SemiSparseFullAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t> >, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SemiSparseFullAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
//    tatami_test::test_simple_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

TEST_P(SemiSparseFullAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    tatami_test::test_simple_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    tatami_test::test_simple_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

INSTANTIATE_TEST_CASE_P(
    SemiCompressedSparseMatrix,
    SemiSparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

class SemiSparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SemiSparseSlicedAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    tatami_test::test_sliced_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(SemiSparseSlicedAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    tatami_test::test_sliced_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    tatami_test::test_sliced_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    SemiCompressedSparseMatrix,
    SemiSparseSlicedAccessTest,
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

class SemiSparseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SemiSparseIndexedAccessTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    tatami_test::test_indexed_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(SemiSparseIndexedAccessTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    tatami_test::test_indexed_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    tatami_test::test_indexed_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    SemiCompressedSparseMatrix,
    SemiSparseIndexedAccessTest,
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

TEST(SemiCompressedSparseMatrix, ExtremeCounting) {
    // Checking for correct counting of duplicates at the extremes.
    std::vector<int> indices {
        2, 2, 2, 8, 11, 11, 12,
        3, 4, 4, 7, 10, 10, 12, 15, 18, 18, 18, 18,
        0, 0, 0, 0,  1,  4,  8,  8,  8, 14, 14
    };

    std::vector<size_t> indptrs {
        0,
        7,
        19,
        30
    };

    tatami::SemiCompressedSparseColumnMatrix<double, int> mat(19, 3, indices, indptrs);

    auto wrk = mat.dense_row();
    {
        std::vector<double> temp{ 0, 0, 4 };
        EXPECT_EQ(wrk->fetch(0), temp);

        temp = std::vector<double>{ 0, 4, 0 };
        EXPECT_EQ(wrk->fetch(18), temp);

        temp = std::vector<double>{ 0, 0, 4 }; // jumping back to the start.
        EXPECT_EQ(wrk->fetch(0), temp);
    }
}
