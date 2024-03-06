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
    tatami_test::test_simple_column_access(&mat, &dense);
    tatami_test::test_simple_row_access(&mat, &dense);

    // Same for row-major.
    indptr.resize(11);
    tatami::SemiCompressedSparseRowMatrix<double, int> rmat(10, 20, indices, indptr);
    EXPECT_TRUE(rmat.sparse());
    EXPECT_EQ(rmat.nrow(), 10);
    EXPECT_EQ(rmat.ncol(), 20);
    tatami_test::test_simple_column_access(&rmat, &dense);
    tatami_test::test_simple_row_access(&rmat, &dense);
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
}

/*************************************
 *************************************/

class SemiSparseFullAccessTest : public ::testing::TestWithParam<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t> >, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SemiSparseFullAccessTest, Full) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters param;
    param.use_row = std::get<0>(tparam);
    param.use_oracle = std::get<1>(tparam);
    param.order = std::get<2>(tparam);
    param.jump = std::get<3>(tparam);

    tatami_test::test_full_access(param, sparse_column.get(), dense.get()); 
    tatami_test::test_full_access(param, sparse_row.get(), dense.get());
}

INSTANTIATE_TEST_SUITE_P(
    SemiCompressedSparseMatrix,
    SemiSparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

class SemiSparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t, std::vector<double> > >, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SemiSparseSlicedAccessTest, Sliced) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters param;
    param.use_row = std::get<0>(tparam);
    param.use_oracle = std::get<1>(tparam);
    param.order = std::get<2>(tparam);
    param.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (param.use_row ? ncol : nrow);
    size_t FIRST = interval_info[0] * len, LAST = interval_info[1] * len;

    tatami_test::test_block_access(param, sparse_column.get(), dense.get(), FIRST, LAST);
    tatami_test::test_block_access(param, sparse_row.get(), dense.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    SemiCompressedSparseMatrix,
    SemiSparseSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
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

class SemiSparseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<bool, bool, tatami_test::TestAccessOrder, size_t, std::vector<double> > >, public SemiSparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SemiSparseIndexedAccessTest, Indexed) {
    auto tparam = GetParam(); 

    tatami_test::TestAccessParameters param;
    param.use_row = std::get<0>(tparam);
    param.use_oracle = std::get<1>(tparam);
    param.order = std::get<2>(tparam);
    param.jump = std::get<3>(tparam);

    auto interval_info = std::get<4>(tparam);
    auto len = (param.use_row ? ncol : nrow);
    size_t FIRST = interval_info[0] * len, STEP = interval_info[1] * len;

    tatami_test::test_indexed_access(param, sparse_column.get(), dense.get(), FIRST, STEP);
    tatami_test::test_indexed_access(param, sparse_row.get(), dense.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    SemiCompressedSparseMatrix,
    SemiSparseIndexedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), 
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
