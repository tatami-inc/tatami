#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/CompressedSparseMatrix.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

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
    tatami_test::test_simple_row_access(&mat, &dense);
    tatami_test::test_simple_column_access(&mat, &dense);

    // Same for row-major.
    indptr.resize(11);
    tatami::CompressedSparseRowMatrix<double, int> rmat(10, 20, values, indices, indptr);
    EXPECT_TRUE(rmat.sparse());
    EXPECT_EQ(rmat.nrow(), 10);
    EXPECT_EQ(rmat.ncol(), 20);
    tatami_test::test_simple_row_access(&rmat, &dense);
    tatami_test::test_simple_column_access(&rmat, &dense);
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

TEST(CompressedSparseMatrix, OddTypes) {
    // Checking for compilation warnings here when the interface and storage types are different.
    {
        std::vector<uint8_t> values;
        std::vector<uint16_t> indices;
        std::vector<uint64_t> indptr(11);
        tatami::CompressedSparseRowMatrix<double, int, decltype(values), decltype(indices), decltype(indptr)> rmat(10, 20, values, indices, indptr);
    }

    {
        std::vector<uint8_t> values;
        std::vector<uint32_t> indices; // Check for signed/unsigned comparisons with the interface index type.
        std::vector<uint64_t> indptr(11);
        tatami::CompressedSparseRowMatrix<double, int32_t, decltype(values), decltype(indices), decltype(indptr)> rmat(10, 20, values, indices, indptr);
    }
}

/*************************************
 *************************************/

class SparseUtils {
protected:
    inline static size_t nrow = 200, ncol = 100;
    inline static std::shared_ptr<tatami::NumericMatrix> dense, sparse_row, sparse_column;

    static void assemble() {
        if (dense) {
            return;
        }
        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, tatami_test::simulate_sparse_vector<double>(nrow * ncol, 0.05)));
        sparse_row = tatami::convert_to_compressed_sparse<true>(dense.get());
        sparse_column = tatami::convert_to_compressed_sparse<false>(dense.get());
    }
};

class SparseTest : public ::testing::Test, public SparseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_F(SparseTest, Basic) {
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

    // Check that we actually get sparse output, and it doesn't just defer to
    // the dense methods (i.e., the number of extracted values < NR * NC).
    {
        auto sext = sparse_column->sparse_column();
        std::vector<double> values(NR);
        std::vector<int> indices(NR);
        int collected = 0;
        for (size_t c = 0; c < NC; ++c) {
            auto out = sext->fetch(c, values.data(), indices.data());
            collected += out.number;
        }
        EXPECT_TRUE(collected < NR * NC);
    }

    {
        auto sext = sparse_column->sparse_row();
        std::vector<double> values(NC);
        std::vector<int> indices(NC);
        int collected = 0;
        for (size_t r = 0; r < NR; ++r) {
            auto out = sext->fetch(r, values.data(), indices.data());
            collected += out.number;
        }
        EXPECT_TRUE(collected < NR * NC);
    }
}

/*************************************
 *************************************/

class SparseFullAccessTest : public ::testing::TestWithParam<tatami_test::StandardTestAccessParameters>, public SparseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(SparseFullAccessTest, Full) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(tparam);
    tatami_test::test_full_access(params, sparse_column.get(), dense.get());
    tatami_test::test_full_access(params, sparse_row.get(), dense.get());
}

INSTANTIATE_TEST_SUITE_P(
    CompressedSparseMatrix,
    SparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // row extraction.
        ::testing::Values(true, false), // an oracle.
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM),
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory; trying out longer jumps to check secondaries.
    )
);

/*************************************
 *************************************/

class SparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, double> > >, public SparseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(SparseSlicedAccessTest, Sliced) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(std::get<0>(tparam));

    auto interval_info = std::get<1>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, LAST = interval_info.second * len;

    tatami_test::test_block_access(params, sparse_column.get(), dense.get(), FIRST, LAST);
    tatami_test::test_block_access(params, sparse_row.get(), dense.get(), FIRST, LAST);
}

INSTANTIATE_TEST_SUITE_P(
    CompressedSparseMatrix,
    SparseSlicedAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 0.51),
            std::make_pair(0.25, 0.9), 
            std::make_pair(0.63, 1.0)
        )
    )
);

/*************************************
 *************************************/

class SparseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<tatami_test::StandardTestAccessParameters, std::pair<double, int> > >, public SparseUtils {
protected:
    static void SetUpTestSuite() {
        assemble();
    }
};

TEST_P(SparseIndexedAccessTest, Indexed) {
    auto tparam = GetParam(); 
    auto params = tatami_test::convert_access_parameters(std::get<0>(tparam));

    auto interval_info = std::get<1>(tparam);
    auto len = (params.use_row ? ncol : nrow);
    size_t FIRST = interval_info.first * len, STEP = interval_info.second;

    tatami_test::test_indexed_access(params, sparse_column.get(), dense.get(), FIRST, STEP);
    tatami_test::test_indexed_access(params, sparse_row.get(), dense.get(), FIRST, STEP);
}

INSTANTIATE_TEST_SUITE_P(
    CompressedSparseMatrix,
    SparseIndexedAccessTest,
    ::testing::Combine(
        tatami_test::standard_test_access_parameter_combinations(),
        ::testing::Values(
            std::make_pair(0.0, 5),
            std::make_pair(0.2, 10), 
            std::make_pair(0.7, 3)
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
            auto extracted = tatami_test::fetch(wrk.get(), r, ncol);
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
    auto sparse_column = tatami::convert_to_compressed_sparse<false>(&dense);
    tatami_test::TestAccessParameters param;
    param.use_row = true;
    param.order = tatami_test::FORWARD;
    tatami_test::test_full_access(param, sparse_column.get(), &dense);
    param.order = tatami_test::REVERSE;
    tatami_test::test_full_access(param, sparse_column.get(), &dense);
}
