#include <gtest/gtest.h>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

class ConvertToCompressedSparseTest : public ::testing::TestWithParam<std::tuple<bool, int> > {};

TEST_P(ConvertToCompressedSparseTest, RowToRow) {
    auto param = GetParam();
    auto two_pass = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    size_t NR = 50, NC = 20;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto mat = std::make_shared<tatami::DenseMatrix<true, double, int> >(NR, NC, vec);

    auto converted = tatami::convert_to_compressed_sparse<true>(mat.get(), two_pass, nthreads);
    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->sparse());
    EXPECT_TRUE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_compressed_sparse<true, int, size_t>(mat.get(), two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_TRUE(converted2->prefer_rows());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NC), expected2);
    }
}

TEST_P(ConvertToCompressedSparseTest, ColumnToColumn) {
    auto param = GetParam();
    auto two_pass = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    size_t NR = 30, NC = 50;
    auto trip = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.1); // check sparse->sparse conversion with matching preferred dimension.
    auto mat = std::make_shared<tatami::CompressedSparseMatrix<false, double, int> >(NR, NC, trip.value, trip.index, trip.ptr);

    auto converted = tatami::convert_to_compressed_sparse<false>(mat.get(), two_pass, nthreads);
    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->sparse());
    EXPECT_FALSE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_compressed_sparse<false, int, size_t>(mat.get(), two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_FALSE(converted2->prefer_rows());

    auto wrk = mat->dense_column();
    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto expected = tatami_test::fetch(wrk.get(), static_cast<int>(i), NR);
        std::vector<int> expected2(expected.begin(), expected.end());
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NR), expected2);
    }
}

TEST_P(ConvertToCompressedSparseTest, RowToColumn) {
    auto param = GetParam();
    auto two_pass = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    size_t NR = 70, NC = 50;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.15);
    auto mat = std::make_shared<tatami::DenseMatrix<true, double, int> >(NR, NC, vec);

    auto converted = tatami::convert_to_compressed_sparse<false>(mat.get(), two_pass, nthreads);
    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->sparse());
    EXPECT_FALSE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_compressed_sparse<false, int, size_t>(mat.get(), two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_FALSE(converted2->prefer_rows());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NC), expected2);
    }
}

TEST_P(ConvertToCompressedSparseTest, ColumnToRow) {
    auto param = GetParam();
    auto two_pass = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    size_t NR = 20, NC = 50;
    auto trip = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.15); // check sparse->sparse conversion with non-matching preferred dimension.
    auto mat = std::make_shared<tatami::CompressedSparseMatrix<false, double, int> >(NR, NC, trip.value, trip.index, trip.ptr);

    auto converted = tatami::convert_to_compressed_sparse<true>(mat.get(), two_pass, nthreads);
    EXPECT_EQ(converted->nrow(), NR);
    EXPECT_EQ(converted->ncol(), NC);
    EXPECT_TRUE(converted->sparse());
    EXPECT_TRUE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), mat.get());
    tatami_test::test_simple_column_access(converted.get(), mat.get());

    auto converted2 = tatami::convert_to_compressed_sparse<true, int, size_t>(mat.get(), two_pass, nthreads); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_TRUE(converted2->prefer_rows());

    auto wrk = mat->dense_column();
    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto expected = tatami_test::fetch(wrk.get(), static_cast<int>(i), NR);
        std::vector<int> expected2(expected.begin(), expected.end());
        EXPECT_EQ(tatami_test::fetch(wrk2.get(), i, NR), expected2);
    }
}

TEST_P(ConvertToCompressedSparseTest, Automatic) {
    auto param = GetParam();
    auto two_pass = std::get<0>(param);
    auto nthreads = std::get<1>(param);

    size_t NR = 70, NC = 50;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.23);

    {
        tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_compressed_sparse(&mat, -1, two_pass, nthreads);
        EXPECT_FALSE(converted->prefer_rows());
    }

    {
        tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_compressed_sparse(&mat, -1, two_pass, nthreads);
        EXPECT_TRUE(converted->prefer_rows());
    }
}

INSTANTIATE_TEST_SUITE_P(
    ConvertToCompressedSparse,
    ConvertToCompressedSparseTest,
    ::testing::Combine(
        ::testing::Values(false, true),
        ::testing::Values(1, 3)
    )
);
