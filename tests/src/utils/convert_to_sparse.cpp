#include <gtest/gtest.h>

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(ConvertToSparse, RowToRow) {
    size_t NR = 50, NC = 20;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_sparse<true>(&mat);
    EXPECT_TRUE(converted->sparse());
    EXPECT_TRUE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_sparse<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_TRUE(converted2->prefer_rows());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }

    // Parallelized.
    auto convertedP = tatami::convert_to_sparse<true>(&mat, 3);
    tatami_test::test_simple_row_access(convertedP.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(convertedP.get(), &mat, true, 1);
}

TEST(ConvertToSparse, ColumnToColumn) {
    size_t NR = 30, NC = 50;
    auto trip = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.1); // check sparse->sparse conversion with matching preferred dimension.
    tatami::CompressedSparseMatrix<false, double, int> mat(NR, NC, trip.value, trip.index, trip.ptr);

    auto converted = tatami::convert_to_sparse<false>(&mat);
    EXPECT_TRUE(converted->sparse());
    EXPECT_FALSE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_sparse<false, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_FALSE(converted2->prefer_rows());

    auto wrk = mat.dense_column();
    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto expected = wrk->fetch(i);
        std::vector<int> expected2(expected.begin(), expected.end());
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }
}

TEST(ConvertToSparse, RowToColumn) {
    size_t NR = 70, NC = 50;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.15);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_sparse<false>(&mat);
    EXPECT_TRUE(converted->sparse());
    EXPECT_FALSE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_sparse<false, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_FALSE(converted2->prefer_rows());

    auto wrk2 = converted2->dense_row();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }

    // Parallelized.
    auto convertedP = tatami::convert_to_sparse<false>(&mat, 3);
    tatami_test::test_simple_row_access(convertedP.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(convertedP.get(), &mat, true, 1);
}

TEST(ConvertToSparse, ColumnToRow) {
    size_t NR = 20, NC = 50;
    auto trip = tatami_test::simulate_sparse_compressed<double>(NC, NR, 0.15); // check sparse->sparse conversion with non-matching preferred dimension.
    tatami::CompressedSparseMatrix<false, double, int> mat(NR, NC, trip.value, trip.index, trip.ptr);

    auto converted = tatami::convert_to_sparse<true>(&mat);
    EXPECT_TRUE(converted->sparse());
    EXPECT_TRUE(converted->prefer_rows());
    tatami_test::test_simple_row_access(converted.get(), &mat, true, 1);
    tatami_test::test_simple_column_access(converted.get(), &mat, true, 1);

    auto converted2 = tatami::convert_to_sparse<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->sparse());
    EXPECT_TRUE(converted2->prefer_rows());

    auto wrk = mat.dense_column();
    auto wrk2 = converted2->dense_column();
    for (size_t i = 0; i < NC; ++i) {
        auto expected = wrk->fetch(i);
        std::vector<int> expected2(expected.begin(), expected.end());
        EXPECT_EQ(wrk2->fetch(i), expected2);
    }
}

TEST(ConvertToSparse, Automatic) {
    size_t NR = 70, NC = 50;
    auto vec = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.23);

    {
        tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_sparse(&mat, -1);
        EXPECT_FALSE(converted->prefer_rows());
    }

    {
        tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_sparse(&mat, -1);
        EXPECT_TRUE(converted->prefer_rows());
    }
}
