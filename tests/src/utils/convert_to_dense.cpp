#include <gtest/gtest.h>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"

#include "../_tests/simulate_vector.h"
#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"

TEST(ConvertToDense, RowToRow) {
    size_t NR = 20, NC = 10;
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<true>(&mat);
    EXPECT_TRUE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->new_row_workspace();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(converted2->row(i, wrk2.get()), expected2);
    }
}

TEST(ConvertToDense, ColumnToColumn) {
    size_t NR = 20, NC = 10;
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<false>(&mat);
    EXPECT_FALSE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<false, int, size_t>(&mat); // works for a different type.
    EXPECT_FALSE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->new_column_workspace();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(converted2->column(i, wrk2.get()), expected2);
    }
}

TEST(ConvertToDense, RowToColumn) {
    size_t NR = 70, NC = 50;
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<false>(&mat);
    EXPECT_FALSE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<false, int, size_t>(&mat); // works for a different type.
    EXPECT_FALSE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->new_row_workspace();
    for (size_t i = 0; i < NR; ++i) {
        auto start = vec.begin() + i * NC;
        std::vector<int> expected2(start, start + NC);
        EXPECT_EQ(converted2->row(i, wrk2.get()), expected2);
    }
}

TEST(ConvertToDense, ColumnToRow) {
    size_t NR = 70, NC = 50;
    auto vec = simulate_dense_vector<double>(NR * NC);
    tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);

    auto converted = tatami::convert_to_dense<true>(&mat);
    EXPECT_TRUE(converted->prefer_rows());
    EXPECT_FALSE(converted->sparse());
    test_simple_row_access(converted.get(), &mat);
    test_simple_column_access(converted.get(), &mat);

    auto converted2 = tatami::convert_to_dense<true, int, size_t>(&mat); // works for a different type.
    EXPECT_TRUE(converted2->prefer_rows());
    EXPECT_FALSE(converted2->sparse());

    auto wrk2 = converted2->new_column_workspace();
    for (size_t i = 0; i < NC; ++i) {
        auto start = vec.begin() + i * NR;
        std::vector<int> expected2(start, start + NR);
        EXPECT_EQ(converted2->column(i, wrk2.get()), expected2);
    }
}

TEST(ConvertToDense, Automatic) {
    size_t NR = 70, NC = 50;
    auto vec = simulate_dense_vector<double>(NR * NC);

    {
        tatami::DenseMatrix<false, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_dense(&mat, -1);
        EXPECT_FALSE(converted->prefer_rows());
    }

    {
        tatami::DenseMatrix<true, double, int> mat(NR, NC, vec);
        auto converted = tatami::convert_to_dense(&mat, -1);
        EXPECT_TRUE(converted->prefer_rows());
    }
}
