#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/dense/convert_to_dense.hpp"
#include "tatami/sparse/convert_to_compressed_sparse.hpp"
#include "tatami/stats/counts.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(ComputingDimCounts, RowNaNCounts) {
    size_t NR = 99, NC = 152;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    for (auto& x : dump) {
        if (x == 0) {
            x = std::numeric_limits<double>::quiet_NaN();
        }
    }

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_compressed_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_compressed_sparse<false>(dense_row.get());

    std::vector<int> ref(NR);
    for (size_t r = 0; r < NR; ++r) {
        for (size_t c = 0; c < NC; ++c) {
            ref[r] += std::isnan(dump[c + r * NC]);
        }
    }

    EXPECT_EQ(ref, tatami::row_nan_counts(dense_row.get()));
    EXPECT_EQ(ref, tatami::row_nan_counts(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_nan_counts(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_nan_counts(sparse_column.get()));

    // Checking same results from parallel code.
    EXPECT_EQ(ref, tatami::row_nan_counts(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_nan_counts(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::row_nan_counts(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_nan_counts(sparse_column.get(), 3));
}

TEST(ComputingDimCounts, ColumNaNCount) {
    size_t NR = 79, NC = 62;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    for (auto& x : dump) {
        if (x == 0) {
            x = std::numeric_limits<double>::quiet_NaN();
        }
    }

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_compressed_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_compressed_sparse<false>(dense_row.get());

    std::vector<int> ref(NC);
    for (size_t c = 0; c < NC; ++c) {
        for (size_t r = 0; r < NR; ++r) {
            ref[c] += std::isnan(dump[c + r * NC]);
        }
    }

    EXPECT_EQ(ref, tatami::column_nan_counts(dense_row.get()));
    EXPECT_EQ(ref, tatami::column_nan_counts(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_nan_counts(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_nan_counts(sparse_column.get()));

    // Checking same results from parallel code.
    EXPECT_EQ(ref, tatami::column_nan_counts(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_nan_counts(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_nan_counts(sparse_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_nan_counts(sparse_column.get(), 3));
}

TEST(ComputingDimCounts, RowZeroCounts) {
    size_t NR = 99, NC = 152;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_compressed_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_compressed_sparse<false>(dense_row.get());

    std::vector<int> ref(NR);
    for (size_t r = 0; r < NR; ++r) {
        for (size_t c = 0; c < NC; ++c) {
            ref[r] += dump[c + r * NC] == 0;
        }
    }

    EXPECT_EQ(ref, tatami::row_zero_counts(dense_row.get()));
    EXPECT_EQ(ref, tatami::row_zero_counts(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_zero_counts(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_zero_counts(sparse_column.get()));

    // Checking same results from parallel code.
    EXPECT_EQ(ref, tatami::row_zero_counts(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_zero_counts(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::row_zero_counts(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_zero_counts(sparse_column.get(), 3));
}

TEST(ComputingDimCounts, ColumZeroCount) {
    size_t NR = 79, NC = 62;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_compressed_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_compressed_sparse<false>(dense_row.get());

    std::vector<int> ref(NC);
    for (size_t c = 0; c < NC; ++c) {
        for (size_t r = 0; r < NR; ++r) {
            ref[c] += dump[c + r * NC] == 0;
        }
    }

    EXPECT_EQ(ref, tatami::column_zero_counts(dense_row.get()));
    EXPECT_EQ(ref, tatami::column_zero_counts(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_zero_counts(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_zero_counts(sparse_column.get()));

    // Checking same results from parallel code.
    EXPECT_EQ(ref, tatami::column_zero_counts(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_zero_counts(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_zero_counts(sparse_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_zero_counts(sparse_column.get(), 3));
}
