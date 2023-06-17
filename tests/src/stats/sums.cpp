#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/sums.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(ComputingDimSums, RowSums) {
    size_t NR = 99, NC = 152;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(NR);
    for (size_t r = 0; r < NR; ++r) {
        for (size_t c = 0; c < NC; ++c) {
            ref[r] += dump[c + r * NC];
        }
    }

    EXPECT_EQ(ref, tatami::row_sums(dense_row.get()));
    EXPECT_EQ(ref, tatami::row_sums(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_sums(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_sums(sparse_column.get()));

    // Checking same results from parallel code.
    EXPECT_EQ(ref, tatami::row_sums(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_sums(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::row_sums(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_sums(sparse_column.get(), 3));
}

TEST(ComputingDimSums, ColumnSums) {
    size_t NR = 79, NC = 62;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(NC);
    for (size_t c = 0; c < NC; ++c) {
        for (size_t r = 0; r < NR; ++r) {
            ref[c] += dump[c + r * NC];
        }
    }

    EXPECT_EQ(ref, tatami::column_sums(dense_row.get()));
    EXPECT_EQ(ref, tatami::column_sums(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_sums(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_sums(sparse_column.get()));

    // Checking same results from parallel code.
    EXPECT_EQ(ref, tatami::column_sums(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_sums(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_sums(sparse_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_sums(sparse_column.get(), 3));
}

TEST(ComputingDimSums, CrankyOracle) {
    size_t NR = 199, NC = 102;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto raw_dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_row = tatami_test::make_CrankyMatrix(raw_dense);
    auto dense_column = tatami_test::make_CrankyMatrix(tatami::convert_to_dense<false>(raw_dense.get()));

    {
        auto ref = tatami::column_sums(raw_dense.get());
        EXPECT_EQ(ref, tatami::column_sums(dense_row.get()));
        EXPECT_EQ(ref, tatami::column_sums(dense_row.get(), 2)); // Works correctly when parallelized.
    }

    {
        auto ref = tatami::row_sums(raw_dense.get());
        EXPECT_EQ(ref, tatami::row_sums(dense_row.get()));
        EXPECT_EQ(ref, tatami::row_sums(dense_row.get(), 3));
    }
}
