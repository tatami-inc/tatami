#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/base/dense/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/sums.hpp"

#include "../_tests/test_column_access.h"
#include "../_tests/test_row_access.h"
#include "../_tests/simulate_vector.h"

TEST(ComputingDimSums, RowSums) {
    size_t NR = 99, NC = 152;
    auto dump = simulate_sparse_vector<double>(NR * NC, 0.1);
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
    auto dump = simulate_sparse_vector<double>(NR * NC, 0.1);
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

TEST(ComputingDimSums, Configuration) {
    typedef tatami::stats::SumFactory<double> SumFact;

    EXPECT_TRUE(tatami::stats::has_sparse_running<SumFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_running_parallel<SumFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running<SumFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running_parallel<SumFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_direct<SumFact>::value);

    typedef decltype(std::declval<SumFact>().dense_direct()) SumDense;
    const bool ndc = tatami::stats::has_nonconst_dense_compute<SumDense, double, int>::value;
    EXPECT_FALSE(ndc);
    typedef decltype(std::declval<SumFact>().sparse_direct()) SumSparse;
    const bool nsc = tatami::stats::has_nonconst_sparse_compute<SumSparse, double, int>::value;
    EXPECT_FALSE(nsc);
}
