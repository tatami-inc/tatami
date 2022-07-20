#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/sums.hpp"

#include "../data/data.h"

TEST(ComputingDimsums, RowSums) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(sparse_nrow);
    for (size_t r = 0; r < sparse_nrow; ++r) {
        for (size_t c = 0; c < sparse_ncol; ++c) {
            ref[r] += sparse_matrix[c + r * sparse_ncol];
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

TEST(ComputingDimsums, ColumnSums) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(sparse_ncol);
    for (size_t c = 0; c < sparse_ncol; ++c) {
        for (size_t r = 0; r < sparse_nrow; ++r) {
            ref[c] += sparse_matrix[c + r * sparse_ncol];
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

TEST(ComputingDimsums, Configuration) {
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
    const tatami::SparseCopyMode nscc = tatami::stats::nonconst_sparse_compute_copy_mode<SumSparse>::value;
    EXPECT_EQ(nscc, tatami::SPARSE_COPY_BOTH); // just a negative control.
}
