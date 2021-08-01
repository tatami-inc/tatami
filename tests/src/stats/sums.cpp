#include <gtest/gtest.h>

#include <vector>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/sums.hpp"

#include "../data/data.h"

TEST(ComputingDimsums, RowSums) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_sparse(dense_row.get(), false);

    auto ref = tatami::row_sums(dense_row.get());
    std::vector<double> expected(sparse_nrow);
    for (size_t r = 0; r < sparse_nrow; ++r) {
        for (size_t c = 0; c < sparse_ncol; ++c) {
            expected[r] += sparse_matrix[c + r * sparse_ncol];
        }
    }
    EXPECT_EQ(ref, expected);

    EXPECT_EQ(ref, tatami::row_sums(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_sums(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_sums(sparse_column.get()));
}

TEST(ComputingDimsums, ColumnSums) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_sparse(dense_row.get(), false);

    auto ref = tatami::column_sums(dense_row.get());
    std::vector<double> expected(sparse_ncol);
    for (size_t c = 0; c < sparse_ncol; ++c) {
        for (size_t r = 0; r < sparse_nrow; ++r) {
            expected[c] += sparse_matrix[c + r * sparse_ncol];
        }
    }
    EXPECT_EQ(ref, expected);

    EXPECT_EQ(ref, tatami::column_sums(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_sums(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_sums(sparse_column.get()));
}
