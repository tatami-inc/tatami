#include <gtest/gtest.h>

#include "../data.h"
#include "../load_sparse.h"
#include "tatami/tatami.h"
#include <vector>

template<class V>
std::unique_ptr<tatami::numeric_matrix> create_dense_column_major(size_t nr, size_t nc, const V& source) {
    auto copy = source;
    auto cIt = copy.begin();
    for (size_t c = 0; c < nc; ++c) {
        for (size_t r = 0; r < nr; ++r, ++cIt) {
            *cIt = source[c + r * nc];            
        }
    }
    return std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseColumnMatrix<double>(nr, nc, copy));
}

TEST(ComputingDimsums, RowSums) {
    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = create_dense_column_major(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_row = load_matrix_as_sparse_row_matrix(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_column = load_matrix_as_sparse_column_matrix(sparse_nrow, sparse_ncol, sparse_matrix);

    auto ref = tatami::rowsums(dense_row.get());
    EXPECT_EQ(ref, tatami::rowsums(dense_column.get()));
    EXPECT_EQ(ref, tatami::rowsums(sparse_row.get()));
    EXPECT_EQ(ref, tatami::rowsums(sparse_column.get()));
}

TEST(ComputingDimsums, ColumnSums) {
    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = create_dense_column_major(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_row = load_matrix_as_sparse_row_matrix(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_column = load_matrix_as_sparse_column_matrix(sparse_nrow, sparse_ncol, sparse_matrix);

    auto ref = tatami::colsums(dense_row.get());
    EXPECT_EQ(ref, tatami::colsums(dense_column.get()));
    EXPECT_EQ(ref, tatami::colsums(sparse_row.get()));
    EXPECT_EQ(ref, tatami::colsums(sparse_column.get()));
}
