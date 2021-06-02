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

template<class L, class R>
void compare_double_vectors (const L& left, const R& right) {
    ASSERT_EQ(left.size(), right.size());
    for (size_t i = 0; i < left.size(); ++i) {
        ASSERT_FLOAT_EQ(left[i], right[i]);
    }
    return;
}

TEST(ComputingDimVariances, RowVariances) {
    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = create_dense_column_major(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_row = load_matrix_as_sparse_row_matrix(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_column = load_matrix_as_sparse_column_matrix(sparse_nrow, sparse_ncol, sparse_matrix);

    auto ref = tatami::row_variances(dense_row.get());
    compare_double_vectors(ref, tatami::row_variances(dense_column.get()));
    compare_double_vectors(ref, tatami::row_variances(sparse_row.get()));
    compare_double_vectors(ref, tatami::row_variances(sparse_column.get()));
}

TEST(ComputingDimVariances, ColumnVariances) {
    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = create_dense_column_major(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_row = load_matrix_as_sparse_row_matrix(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_column = load_matrix_as_sparse_column_matrix(sparse_nrow, sparse_ncol, sparse_matrix);

    auto ref = tatami::column_variances(dense_row.get());
    compare_double_vectors(ref, tatami::column_variances(dense_column.get()));
    compare_double_vectors(ref, tatami::column_variances(sparse_row.get()));
    compare_double_vectors(ref, tatami::column_variances(sparse_column.get()));
}
