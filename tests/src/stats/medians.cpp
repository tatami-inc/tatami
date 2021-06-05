#include <gtest/gtest.h>

#include "../data.h"
#include "../load_sparse.h"
#include "tatami/stats/medians.hpp"
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

TEST(ComputingDimMedians, SparseMedians) {
    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = create_dense_column_major(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_row = load_matrix_as_sparse_row_matrix(sparse_nrow, sparse_ncol, sparse_matrix);
    auto sparse_column = load_matrix_as_sparse_column_matrix(sparse_nrow, sparse_ncol, sparse_matrix);

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), sparse_nrow);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), sparse_ncol);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

/* Lots of additional checks necessary to account for the 
 * many conditional branches in the sparse case.
 */

TEST(ComputingDimMedians, TriangularMedians) {
    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, triangular_matrix_odd));
    auto dense_column = create_dense_column_major(triangular_order_odd, triangular_order_odd, triangular_matrix_odd);
    auto sparse_row = load_matrix_as_sparse_row_matrix(triangular_order_odd, triangular_order_odd, triangular_matrix_odd);
    auto sparse_column = load_matrix_as_sparse_column_matrix(triangular_order_odd, triangular_order_odd, triangular_matrix_odd);

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_odd);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_odd);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST(ComputingDimMedians, TriangularMediansZero) {
    auto copy = triangular_matrix_odd;
    for (auto& s : copy) { s = 0; }

    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, copy));
    auto dense_column = create_dense_column_major(triangular_order_odd, triangular_order_odd, copy);
    auto sparse_row = load_matrix_as_sparse_row_matrix(triangular_order_odd, triangular_order_odd, copy);
    auto sparse_column = load_matrix_as_sparse_column_matrix(triangular_order_odd, triangular_order_odd, copy);

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref, std::vector<double>(triangular_order_odd));
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref, std::vector<double>(triangular_order_odd));
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST(ComputingDimMedians, TriangularMediansNegative) {
    auto copy = triangular_matrix_odd;
    for (auto& s : copy) { s *= -1; }

    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, copy));
    auto dense_column = create_dense_column_major(triangular_order_odd, triangular_order_odd, copy);
    auto sparse_row = load_matrix_as_sparse_row_matrix(triangular_order_odd, triangular_order_odd, copy);
    auto sparse_column = load_matrix_as_sparse_column_matrix(triangular_order_odd, triangular_order_odd, copy);

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_odd);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_odd);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST(ComputingDimMedians, TriangularMediansMixed) {
    auto copy = triangular_matrix_odd;
    double mult = 1;
    for (auto& s : copy) { 
        s *= mult;
        mult *= -1;
    }

    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, copy));
    auto dense_column = create_dense_column_major(triangular_order_odd, triangular_order_odd, copy);
    auto sparse_row = load_matrix_as_sparse_row_matrix(triangular_order_odd, triangular_order_odd, copy);
    auto sparse_column = load_matrix_as_sparse_column_matrix(triangular_order_odd, triangular_order_odd, copy);

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_odd);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_odd);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST(ComputingDimMedians, TriangularMediansEven) {
    auto copy = triangular_matrix_even;
    double mult = 1;
    for (auto& s : copy) { 
        s *= mult;
        mult *= -1;
    }

    auto dense_row = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(triangular_order_even, triangular_order_even, copy));
    auto dense_column = create_dense_column_major(triangular_order_even, triangular_order_even, copy);
    auto sparse_row = load_matrix_as_sparse_row_matrix(triangular_order_even, triangular_order_even, copy);
    auto sparse_column = load_matrix_as_sparse_column_matrix(triangular_order_even, triangular_order_even, copy);

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_even);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), triangular_order_even);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST(ComputingDimMedians, RowMediansNaN) {
    auto copy = sparse_matrix;
    copy.resize(0);
    auto dense = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, 0, copy));

    auto cref = tatami::column_medians(dense.get());
    EXPECT_EQ(cref.size(), 0);

    auto rref = tatami::row_medians(dense.get());
    EXPECT_TRUE(std::isnan(rref.front()));
    EXPECT_TRUE(std::isnan(rref.back()));
}
