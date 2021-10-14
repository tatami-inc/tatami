#include <gtest/gtest.h>

#include <vector>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/variances.hpp"

#include "../data/data.h"

template<class L, class R>
void compare_double_vectors (const L& left, const R& right) {
    ASSERT_EQ(left.size(), right.size());
    for (size_t i = 0; i < left.size(); ++i) {
        EXPECT_FLOAT_EQ(left[i], right[i]);
    }
    return;
}

TEST(ComputingDimVariances, RowVariances) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_sparse(dense_row.get(), false);

    // Doing the difference of squares as a quick-and-dirty reference.
    auto ref = tatami::row_variances(dense_row.get());
    std::vector<double> expected(sparse_nrow), expectedm(sparse_nrow);
    for (size_t r = 0; r < sparse_nrow; ++r) {
        for (size_t c = 0; c < sparse_ncol; ++c) {
            double x = sparse_matrix[c + r * sparse_ncol];
            expectedm[r] += x;
            expected[r] += x * x;
        }
        expectedm[r] /= sparse_ncol;
        expected[r] /= sparse_ncol;
        expected[r] -= expectedm[r] * expectedm[r];
        expected[r] *= sparse_ncol;
        expected[r] /= sparse_ncol - 1;
    }
    compare_double_vectors(ref, expected);

    compare_double_vectors(ref, tatami::row_variances(dense_column.get()));
    compare_double_vectors(ref, tatami::row_variances(sparse_row.get()));
    compare_double_vectors(ref, tatami::row_variances(sparse_column.get()));
}

TEST(ComputingDimVariances, ColumnVariances) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_sparse(dense_row.get(), false);

    // Doing the difference of squares as a quick-and-dirty reference.
    auto ref = tatami::column_variances(dense_row.get());
    std::vector<double> expected(sparse_ncol), expectedm(sparse_ncol);
    for (size_t c = 0; c < sparse_ncol; ++c) {
        for (size_t r = 0; r < sparse_nrow; ++r) {
            double x = sparse_matrix[c + r * sparse_ncol];
            expectedm[c] += x;
            expected[c] += x * x;
        }
        expectedm[c] /= sparse_nrow;
        expected[c] /= sparse_nrow;
        expected[c] -= expectedm[c] * expectedm[c];
        expected[c] *= sparse_nrow;
        expected[c] /= sparse_nrow - 1;
    }
    compare_double_vectors(ref, expected);

    compare_double_vectors(ref, tatami::column_variances(dense_column.get()));
    compare_double_vectors(ref, tatami::column_variances(sparse_row.get()));
    compare_double_vectors(ref, tatami::column_variances(sparse_column.get()));
}

TEST(ComputingDimVariances, RowVariancesNaN) {
    auto copy = sparse_matrix;
    copy.resize(0);
    auto dense = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, 0, copy));

    auto cref = tatami::column_variances(dense.get());
    EXPECT_EQ(cref.size(), 0);
    
    auto rref = tatami::row_variances(dense.get());
    EXPECT_TRUE(std::isnan(rref.front()));
    EXPECT_TRUE(std::isnan(rref.back()));
}

TEST(ComputingDimVariances, Configuration) {
    typedef tatami::stats::VarianceFactory<double> VarFact;
    EXPECT_TRUE(tatami::stats::has_sparse_running<VarFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_running_parallel<VarFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running<VarFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running_parallel<VarFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_direct<VarFact>::value);

    typedef decltype(std::declval<VarFact>().dense_direct()) VarDense;
    const bool ndc = tatami::stats::has_nonconst_dense_compute<VarDense, double, int>::value;
    EXPECT_FALSE(ndc);
    typedef decltype(std::declval<VarFact>().sparse_direct()) VarSparse;
    const bool nsc = tatami::stats::has_nonconst_sparse_compute<VarSparse, double, int>::value;
    EXPECT_FALSE(nsc);
    const tatami::SparseCopyMode nscc = tatami::stats::nonconst_sparse_compute_copy_mode<VarSparse>::value;
    EXPECT_EQ(nscc, tatami::SPARSE_COPY_BOTH); // just a negative control.
}
