#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/medians.hpp"

#include "../data/data.h"

TEST(ComputingDimMedians, SparseMedians) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto rref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(rref.size(), sparse_nrow);
    EXPECT_EQ(rref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(rref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(rref, tatami::row_medians(sparse_column.get()));

    auto cref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(cref.size(), sparse_ncol);
    EXPECT_EQ(cref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(cref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(cref, tatami::column_medians(sparse_column.get()));

    // Checking that the parallel code is the same.
    EXPECT_EQ(rref, tatami::row_medians(dense_row.get(), 3));
    EXPECT_EQ(rref, tatami::row_medians(dense_column.get(), 3));
    EXPECT_EQ(rref, tatami::row_medians(sparse_row.get(), 3));
    EXPECT_EQ(rref, tatami::row_medians(sparse_column.get(), 3));

    EXPECT_EQ(cref, tatami::column_medians(dense_row.get(), 3));
    EXPECT_EQ(cref, tatami::column_medians(dense_column.get(), 3));
    EXPECT_EQ(cref, tatami::column_medians(sparse_row.get(), 3));
    EXPECT_EQ(cref, tatami::column_medians(sparse_column.get(), 3));
}

/* Lots of additional checks necessary to account for the 
 * many conditional branches in the sparse case.
 */

TEST(ComputingDimMedians, TriangularMedians) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, triangular_matrix_odd));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

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

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, copy));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

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

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, copy));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

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

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(triangular_order_odd, triangular_order_odd, copy));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

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

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(triangular_order_even, triangular_order_even, copy));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

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
    auto dense = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, 0, copy));

    auto cref = tatami::column_medians(dense.get());
    EXPECT_EQ(cref.size(), 0);

    auto rref = tatami::row_medians(dense.get());
    EXPECT_TRUE(std::isnan(rref.front()));
    EXPECT_TRUE(std::isnan(rref.back()));
}

TEST(ComputingDimMedians, Configuration) {
    typedef tatami::stats::MedianFactory<double> MedFact;
    EXPECT_FALSE(tatami::stats::has_sparse_running<MedFact>::value);
    EXPECT_FALSE(tatami::stats::has_sparse_running_parallel<MedFact>::value);
    EXPECT_FALSE(tatami::stats::has_dense_running<MedFact>::value);
    EXPECT_FALSE(tatami::stats::has_dense_running_parallel<MedFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_direct<MedFact>::value);

    typedef decltype(std::declval<MedFact>().dense_direct()) MedDense;
    const bool ndc = tatami::stats::has_nonconst_dense_compute<MedDense, double, int>::value;
    EXPECT_TRUE(ndc);
    typedef decltype(std::declval<MedFact>().sparse_direct()) MedSparse;
    const bool nsc = tatami::stats::has_nonconst_sparse_compute<MedSparse, double, int>::value;
    EXPECT_TRUE(nsc);
    const tatami::SparseCopyMode nscc = tatami::stats::nonconst_sparse_compute_copy_mode<MedSparse>::value;
    EXPECT_EQ(nscc, tatami::SPARSE_COPY_VALUE);
}
