#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/ranges.hpp"

#include "../data/data.h"

TEST(ComputingDimMins, RowMins) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(sparse_nrow);
    for (size_t r = 0; r < sparse_nrow; ++r) {
        ref[r] = sparse_matrix[r * sparse_ncol];
        for (size_t c = 1; c < sparse_ncol; ++c) {
            ref[r] = std::min(ref[r], sparse_matrix[c + r * sparse_ncol]);
        }
    }

    EXPECT_EQ(ref, tatami::row_mins(dense_row.get()));
    EXPECT_EQ(ref, tatami::row_mins(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_mins(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_mins(sparse_column.get()));

    // Same results from parallel code.
    EXPECT_EQ(ref, tatami::row_mins(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_mins(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::row_mins(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_mins(sparse_column.get(), 3));
}

TEST(ComputingDimMins, ColumnMins) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(sparse_ncol);
    for (size_t c = 0; c < sparse_ncol; ++c) {
        ref[c] = sparse_matrix[c];
        for (size_t r = 1; r < sparse_nrow; ++r) {
            ref[c] = std::min(ref[c], sparse_matrix[c + r * sparse_ncol]);
        }
    }

    EXPECT_EQ(ref, tatami::column_mins(dense_row.get()));
    EXPECT_EQ(ref, tatami::column_mins(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_mins(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_mins(sparse_column.get()));

    // Same results from parallel code.
    EXPECT_EQ(ref, tatami::column_mins(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::column_mins(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_mins(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::column_mins(sparse_column.get(), 3));
}

TEST(ComputingDimMins, Configuration) {
    typedef tatami::stats::MinFactory<double> MinFact;

    EXPECT_TRUE(tatami::stats::has_sparse_running<MinFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_running_parallel<MinFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running<MinFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running_parallel<MinFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_direct<MinFact>::value);

    typedef decltype(std::declval<MinFact>().dense_direct()) MinDense;
    const bool ndc = tatami::stats::has_nonconst_dense_compute<MinDense, double, int>::value;
    EXPECT_FALSE(ndc);
    typedef decltype(std::declval<MinFact>().sparse_direct()) MinSparse;
    const bool nsc = tatami::stats::has_nonconst_sparse_compute<MinSparse, double, int>::value;
    EXPECT_FALSE(nsc);
    const tatami::SparseCopyMode nscc = tatami::stats::nonconst_sparse_compute_copy_mode<MinSparse>::value;
    EXPECT_EQ(nscc, tatami::SPARSE_COPY_BOTH); // just a negative control.
}

/********************************************/

TEST(ComputingDimMaxs, RowMaxs) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(sparse_nrow);
    for (size_t r = 0; r < sparse_nrow; ++r) {
        ref[r] = sparse_matrix[r * sparse_ncol];
        for (size_t c = 1; c < sparse_ncol; ++c) {
            ref[r] = std::max(ref[r], sparse_matrix[c + r * sparse_ncol]);
        }
    }

    EXPECT_EQ(ref, tatami::row_maxs(dense_row.get()));
    EXPECT_EQ(ref, tatami::row_maxs(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_maxs(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_maxs(sparse_column.get()));

    // Same results from parallel code.
    EXPECT_EQ(ref, tatami::row_maxs(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_maxs(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::row_maxs(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_maxs(sparse_column.get(), 3));
}

TEST(ComputingDimMaxs, ColumnMaxs) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(sparse_ncol);
    for (size_t c = 0; c < sparse_ncol; ++c) {
        ref[c] = sparse_matrix[c];
        for (size_t r = 1; r < sparse_nrow; ++r) {
            ref[c] = std::max(ref[c], sparse_matrix[c + r * sparse_ncol]);
        }
    }

    EXPECT_EQ(ref, tatami::column_maxs(dense_row.get()));
    EXPECT_EQ(ref, tatami::column_maxs(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_maxs(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_maxs(sparse_column.get()));

    // Same results from parallel code.
    EXPECT_EQ(ref, tatami::column_maxs(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::column_maxs(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_maxs(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::column_maxs(sparse_column.get(), 3));
}

/********************************************/

TEST(ComputingDimRanges, RowRanges) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::row_ranges(dense_row.get());
    EXPECT_EQ(ref.first, tatami::row_mins(dense_row.get()));
    EXPECT_EQ(ref.second, tatami::row_maxs(dense_row.get()));

    EXPECT_EQ(ref, tatami::row_ranges(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_ranges(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_ranges(sparse_column.get()));

    // Same results from parallel code.
    EXPECT_EQ(ref, tatami::row_ranges(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_ranges(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::row_ranges(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::row_ranges(sparse_column.get(), 3));
}

TEST(ComputingDimRanges, ColumnRanges) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::column_ranges(dense_row.get());
    EXPECT_EQ(ref.first, tatami::column_mins(dense_row.get()));
    EXPECT_EQ(ref.second, tatami::column_maxs(dense_row.get()));

    EXPECT_EQ(ref, tatami::column_ranges(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_ranges(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_ranges(sparse_column.get()));

    // Same results from parallel code.
    EXPECT_EQ(ref, tatami::column_ranges(dense_row.get(), 3));
    EXPECT_EQ(ref, tatami::column_ranges(dense_column.get(), 3));
    EXPECT_EQ(ref, tatami::column_ranges(sparse_row.get(), 3));
    EXPECT_EQ(ref, tatami::column_ranges(sparse_column.get(), 3));
}

TEST(ComputingDimRanges, Configuration) {
    typedef tatami::stats::RangeFactory<double> RangeFact;

    EXPECT_TRUE(tatami::stats::has_sparse_running<RangeFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_running_parallel<RangeFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running<RangeFact>::value);
    EXPECT_TRUE(tatami::stats::has_dense_running_parallel<RangeFact>::value);
    EXPECT_TRUE(tatami::stats::has_sparse_direct<RangeFact>::value);

    typedef decltype(std::declval<RangeFact>().dense_direct()) RangeDense;
    const bool ndc = tatami::stats::has_nonconst_dense_compute<RangeDense, double, int>::value;
    EXPECT_FALSE(ndc);
    typedef decltype(std::declval<RangeFact>().sparse_direct()) RangeSparse;
    const bool nsc = tatami::stats::has_nonconst_sparse_compute<RangeSparse, double, int>::value;
    EXPECT_FALSE(nsc);
    const tatami::SparseCopyMode nscc = tatami::stats::nonconst_sparse_compute_copy_mode<RangeSparse>::value;
    EXPECT_EQ(nscc, tatami::SPARSE_COPY_BOTH); // just a negative control.
}

/********************************************/

TEST(ComputingDimRanges, AllZeros) {
    // Testing for correct sparse behavior with all-zeros.
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(10, 20, std::vector<double>(200)));
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto cref = std::make_pair(std::vector<double>(20), std::vector<double>(20));
    EXPECT_EQ(cref, tatami::column_ranges(dense_row.get()));
    EXPECT_EQ(cref, tatami::column_ranges(sparse_row.get()));
    EXPECT_EQ(cref, tatami::column_ranges(sparse_column.get()));

    auto rref = std::make_pair(std::vector<double>(10), std::vector<double>(10));
    EXPECT_EQ(rref, tatami::row_ranges(dense_row.get()));
    EXPECT_EQ(rref, tatami::row_ranges(sparse_row.get()));
    EXPECT_EQ(rref, tatami::row_ranges(sparse_row.get()));
}

TEST(ComputingDimRanges, NoZeros) {
    // Testing for correct behavior with no zeros.
    std::vector<double> stuff(200);
    for (size_t s = 0; s < stuff.size(); ++s) {
        stuff[s] = s + 1;
    }

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(10, 20, stuff));
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto cref = tatami::column_ranges(dense_row.get());
    EXPECT_EQ(cref, tatami::column_ranges(sparse_row.get()));
    EXPECT_EQ(cref, tatami::column_ranges(sparse_column.get()));

    auto rref = tatami::row_ranges(dense_row.get());
    EXPECT_EQ(rref, tatami::row_ranges(sparse_row.get()));
    EXPECT_EQ(rref, tatami::row_ranges(sparse_row.get()));
}

