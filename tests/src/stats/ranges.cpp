#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/ranges.hpp"

#include "tatami_test/tatami_test.hpp"

class ComputingDimExtremesTest : public ::testing::TestWithParam<std::pair<double, double> > {
protected:
    static std::vector<double> simulate(size_t NR, size_t NC, std::pair<double, double> limits) {
        return tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1, limits.first, limits.second);
    }
};

TEST_P(ComputingDimExtremesTest, RowMins) {
    size_t NR = 39, NC = 42;
    auto dump = simulate(NR, NC, GetParam());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(NR);
    for (size_t r = 0; r < NR; ++r) {
        ref[r] = dump[r * NC];
        for (size_t c = 1; c < NC; ++c) {
            ref[r] = std::min(ref[r], dump[c + r * NC]);
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

TEST_P(ComputingDimExtremesTest, ColumnMins) {
    size_t NR = 55, NC = 122;
    auto dump = simulate(NR, NC, GetParam());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(NC);
    for (size_t c = 0; c < NC; ++c) {
        ref[c] = dump[c];
        for (size_t r = 1; r < NR; ++r) {
            ref[c] = std::min(ref[c], dump[c + r * NC]);
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

/********************************************/

TEST_P(ComputingDimExtremesTest, RowMaxs) {
    size_t NR = 125, NC = 32;
    auto dump = simulate(NR, NC, GetParam());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(NR);
    for (size_t r = 0; r < NR; ++r) {
        ref[r] = dump[r * NC];
        for (size_t c = 1; c < NC; ++c) {
            ref[r] = std::max(ref[r], dump[c + r * NC]);
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

TEST_P(ComputingDimExtremesTest, ColumnMaxs) {
    size_t NR = 115, NC = 42;
    auto dump = simulate(NR, NC, GetParam());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<double> ref(NC);
    for (size_t c = 0; c < NC; ++c) {
        ref[c] = dump[c];
        for (size_t r = 1; r < NR; ++r) {
            ref[c] = std::max(ref[c], dump[c + r * NC]);
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

TEST_P(ComputingDimExtremesTest, RowRanges) {
    size_t NR = 75, NC = 62;
    auto dump = simulate(NR, NC, GetParam());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
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

TEST_P(ComputingDimExtremesTest, ColumnRanges) {
    size_t NR = 111, NC = 52;
    auto dump = simulate(NR, NC, GetParam());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
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

/********************************************/

INSTANTIATE_TEST_CASE_P(
    ComputingDimExtremes,
    ComputingDimExtremesTest,
    ::testing::Values(
        std::make_pair(0.1, 10.0),   // only above.
        std::make_pair(-10.0, -0.4), // only below
        std::make_pair(-5.0, 5.0)    // mix of values above and below zero.
    ) 
);

/********************************************/

TEST(ComputingDimExtremes, AllZeros) {
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

TEST(ComputingDimExtremes, NoZeros) {
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

TEST(ComputingDimExtremes, Empty) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(10, 0, std::vector<double>()));
    EXPECT_EQ(tatami::column_mins(dense_row.get()).size(), 0);
    EXPECT_EQ(tatami::row_mins(dense_row.get()), std::vector<double>(10));
}

/********************************************/

TEST(ComputingDimExtremes, CrankyOracle) {
    size_t NR = 199, NC = 252;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto raw_dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_row = tatami_test::make_CrankyMatrix(raw_dense);
    auto dense_column = tatami_test::make_CrankyMatrix(tatami::convert_to_dense<false>(raw_dense.get()));

    {
        auto ref = tatami::column_mins(raw_dense.get());
        EXPECT_EQ(ref, tatami::column_mins(dense_row.get()));
        EXPECT_EQ(ref, tatami::column_mins(dense_row.get(), 2)); // works correctly with parallelization.
    }

    {
        auto ref = tatami::row_mins(raw_dense.get());
        EXPECT_EQ(ref, tatami::row_mins(dense_row.get()));
        EXPECT_EQ(ref, tatami::row_mins(dense_row.get(), 2));
    }
}
