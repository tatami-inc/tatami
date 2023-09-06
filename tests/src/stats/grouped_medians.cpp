#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/subset/make_DelayedSubset.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/grouped_medians.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(GroupedMedians, ByRow) {
    size_t NR = 99, NC = 155;
    auto dense_row = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.5)));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<int> cgroups(NC);
    int ngroup = 3; 
    std::vector<std::vector<int> > subsets(ngroup);
    for (size_t c = 0; c < NC; ++c) {
        cgroups[c] = c % ngroup;
        subsets[cgroups[c]].push_back(c);
    }

    auto rref = tatami::row_medians_by_group(dense_row.get(), cgroups.data());
    EXPECT_EQ(rref.size(), NR * ngroup);

    std::vector<double> expected(rref.size());
    for (int g = 0; g < ngroup; ++g) {
        auto sub = tatami::make_DelayedSubset<1>(dense_row, subsets[g]);
        auto per_group = tatami::row_medians(sub.get());
        for (size_t r = 0; r < NR; ++r) {
            expected[g + r * ngroup] = per_group[r];
        }
    }
    EXPECT_EQ(expected, rref);

    EXPECT_EQ(rref, tatami::row_medians_by_group(dense_column.get(), cgroups.data()));
    EXPECT_EQ(rref, tatami::row_medians_by_group(sparse_row.get(), cgroups.data()));
    EXPECT_EQ(rref, tatami::row_medians_by_group(sparse_column.get(), cgroups.data()));

    // Checking that the parallel code is the same.
    EXPECT_EQ(rref, tatami::row_medians_by_group(dense_row.get(), cgroups.data(), 3));
    EXPECT_EQ(rref, tatami::row_medians_by_group(dense_column.get(), cgroups.data(), 3));
    EXPECT_EQ(rref, tatami::row_medians_by_group(sparse_row.get(), cgroups.data(), 3));
    EXPECT_EQ(rref, tatami::row_medians_by_group(sparse_column.get(), cgroups.data(), 3));
}

TEST(GroupedMedians, ByColumn) {
    size_t NR = 56, NC = 179;
    auto dense_row = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.5)));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    std::vector<int> rgroups(NR);
    int ngroup = 7; 
    std::vector<std::vector<int> > subsets(ngroup);
    for (size_t r = 0; r < NR; ++r) {
        rgroups[r] = r % ngroup;
        subsets[rgroups[r]].push_back(r);
    }

    auto cref = tatami::column_medians_by_group(dense_row.get(), rgroups.data());
    EXPECT_EQ(cref.size(), NC * ngroup);

    std::vector<double> expected(cref.size());
    for (int g = 0; g < ngroup; ++g) {
        auto sub = tatami::make_DelayedSubset<0>(dense_row, subsets[g]);
        auto per_group = tatami::column_medians(sub.get());
        for (size_t c = 0; c < NC; ++c) {
            expected[g + c * ngroup] = per_group[c];
        }
    }
    EXPECT_EQ(expected, cref);

    EXPECT_EQ(cref, tatami::column_medians_by_group(dense_column.get(), rgroups.data()));
    EXPECT_EQ(cref, tatami::column_medians_by_group(sparse_row.get(), rgroups.data()));
    EXPECT_EQ(cref, tatami::column_medians_by_group(sparse_column.get(), rgroups.data()));

    // Checking that the parallel code is the same.
    EXPECT_EQ(cref, tatami::column_medians_by_group(dense_row.get(), rgroups.data(), 3));
    EXPECT_EQ(cref, tatami::column_medians_by_group(dense_column.get(), rgroups.data(), 3));
    EXPECT_EQ(cref, tatami::column_medians_by_group(sparse_row.get(), rgroups.data(), 3));
    EXPECT_EQ(cref, tatami::column_medians_by_group(sparse_column.get(), rgroups.data(), 3));
}

TEST(GroupedMedians, EdgeCases) {
    tatami::DenseRowMatrix<double, int> empty1(0, 10, std::vector<double>());
    tatami::DenseRowMatrix<double, int> empty2(10, 0, std::vector<double>());

    std::vector<int> grouping { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0 };
    EXPECT_TRUE(tatami::row_medians_by_group(&empty1, grouping.data()).empty());
    EXPECT_TRUE(tatami::column_medians_by_group(&empty2, grouping.data()).empty());

    grouping.clear();
    EXPECT_TRUE(tatami::row_medians_by_group(&empty2, grouping.data()).empty());
    EXPECT_TRUE(tatami::column_medians_by_group(&empty1, grouping.data()).empty());
}

TEST(GroupedMedians, CrankyOracle) {
    size_t NR = 199, NC = 20;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.5);

    auto raw_dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_row = tatami_test::make_CrankyMatrix(raw_dense);
    auto dense_column = tatami_test::make_CrankyMatrix(tatami::convert_to_dense<false>(raw_dense.get()));

    auto raw_sparse = tatami::convert_to_sparse<true>(raw_dense.get());
    auto sparse_row = tatami_test::make_CrankyMatrix(raw_sparse);
    auto sparse_column = tatami_test::make_CrankyMatrix(tatami::convert_to_sparse<false>(raw_sparse.get()));

    {
        std::vector<int> grouping(NR);
        for (size_t i = 0; i < NR; ++i) {
            grouping[i] = i % 3;
        }

        auto ref = tatami::column_medians_by_group(raw_dense.get(), grouping.data());

        EXPECT_EQ(ref, tatami::column_medians_by_group(dense_row.get(), grouping.data()));
        EXPECT_EQ(ref, tatami::column_medians_by_group(dense_column.get(), grouping.data()));
        EXPECT_EQ(ref, tatami::column_medians_by_group(sparse_row.get(), grouping.data()));
        EXPECT_EQ(ref, tatami::column_medians_by_group(sparse_column.get(), grouping.data()));

        // Works correctly when parallelized.
        EXPECT_EQ(ref, tatami::column_medians_by_group(dense_row.get(), grouping.data(), 2)); 
        EXPECT_EQ(ref, tatami::column_medians_by_group(dense_column.get(), grouping.data(), 2)); 
        EXPECT_EQ(ref, tatami::column_medians_by_group(sparse_row.get(), grouping.data(), 2)); 
        EXPECT_EQ(ref, tatami::column_medians_by_group(sparse_column.get(), grouping.data(), 2)); 
    }

    {
        std::vector<int> grouping(NC);
        for (size_t i = 0; i < NC; ++i) {
            grouping[i] = i % 3;
        }

        auto ref = tatami::row_medians_by_group(raw_dense.get(), grouping.data());

        EXPECT_EQ(ref, tatami::row_medians_by_group(dense_row.get(), grouping.data()));
        EXPECT_EQ(ref, tatami::row_medians_by_group(dense_column.get(), grouping.data()));
        EXPECT_EQ(ref, tatami::row_medians_by_group(sparse_row.get(), grouping.data()));
        EXPECT_EQ(ref, tatami::row_medians_by_group(sparse_column.get(), grouping.data()));

        // Works correctly when parallelized.
        EXPECT_EQ(ref, tatami::row_medians_by_group(dense_row.get(), grouping.data(), 2)); 
        EXPECT_EQ(ref, tatami::row_medians_by_group(dense_column.get(), grouping.data(), 2)); 
        EXPECT_EQ(ref, tatami::row_medians_by_group(sparse_row.get(), grouping.data(), 2)); 
        EXPECT_EQ(ref, tatami::row_medians_by_group(sparse_column.get(), grouping.data(), 2)); 
    }
}
