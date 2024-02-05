#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/medians.hpp"

#include "tatami_test/tatami_test.hpp"

TEST(ComputeMedians, Dense) {
    std::vector<int> vec { 2, 1, 4, 5, 3 };
    EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size()), 3);
    EXPECT_EQ(tatami::stats::compute_median(vec.data() + 1, vec.size() - 1), 3.5);
}

TEST(ComputeMedians, Sparse) {
    {
        std::vector<int> vec { 2, 1, 4, 5, 3 };
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 5), 3);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 11), 0);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 10), 0.5);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 9), 1);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 8), 1.5);
    }

    {
        std::vector<int> vec { -2, -1, -4, -5, -3 };
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 5), -3);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 11), 0);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 10), -0.5);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 9), -1);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 8), -1.5);
    }

    // Various mixed flavors.
    {
        std::vector<double> vec { 2.5, -1, 4, -5, 3 };
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 5), 2.5);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 11), 0);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 10), 0);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 6), 1.25);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 7), 0);
    }

    {
        std::vector<double> vec { -2.5, 1, -4, 5, -3 };
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 5), -2.5);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 11), 0);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 10), 0);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 6), -1.25);
        EXPECT_EQ(tatami::stats::compute_median(vec.data(), vec.size(), 7), 0);
    }
}


TEST(ComputingDimMedians, SparseMedians) {
    size_t NR = 111, NC = 222;

    // We use a density of 0.5 so that we some of the median calculations will
    // need to use the structural zeros.  We also put all non-zero values on
    // one side of zero, otherwise the structural zeros will dominate the
    // median; in this case, we choose all-positive values.
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.5, 1, 10)));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto rref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(rref.size(), NR);
    EXPECT_EQ(rref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(rref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(rref, tatami::row_medians(sparse_column.get()));

    auto cref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(cref.size(), NC);
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

TEST(ComputingDimMedians, AllZero) {
    size_t NR = 55, NC = 22;
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, std::vector<double>(NR * NC)));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref, std::vector<double>(NR));
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref, std::vector<double>(NC));
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

/* Lots of additional checks necessary to account for the 
 * many conditional branches in the sparse case. We create
 * a triangular matrix to ensure that we get some good coverage
 * of all possible zero/non-zero combinations.
 */

class MedianTriangularTest : public ::testing::TestWithParam<int> {
protected:
    void triangularize(size_t order, std::vector<double>& values) {
        for (size_t r = 0; r < order; ++r) {
            for (size_t c = r + 1; c < order; ++c) {
                values[r * order + c] = 0; // wiping out the upper triangular.
            }
        }
    }
};

TEST_P(MedianTriangularTest, Positive) {
    size_t order = GetParam();
    auto dump = tatami_test::simulate_dense_vector<double>(order * order, 0.1, 1);
    triangularize(order, dump);

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(order, order, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), order);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), order);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST_P(MedianTriangularTest, Negative) {
    // Seeing what happens if all non-zeros are less than zero.
    size_t order = GetParam();
    auto dump = tatami_test::simulate_dense_vector<double>(order * order, -2, -0.1);
    triangularize(order, dump);

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(order, order, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), order);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), order);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

TEST_P(MedianTriangularTest, Mixed) {
    // Mixing up the ratios of non-zeros on both sides of zero.
    size_t order = GetParam();
    auto dump = tatami_test::simulate_dense_vector<double>(order * order, -2, 2);
    triangularize(order, dump);

    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(order, order, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::row_medians(dense_row.get());
    EXPECT_EQ(ref.size(), order);
    EXPECT_EQ(ref, tatami::row_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::row_medians(sparse_column.get()));

    ref = tatami::column_medians(dense_row.get());
    EXPECT_EQ(ref.size(), order);
    EXPECT_EQ(ref, tatami::column_medians(dense_column.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_row.get()));
    EXPECT_EQ(ref, tatami::column_medians(sparse_column.get()));
}

INSTANTIATE_TEST_SUITE_P(
    ComputingDimMedians,
    MedianTriangularTest,
    ::testing::Values(13, 22, 51, 80) // mix of even and odd numbers
);

TEST(ComputingDimMedians, RowMediansNaN) {
    auto dense = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(111, 0, std::vector<double>()));

    auto cref = tatami::column_medians(dense.get());
    EXPECT_EQ(cref.size(), 0);

    auto rref = tatami::row_medians(dense.get());
    EXPECT_TRUE(rref.size() > 0);
    EXPECT_TRUE(std::isnan(rref.front()));
    EXPECT_TRUE(std::isnan(rref.back()));
}

TEST(ComputingDimMedians, DirtyOutput) {
    size_t NR = 99, NC = 152;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.5, 1, 10); // see comments above about why we use 0.5.
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    auto ref = tatami::row_medians(dense_row.get());

    // Works when the input vector is a bit dirty.
    std::vector<double> dirty(NR, -1);
    tatami::row_medians(dense_row.get(), dirty.data());
    EXPECT_EQ(ref, dirty);

    std::fill(dirty.begin(), dirty.end(), -1);
    tatami::row_medians(dense_column.get(), dirty.data());
    EXPECT_EQ(ref, dirty);

    std::fill(dirty.begin(), dirty.end(), -1);
    tatami::row_medians(sparse_row.get(), dirty.data());
    EXPECT_EQ(ref, dirty);

    std::fill(dirty.begin(), dirty.end(), -1);
    tatami::row_medians(sparse_column.get(), dirty.data());
    EXPECT_EQ(ref, dirty);
}
