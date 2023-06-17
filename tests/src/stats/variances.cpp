#include <gtest/gtest.h>

#include <vector>

#ifdef CUSTOM_PARALLEL_TEST
// Put this before any tatami apply imports.
#include "custom_parallel.h"
#endif

#include "tatami/dense/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/variances.hpp"

#include "tatami_test/tatami_test.hpp"

template<class L, class R>
void compare_double_vectors (const L& left, const R& right) {
    ASSERT_EQ(left.size(), right.size());
    for (size_t i = 0; i < left.size(); ++i) {
        EXPECT_FLOAT_EQ(left[i], right[i]);
    }
    return;
}

TEST(ComputingDimVariances, RowVariances) {
    size_t NR = 109, NC = 82;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    // Doing the difference of squares as a quick-and-dirty reference.
    std::vector<double> ref(NR), expectedm(NR);
    for (size_t r = 0; r < NR; ++r) {
        for (size_t c = 0; c < NC; ++c) {
            double x = dump[c + r * NC];
            expectedm[r] += x;
            ref[r] += x * x;
        }
        expectedm[r] /= NC;
        ref[r] /= NC;
        ref[r] -= expectedm[r] * expectedm[r];
        ref[r] *= NC;
        ref[r] /= NC - 1;
    }

    compare_double_vectors(ref, tatami::row_variances(dense_row.get()));
    compare_double_vectors(ref, tatami::row_variances(dense_column.get()));
    compare_double_vectors(ref, tatami::row_variances(sparse_row.get()));
    compare_double_vectors(ref, tatami::row_variances(sparse_column.get()));

    // Same results from parallel code.
    compare_double_vectors(ref, tatami::row_variances(dense_row.get(), 3));
    compare_double_vectors(ref, tatami::row_variances(dense_column.get(), 3));
    compare_double_vectors(ref, tatami::row_variances(sparse_row.get(), 3));
    compare_double_vectors(ref, tatami::row_variances(sparse_column.get(), 3));
}

TEST(ComputingDimVariances, ColumnVariances) {
    size_t NR = 99, NC = 92;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_column = tatami::convert_to_dense<false>(dense_row.get());
    auto sparse_row = tatami::convert_to_sparse<true>(dense_row.get());
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    // Doing the difference of squares as a quick-and-dirty reference.
    std::vector<double> ref(NC), expectedm(NC);
    for (size_t c = 0; c < NC; ++c) {
        for (size_t r = 0; r < NR; ++r) {
            double x = dump[c + r * NC];
            expectedm[c] += x;
            ref[c] += x * x;
        }
        expectedm[c] /= NR;
        ref[c] /= NR;
        ref[c] -= expectedm[c] * expectedm[c];
        ref[c] *= NR;
        ref[c] /= NR - 1;
    }

    compare_double_vectors(ref, tatami::column_variances(dense_row.get()));
    compare_double_vectors(ref, tatami::column_variances(dense_column.get()));
    compare_double_vectors(ref, tatami::column_variances(sparse_row.get()));
    compare_double_vectors(ref, tatami::column_variances(sparse_column.get()));

    // Same results from parallel code.
    compare_double_vectors(ref, tatami::column_variances(dense_row.get(), 3));
    compare_double_vectors(ref, tatami::column_variances(dense_column.get(), 3));
    compare_double_vectors(ref, tatami::column_variances(sparse_row.get(), 3));
    compare_double_vectors(ref, tatami::column_variances(sparse_column.get(), 3));
}

TEST(ComputingDimVariances, RowVariancesNaN) {
    auto dense = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(111, 0, std::vector<double>()));

    auto cref = tatami::column_variances(dense.get());
    EXPECT_EQ(cref.size(), 0);
    
    auto rref = tatami::row_variances(dense.get());
    EXPECT_TRUE(rref.size() > 0);
    EXPECT_TRUE(std::isnan(rref.front()));
    EXPECT_TRUE(std::isnan(rref.back()));
}

TEST(RunningVariances, SensibleZeros) {
    size_t NR = 55, NC = 52;
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1)));
    auto sparse_column = tatami::convert_to_sparse<false>(dense_row.get());

    // We force the first (non-zero) value to be zero, and we check that 
    // the number of non-zeros is still correctly reported.
    std::vector<double> running_vars(NR);
    std::vector<double> running_means(NR);
    std::vector<int> running_nzeros(NR);
    {
        std::vector<int> ibuffer(NR);
        std::vector<double> vbuffer(NR);
        std::vector<int> ref_nzeros(NR);
        auto wrk = sparse_column->sparse_column();

        for (int c = 0; c < static_cast<int>(NC); ++c) {
            auto range = wrk->fetch_copy(c, vbuffer.data(), ibuffer.data());
            vbuffer[0] = 0; 
            tatami::stats::variances::compute_running(range, running_means.data(), running_vars.data(), running_nzeros.data(), c);
            for (size_t r = 1; r < range.number; ++r) {
                ref_nzeros[range.index[r]] += (range.value[r] != 0);
            }
        }
        tatami::stats::variances::finish_running(NR, running_means.data(), running_vars.data(), running_nzeros.data(), NC);

        EXPECT_EQ(ref_nzeros, running_nzeros);
    }

    // Unless we set skip_zeros = false.
    std::vector<double> running_vars2(NR);
    std::vector<double> running_means2(NR);
    std::vector<int> running_nzeros2(NR);
    {
        std::vector<int> ibuffer(NR);
        std::vector<double> vbuffer(NR);
        std::vector<int> ref_nzeros(NR);
        auto wrk = sparse_column->sparse_column();

        for (int c = 0; c < static_cast<int>(NC); ++c) {
            auto range = wrk->fetch_copy(c, vbuffer.data(), ibuffer.data());
            vbuffer[0] = 0; 
            tatami::stats::variances::compute_running(range, running_means2.data(), running_vars2.data(), running_nzeros2.data(), c, false);
            for (size_t r = 0; r < range.number; ++r) {
                ++ref_nzeros[range.index[r]];
            }
        }
        tatami::stats::variances::finish_running(NR, running_means2.data(), running_vars2.data(), running_nzeros2.data(), NC);

        EXPECT_EQ(ref_nzeros, running_nzeros2);
        for (size_t i = 0; i < NR; ++i) {
            EXPECT_FLOAT_EQ(running_means2[i], running_means[i]);
            EXPECT_FLOAT_EQ(running_vars2[i], running_vars[i]);
        }
    }
}

TEST(ComputingDimVariances, CrankyOracle) {
    size_t NR = 155, NC = 172;
    auto dump = tatami_test::simulate_sparse_vector<double>(NR * NC, 0.1);
    auto raw_dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(NR, NC, dump));
    auto dense_row = tatami_test::make_CrankyMatrix(raw_dense);

    {
        auto ref = tatami::column_variances(raw_dense.get());
        EXPECT_EQ(ref, tatami::column_variances(dense_row.get()));
        EXPECT_EQ(ref, tatami::column_variances(dense_row.get(), 2)); // works correctly when parallelized.
    }

    {
        auto ref = tatami::row_variances(raw_dense.get());
        EXPECT_EQ(ref, tatami::row_variances(dense_row.get()));
        EXPECT_EQ(ref, tatami::row_variances(dense_row.get(), 3));
    }
}
