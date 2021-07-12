#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedBind.hpp"
#include "tatami/base/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

const double MULT1 = 10, MULT2 = 1.5;

template<class PARAM>
class BindTest: public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse, bound;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column-major.
    }

    void extra_assemble(const PARAM& param) {
        std::shared_ptr<tatami::NumericMatrix> one, two;

        // Multiplying to make sure we're actually extracting from a different submatrix.
        if (std::get<0>(param)) {
            one = sparse;
        } else {
            one = dense;            
        }
        one = tatami::make_DelayedIsometricOp(one, tatami::DelayedMultiplyScalarHelper(MULT1));

        if (std::get<1>(param)) {
            two = sparse;
        } else {
            two = dense;            
        }
        two = tatami::make_DelayedIsometricOp(two, tatami::DelayedMultiplyScalarHelper(MULT2));

        if (std::get<2>(param)) {
            bound = tatami::make_DelayedBind<0>(std::vector{ one, two });
        } else {
            bound = tatami::make_DelayedBind<1>(std::vector{ one, two });
        }

        return;
    }

    std::vector<double> harvest_expected_row(size_t i, bool rbind) {
        std::vector<double> expected;
        if (rbind) {
            const size_t NR = dense->nrow();
            expected = this->template extract_dense<true>(dense.get(), i % NR);
            for (auto& e : expected) { e *= (i < NR ? MULT1 : MULT2); }
        } else {
            expected = this->template extract_dense<true>(dense.get(), i);
            auto expected2 = expected;
            for (auto& e : expected) { e *= MULT1; }
            for (auto& e : expected2) { e *= MULT2; }
            expected.insert(expected.end(), expected2.begin(), expected2.end());
        }
        return expected;
    }

    std::vector<double> harvest_expected_column(size_t i, bool rbind) {
        std::vector<double> expected;
        if (rbind) {
            expected = this->template extract_dense<false>(dense.get(), i);
            auto expected2 = expected;
            for (auto& e : expected) { e *= MULT1; }
            for (auto& e : expected2) { e *= MULT2; }
            expected.insert(expected.end(), expected2.begin(), expected2.end());
        } else {
            const size_t NC = dense->ncol();
            expected = this->template extract_dense<false>(dense.get(), i % NC);
            for (auto& e : expected) { e *= (i < NC ? MULT1 : MULT2); }
        }
        return expected;
    }
};

/****************************
 ****************************/

using BindFullAccessTest = BindTest<std::tuple<bool, bool, bool, size_t> >;

TEST_P(BindFullAccessTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);

    size_t NR = dense->nrow();
    size_t NC = dense->ncol();
    if (std::get<2>(param)) {
        EXPECT_EQ(bound->ncol(), dense->ncol());
        EXPECT_EQ(bound->nrow(), 2 * NR);
    } else {
        EXPECT_EQ(bound->nrow(), dense->nrow());
        EXPECT_EQ(bound->ncol(), 2 * NC);
    }

    if (std::get<0>(param)!=std::get<1>(param)) {
        // only true if we're combining dense (row major) with sparse (column-major).
        auto preference = bound->dimension_preference();
        EXPECT_EQ(preference.first, preference.second);
    }

    auto work_bound = bound->new_workspace(true);

    for (size_t i = 0; i < bound->nrow(); i += std::get<3>(param)) {
        auto expected = harvest_expected_row(i, std::get<2>(param));

        auto outputD = extract_dense<true>(bound.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<true>(bound.get(), i, work_bound.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<true>(bound.get(), i);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<true>(bound.get(), i, work_bound.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(BindFullAccessTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    size_t NR = dense->nrow();
    size_t NC = dense->ncol();
    auto work_bound = bound->new_workspace(false);

    for (size_t i = 0; i < bound->ncol(); i += std::get<3>(param)) {
        auto expected = harvest_expected_column(i, std::get<2>(param));

        auto outputD = extract_dense<false>(bound.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<false>(bound.get(), i, work_bound.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<false>(bound.get(), i);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<false>(bound.get(), i, work_bound.get());
        EXPECT_EQ(outputSW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedBind,
    BindFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // use dense or sparse for the first matrix.
        ::testing::Values(true, false), // use dense or sparse for the second matrix.
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(1, 3) // jump, to test the workspace's memory.
    )
);

/****************************
 ****************************/

using BindSubsetAccessTest = BindTest<std::tuple<bool, bool, bool, size_t, std::vector<size_t> > >;

TEST_P(BindSubsetAccessTest, RowAccess) {
    auto param = GetParam();
    extra_assemble(param);

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    auto work_bound = bound->new_workspace(true);

    for (size_t i = 0; i < bound->nrow(); i += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected_raw = harvest_expected_row(i, std::get<2>(param));
        std::vector<double> expected(expected_raw.begin() + start, expected_raw.begin() + end);

        auto outputD = extract_dense<true>(bound.get(), i, start, end);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<true>(bound.get(), i, start, end, work_bound.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<true>(bound.get(), i, start, end);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<true>(bound.get(), i, start, end, work_bound.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(BindSubsetAccessTest, ColumnAccess) {
    auto param = GetParam();
    extra_assemble(param);

    size_t NR = dense->nrow();
    size_t NC = dense->ncol();

    size_t JUMP = std::get<3>(param);
    auto interval_info = std::get<4>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    auto work_bound = bound->new_workspace(true);

    for (size_t i = 0; i < bound->ncol(); i += JUMP, FIRST += SHIFT) {
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected_raw = harvest_expected_column(i, std::get<2>(param));
        std::vector<double> expected(expected_raw.begin() + start, expected_raw.begin() + end);
      
        auto outputD = extract_dense<false>(bound.get(), i, start, end);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<false>(bound.get(), i, start, end, work_bound.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<false>(bound.get(), i, start, end);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<false>(bound.get(), i, start, end, work_bound.get());
        EXPECT_EQ(outputSW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    DelayedBind,
    BindSubsetAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // use dense or sparse for the first matrix.
        ::testing::Values(true, false), // use dense or sparse for the second matrix.
        ::testing::Values(true, false), // bind by row or by column
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);
