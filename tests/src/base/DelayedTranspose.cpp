#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

template<class PARAM>
class TransposeTest: public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse, tdense, tsparse;
protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse = tatami::convert_to_sparse<false>(dense.get()); // column-major.
        tdense = tatami::make_DelayedTranspose(dense);
        tsparse = tatami::make_DelayedTranspose(sparse);
    }
};

using TransposeFullTest = TransposeTest<int>;

TEST_P(TransposeFullTest, RowAccess) {
    EXPECT_EQ(tdense->ncol(), dense->nrow());
    EXPECT_EQ(tdense->nrow(), dense->ncol());
    EXPECT_EQ(!tdense->prefer_rows(), dense->prefer_rows());
    EXPECT_EQ(tdense->sparse(), dense->sparse());
    EXPECT_EQ(tsparse->sparse(), sparse->sparse());

    auto work_dense = tdense->new_workspace(false);
    auto work_sparse = tsparse ->new_workspace(false);

    size_t JUMP = GetParam();
    for (size_t i = 0; i < tdense->nrow(); i+=JUMP) {
        auto expected = extract_dense<false>(dense.get(), i);

        auto outputD = extract_dense<true>(tdense.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<true>(tdense.get(), i, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<true>(tsparse.get(), i);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<true>(tsparse.get(), i, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(TransposeFullTest, ColumnAccess) {
    auto work_dense = tdense->new_workspace(false);
    auto work_sparse = tsparse ->new_workspace(false);

    size_t JUMP = GetParam();
    for (size_t i = 0; i < tdense->ncol(); i+=JUMP) {
        auto expected = extract_dense<true>(dense.get(), i);

        auto outputD = extract_dense<false>(tdense.get(), i);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<false>(tdense.get(), i, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<false>(tsparse.get(), i);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<false>(tsparse.get(), i, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    TransposeTest,
    TransposeFullTest,
    ::testing::Values(1, 4) // jumps (to test workspace memory)
);

using TransposeSubsetTest = TransposeTest<std::tuple<int, std::vector<size_t> > >;

TEST_P(TransposeSubsetTest, RowAccess) {
    auto param = GetParam();
    size_t JUMP = std::get<0>(param);
    auto vec = std::get<1>(param);
    size_t FIRST = vec[0], LEN = vec[1], SHIFT = vec[2];

    auto work_dense = tdense->new_workspace(false);
    auto work_sparse = tsparse ->new_workspace(false);

    for (size_t i = 0; i < tdense->nrow(); i+=JUMP, FIRST += SHIFT) {
        auto interval_info = wrap_intervals(FIRST, FIRST + LEN, tdense->nrow());
        size_t first = interval_info.first, last = interval_info.second;

        auto expected = extract_dense<false>(dense.get(), i, first, last);

        auto outputD = extract_dense<true>(tdense.get(), i, first, last);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<true>(tdense.get(), i, first, last, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<true>(tsparse.get(), i, first, last);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<true>(tsparse.get(), i, first, last, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

TEST_P(TransposeSubsetTest, ColumnAccess) {
    auto param = GetParam();
    size_t JUMP = std::get<0>(param);
    auto vec = std::get<1>(param);
    size_t FIRST = vec[0], LEN = vec[1], SHIFT = vec[2];

    auto work_dense = tdense->new_workspace(false);
    auto work_sparse = tsparse ->new_workspace(false);

    for (size_t i = 0; i < tdense->ncol(); i+=JUMP, FIRST += SHIFT) {
        auto interval_info = wrap_intervals(FIRST, FIRST + LEN, tdense->nrow());
        size_t first = interval_info.first, last = interval_info.second;

        auto expected = extract_dense<true>(dense.get(), i, first, last);

        auto outputD = extract_dense<false>(tdense.get(), i, first, last);
        EXPECT_EQ(outputD, expected);

        auto outputDW = extract_dense<false>(tdense.get(), i, first, last, work_dense.get());
        EXPECT_EQ(outputDW, expected);

        auto outputS = extract_sparse<false>(tsparse.get(), i, first, last);
        EXPECT_EQ(outputS, expected);

        auto outputSW = extract_sparse<false>(tsparse.get(), i, first, last, work_sparse.get());
        EXPECT_EQ(outputSW, expected);
    }
}

INSTANTIATE_TEST_CASE_P(
    TransposeTest,
    TransposeSubsetTest,
    ::testing::Combine(
        ::testing::Values(1, 4), // jumps (to test workspace memory)
        ::testing::Values(
            std::vector<size_t>({ 0, 5, 3 }), // overlapping
            std::vector<size_t>({ 1, 7, 7 }), // non-overlapping
            std::vector<size_t>({ 5, 11, 0 }) // constant
        )
    )
);


