#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/arith_scalar_helpers.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

TEST(CompressedSparseMatrix, ConstructionEmpty) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr);
    EXPECT_TRUE(mat.sparse());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
    EXPECT_EQ(mat.type(), tatami::_double);
}

template<class PARAM>
class SparseTest : public TestCore<::testing::TestWithParam<PARAM> > {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense, sparse_row, sparse_column;
    std::shared_ptr<tatami::workspace> work_sparse_row, work_sparse_column;

protected:
    void SetUp() {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse_row = tatami::convert_to_sparse(dense.get(), true);
        sparse_column = tatami::convert_to_sparse(dense.get(), false);
        return;
    }
};

/*************************************
 *************************************/

using SparseFullAccessTest = SparseTest<std::tuple<bool, size_t> >;

TEST_P(SparseFullAccessTest, ColumnAccess) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();
    EXPECT_EQ(NC, sparse_ncol);
    EXPECT_EQ(NR, sparse_nrow);
    EXPECT_EQ(sparse_row->ncol(), sparse_ncol);
    EXPECT_EQ(sparse_row->nrow(), sparse_nrow);

    EXPECT_FALSE(dense->sparse());
    EXPECT_TRUE(sparse_column->sparse());
    EXPECT_TRUE(sparse_row->sparse());

    EXPECT_FALSE(sparse_column->prefer_rows());
    EXPECT_TRUE(sparse_row->prefer_rows());

    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto work_column = sparse_column->new_workspace(false);
    auto work_column2 = sparse_column->new_workspace(false);
    auto work_row = sparse_row->new_workspace(false);
    auto work_row2 = sparse_row->new_workspace(false);
    EXPECT_EQ(work_column.get(), nullptr);

    for (size_t i = 0; i < NC; i += JUMP) {
        size_t c = (FORWARD ? i : NC - i - 1);
        auto expected = extract_dense<false>(dense.get(), c);

        // Checking the CSC extractors.
        auto outputC = extract_dense<false>(sparse_column.get(), c);
        EXPECT_EQ(outputC, expected);

        auto outputCW = extract_dense<false>(sparse_column.get(), c, work_column.get());
        EXPECT_EQ(outputC, expected);

        auto outputCS = extract_sparse<false>(sparse_column.get(), c);
        EXPECT_EQ(outputCS, expected);

        auto outputCSW = extract_sparse<false>(sparse_column.get(), c, work_column2.get());
        EXPECT_EQ(outputCSW, expected);

        // Checking the CSR extractors.
        auto outputR = extract_dense<false>(sparse_row.get(), c);
        EXPECT_EQ(outputR, expected);

        std::vector<size_t> old_offsets = dynamic_cast<tatami::CompressedSparseRowMatrix<double, int>::compressed_sparse_workspace*>(work_row.get())->current_indptrs;

        auto outputRW = extract_dense<false>(sparse_row.get(), c, work_row.get());
        EXPECT_EQ(outputR, expected);

        if (!FORWARD || c!=0) {
            // Checking that it actually has an effect on the latest indptrs.
            std::vector<size_t> new_offsets = dynamic_cast<tatami::CompressedSparseRowMatrix<double, int>::compressed_sparse_workspace*>(work_row.get())->current_indptrs;
            EXPECT_NE(old_offsets, new_offsets); 
        }

        auto outputRS = extract_sparse<false>(sparse_row.get(), c);
        EXPECT_EQ(outputRS, expected);

        auto outputRSW = extract_sparse<false>(sparse_row.get(), c, work_row2.get());
        EXPECT_EQ(outputRSW, expected);

        // Double-checking the read-only sparse extractors.
        std::vector<double> outval(sparse_column->nrow());
        std::vector<int> outidx(sparse_column->nrow());

        auto x = sparse_column->sparse_column(i, outval.data(), outidx.data());
        EXPECT_TRUE(x.number < NR);
        EXPECT_FALSE(outval.data()==x.value); // points to internal data.
        EXPECT_FALSE(outidx.data()==x.index);

        auto y = sparse_row->sparse_column(i, outval.data(), outidx.data());
        EXPECT_TRUE(y.number < NR);
        EXPECT_TRUE(outval.data()==y.value); // points to buffer.
        EXPECT_TRUE(outidx.data()==y.index);
    }
}

TEST_P(SparseFullAccessTest, RowAccess) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();

    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto work_column = sparse_column->new_workspace(true);
    auto work_column2 = sparse_column->new_workspace(true);
    auto work_row = sparse_row->new_workspace(true);
    auto work_row2 = sparse_row->new_workspace(true);
    EXPECT_EQ(work_row.get(), nullptr);

    for (size_t i = 0; i < NR; i += JUMP) {
        size_t r = (FORWARD ? i : NR - i - 1);
        auto expected = extract_dense<true>(dense.get(), r);

        // Checking the CSC extractors.
        auto outputC = extract_dense<true>(sparse_column.get(), r);
        EXPECT_EQ(outputC, expected);

        auto outputCW = extract_dense<true>(sparse_column.get(), r, work_column.get());
        EXPECT_EQ(outputC, expected);

        auto outputCS = extract_sparse<true>(sparse_column.get(), r);
        EXPECT_EQ(outputCS, expected);

        auto outputCSW = extract_sparse<true>(sparse_column.get(), r, work_column2.get());
        EXPECT_EQ(outputCSW, expected);

        // Checking the CSR extractors.
        auto outputR = extract_dense<true>(sparse_row.get(), r);
        EXPECT_EQ(outputR, expected);

        auto outputRW = extract_dense<true>(sparse_row.get(), r, work_row.get());
        EXPECT_EQ(outputR, expected);

        auto outputRS = extract_sparse<true>(sparse_row.get(), r);
        EXPECT_EQ(outputRS, expected);

        auto outputRSW = extract_sparse<true>(sparse_row.get(), r, work_row2.get());
        EXPECT_EQ(outputRSW, expected);
    }

}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

using SparseSlicedAccessTest = SparseTest<std::tuple<bool, size_t, std::vector<size_t> > >;

TEST_P(SparseSlicedAccessTest, ColumnAccess) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();
    EXPECT_EQ(NC, sparse_ncol);
    EXPECT_EQ(NR, sparse_nrow);
    EXPECT_EQ(sparse_row->ncol(), sparse_ncol);
    EXPECT_EQ(sparse_row->nrow(), sparse_nrow);

    EXPECT_FALSE(dense->sparse());
    EXPECT_TRUE(sparse_column->sparse());
    EXPECT_TRUE(sparse_row->sparse());

    EXPECT_FALSE(sparse_column->prefer_rows());
    EXPECT_TRUE(sparse_row->prefer_rows());

    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    auto work_column = sparse_column->new_workspace(false);
    auto work_column2 = sparse_column->new_workspace(false);
    auto work_row = sparse_row->new_workspace(false);
    auto work_row2 = sparse_row->new_workspace(false);
    EXPECT_EQ(work_column.get(), nullptr);

    for (size_t i = 0; i < NC; i += JUMP, FIRST += SHIFT) {
        size_t c = (FORWARD ? i : NC - i - 1);
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<false>(dense.get(), c, start, end);

        // Checking the CSC extractors.
        auto outputC = extract_dense<false>(sparse_column.get(), c, start, end);
        EXPECT_EQ(outputC, expected);

        auto outputCW = extract_dense<false>(sparse_column.get(), c, start, end, work_column.get());
        EXPECT_EQ(outputC, expected);

        auto outputCS = extract_sparse<false>(sparse_column.get(), c, start, end);
        EXPECT_EQ(outputCS, expected);

        auto outputCSW = extract_sparse<false>(sparse_column.get(), c, start, end, work_column2.get());
        EXPECT_EQ(outputCSW, expected);

        // Checking the CSR extractors.
        auto outputR = extract_dense<false>(sparse_row.get(), c, start, end);
        EXPECT_EQ(outputR, expected);

        auto outputRW = extract_dense<false>(sparse_row.get(), c, start, end, work_row.get());
        EXPECT_EQ(outputR, expected);

        auto outputRS = extract_sparse<false>(sparse_row.get(), c, start, end);
        EXPECT_EQ(outputRS, expected);

        auto outputRSW = extract_sparse<false>(sparse_row.get(), c, start, end, work_row2.get());
        EXPECT_EQ(outputRSW, expected);
    }
}

TEST_P(SparseSlicedAccessTest, RowAccess) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();

    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    auto work_column = sparse_column->new_workspace(true);
    auto work_column2 = sparse_column->new_workspace(true);
    auto work_row = sparse_row->new_workspace(true);
    auto work_row2 = sparse_row->new_workspace(true);
    EXPECT_EQ(work_row.get(), nullptr);

    for (size_t i = 0; i < NR; i += JUMP) {
        size_t r = (FORWARD ? i : NR - i - 1);
        auto interval = wrap_intervals(FIRST, FIRST + LEN, dense->ncol());
        size_t start = interval.first, end = interval.second;

        auto expected = extract_dense<true>(dense.get(), r, start, end);

        // Checking the CSC extractors.
        auto outputC = extract_dense<true>(sparse_column.get(), r, start, end);
        EXPECT_EQ(outputC, expected);

        auto outputCW = extract_dense<true>(sparse_column.get(), r, start, end, work_column.get());
        EXPECT_EQ(outputC, expected);

        auto outputCS = extract_sparse<true>(sparse_column.get(), r, start, end);
        EXPECT_EQ(outputCS, expected);

        auto outputCSW = extract_sparse<true>(sparse_column.get(), r, start, end, work_column2.get());
        EXPECT_EQ(outputCSW, expected);

        // Checking the CSR extractors.
        auto outputR = extract_dense<true>(sparse_row.get(), r, start, end);
        EXPECT_EQ(outputR, expected);

        auto outputRW = extract_dense<true>(sparse_row.get(), r, start, end, work_row.get());
        EXPECT_EQ(outputR, expected);

        auto outputRS = extract_sparse<true>(sparse_row.get(), r, start, end);
        EXPECT_EQ(outputRS, expected);

        auto outputRSW = extract_sparse<true>(sparse_row.get(), r, start, end, work_row2.get());
        EXPECT_EQ(outputRSW, expected);
    }

}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<size_t>({ 0, 8, 3 }), // overlapping shifts
            std::vector<size_t>({ 1, 4, 4 }), // non-overlapping shifts
            std::vector<size_t>({ 3, 10, 0 })
        )
    )
);
