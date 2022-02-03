#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/arith_scalar_helpers.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"
#include "../data/data.h"

TEST(CompressedSparseMatrix, ConstructionEmpty) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr);
    EXPECT_TRUE(mat.sparse());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
}

/*************************************
 *************************************/

class SparseTestMethods {
protected:
    std::shared_ptr<tatami::NumericMatrix> dense, sparse_row, sparse_column;

    void assemble() {
        dense = std::shared_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
        sparse_row = tatami::convert_to_sparse(dense.get(), true);
        sparse_column = tatami::convert_to_sparse(dense.get(), false);
        return;
    }
};

class SparseUtilsTest : public ::testing::Test, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_F(SparseUtilsTest, Basic) {
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

    // Checking that the workspace is null for column extraction of a CSC matrix.
    auto work_column = sparse_column->new_workspace(false);
    EXPECT_EQ(work_column.get(), nullptr);

    // ... and vice versa.
    auto work_row = sparse_row->new_workspace(true);
    EXPECT_EQ(work_row.get(), nullptr);
}

/*************************************
 *************************************/

class SparseFullAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t> >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseFullAccessTest, Basic) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    test_simple_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP);

    test_simple_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

TEST_P(SparseFullAccessTest, Details) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto work_row = sparse_row->new_workspace(false);
    auto fetch_offsets = [&]() -> std::vector<size_t> {
        return dynamic_cast<tatami::CompressedSparseRowMatrix<double, int>::CompressedSparseWorkspace*>(work_row.get())->current_indptrs;
    };

    for (size_t i = 0; i < NC; i += JUMP) {
        size_t c = (FORWARD ? i : NC - i - 1);

        // Checking that CSR extraction actually has an effect on the latest indptrs.
        std::vector<size_t> old_offsets = fetch_offsets();
        sparse_row->column(c, work_row.get());

        if (!FORWARD || c != 0) {
            std::vector<size_t> new_offsets = fetch_offsets();
            EXPECT_NE(old_offsets, new_offsets); 
        }

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

class SparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<size_t> > >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseSlicedAccessTest, Basic) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0], LEN = interval_info[1], SHIFT = interval_info[2];

    test_sliced_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LEN, SHIFT);
    test_sliced_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LEN, SHIFT);

    test_sliced_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LEN, SHIFT);
    test_sliced_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LEN, SHIFT);
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
