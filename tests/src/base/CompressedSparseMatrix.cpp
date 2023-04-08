#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/arith_scalar_helpers.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../_tests/test_row_access.h"
#include "../_tests/test_column_access.h"
#include "../_tests/simulate_vector.h"

TEST(CompressedSparseMatrix, ConstructionEmpty) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr);
    EXPECT_TRUE(mat.sparse());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);

    // Comparing access for an empty matrix.
    tatami::DenseColumnMatrix<double, int> dense(10, 20, std::vector<double>(200));
    test_simple_column_access(&mat, &dense);
    test_simple_row_access(&mat, &dense);
}

/*************************************
 *************************************/

class SparseTestMethods {
protected:
    size_t nrow = 200, ncol = 100;
    std::shared_ptr<tatami::NumericMatrix> dense, sparse_row, sparse_column;

    void assemble() {
        dense.reset(new tatami::DenseRowMatrix<double, int>(nrow, ncol, simulate_sparse_vector<double>(nrow * ncol, 0.05)));
        sparse_row = tatami::convert_to_sparse<true>(dense.get());
        sparse_column = tatami::convert_to_sparse<false>(dense.get());
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
    EXPECT_EQ(NC, ncol);
    EXPECT_EQ(NR, nrow);
    EXPECT_EQ(sparse_row->ncol(), ncol);
    EXPECT_EQ(sparse_row->nrow(), nrow);

    EXPECT_FALSE(dense->sparse());
    EXPECT_TRUE(sparse_column->sparse());
    EXPECT_TRUE(sparse_row->sparse());

    EXPECT_FALSE(sparse_column->prefer_rows());
    EXPECT_TRUE(sparse_row->prefer_rows());

    auto cprefs = sparse_column->dimension_preference();
    EXPECT_TRUE(cprefs.first == 0);
    EXPECT_TRUE(cprefs.second > 0);

    auto rprefs = sparse_row->dimension_preference();
    EXPECT_TRUE(rprefs.first > 0);
    EXPECT_TRUE(rprefs.second == 0);
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

TEST_P(SparseFullAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    test_simple_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    test_simple_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

TEST_P(SparseFullAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    test_simple_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP);
    test_simple_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP);
}

TEST_P(SparseFullAccessTest, Details) {
    size_t NC = sparse_column->ncol(), NR = sparse_column->nrow();
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);

    auto work_row = sparse_row->sparse_column_workspace();
    auto fetch_offsets = [](const auto* ptr) -> std::vector<size_t> {
        return dynamic_cast<const tatami::CompressedSparseRowMatrix<double, int>::CompressedSparseSecondarySparseWorkspace*>(ptr)->core.current_indptrs;
    };
    auto work_col = sparse_column->sparse_column_workspace();
    auto work_row2 = sparse_row->sparse_column_workspace();

    for (size_t i = 0; i < NC; i += JUMP) {
        size_t c = (FORWARD ? i : NC - i - 1);

        // Checking that CSR extraction actually has an effect on the latest indptrs.
        std::vector<size_t> old_offsets = fetch_offsets(work_row.get());
        sparse_row->column(c, work_row.get());
        std::vector<size_t> new_offsets = fetch_offsets(work_row.get());

        if (!FORWARD || c != 0) {
            EXPECT_NE(old_offsets, new_offsets); 
        }

        // Double-checking the read-only sparse extractors.
        std::vector<double> outval(sparse_column->nrow());
        std::vector<int> outidx(sparse_column->nrow());

        auto x = sparse_column->column(c, outval.data(), outidx.data(), work_col.get());
        EXPECT_TRUE(x.number < NR);
        EXPECT_FALSE(outval.data()==x.value); // points to internal data.
        EXPECT_FALSE(outidx.data()==x.index);

        auto y = sparse_row->column(c, outval.data(), outidx.data(), work_row2.get());
        EXPECT_TRUE(y.number < NR);
        EXPECT_TRUE(outval.data()==y.value); // points to buffer.
        EXPECT_TRUE(outidx.data()==y.index);
        EXPECT_EQ(new_offsets, fetch_offsets(work_row2.get()));
    }
}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseFullAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 4, 10, 20) // jump, to test the workspace's memory.
    )
);

/*************************************
 *************************************/

class SparseSlicedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseSlicedAccessTest, Column) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, LAST = interval_info[1] * nrow;

    test_sliced_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

TEST_P(SparseSlicedAccessTest, Row) {
    auto param = GetParam(); 
    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, LAST = interval_info[1] * ncol;

    test_sliced_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
    test_sliced_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, LAST);
}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseSlicedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.51 }),
            std::vector<double>({ 0.25, 0.9 }), 
            std::vector<double>({ 0.63, 1 })
        )
    )
);

/*************************************
 *************************************/

class SparseIndexedAccessTest : public ::testing::TestWithParam<std::tuple<bool, size_t, std::vector<double> > >, public SparseTestMethods {
protected:
    void SetUp() {
        assemble();
        return;
    }
};

TEST_P(SparseIndexedAccessTest, Column) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * nrow, STEP = interval_info[1] * nrow;

    test_indexed_column_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_column_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

TEST_P(SparseIndexedAccessTest, Row) {
    auto param = GetParam(); 

    bool FORWARD = std::get<0>(param);
    size_t JUMP = std::get<1>(param);
    auto interval_info = std::get<2>(param);
    size_t FIRST = interval_info[0] * ncol, STEP = interval_info[1] * ncol;

    test_indexed_row_access(sparse_column.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
    test_indexed_row_access(sparse_row.get(), dense.get(), FORWARD, JUMP, FIRST, STEP);
}

INSTANTIATE_TEST_CASE_P(
    CompressedSparseMatrix,
    SparseIndexedAccessTest,
    ::testing::Combine(
        ::testing::Values(true, false), // iterate forward or back, to test the workspace's memory.
        ::testing::Values(1, 3), // jump, to test the workspace's memory.
        ::testing::Values(
            std::vector<double>({ 0, 0.05 }),
            std::vector<double>({ 0.2, 0.1 }), 
            std::vector<double>({ 0.7, 0.03 })
        )
    )
);

