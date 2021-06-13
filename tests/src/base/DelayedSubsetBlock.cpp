#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedSubset.hpp"
#include "tatami/base/DelayedSubsetBlock.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

class SubsetBlockTestCore : public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense, sparse;
    std::shared_ptr<tatami::workspace> work_dense, work_sparse;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column-major.
        return;
    }

    void create_create_workspaces(bool row) {
        work_dense = dense->new_workspace(row);
        work_sparse = sparse->new_workspace(row);
        return;
    }
};

class SubsetBlockTest : public SubsetBlockTestCore {
protected:
    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

TEST_F(SubsetBlockTest, SubsetRowFullColumnAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 3, 4, 5, 6 };
    std::vector<double> buffer_full(dense->nrow());

    auto ref = tatami::make_DelayedSubset<0>(dense, sub);
    auto dense_block = tatami::make_DelayedSubsetBlock<0>(dense, 3, 7);
    auto sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, 3, 7);

    set_sizes(0, sub.size());
    EXPECT_EQ(sub.size(), dense_block->nrow());
    EXPECT_EQ(dense->ncol(), dense_block->ncol());

    EXPECT_TRUE(dense_block->prefer_rows());
    EXPECT_FALSE(sparse_block->prefer_rows());

    create_create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(ref->column(i, expected.data()));

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_block->column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_block->column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetBlockTest, SubsetRowSlicedColumnAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 12, 13, 14, 15, 16, 17, 18, 19 };
    std::vector<double> buffer_full(dense->nrow());

    auto ref = tatami::make_DelayedSubset<0>(dense, sub);
    auto dense_block = tatami::make_DelayedSubsetBlock<0>(dense, 12, 20);
    auto sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, 12, 20);

    size_t LEN = 4;
    size_t first = 0;

    create_create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, sub.size()));

        wipe_expected();
        fill_expected(ref->column(i, output.data(), first, last));

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_block->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_block->column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);

        first += 13;
        first %= sub.size();
    }
}

TEST_F(SubsetBlockTest, SubsetRowFullRowAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 7, 8, 9, 10, 11, 12 };
    std::vector<double> buffer_full(dense->nrow());

    auto ref = tatami::make_DelayedSubset<0>(dense, sub);
    auto dense_block = tatami::make_DelayedSubsetBlock<0>(dense, 7, 13);
    auto sparse_block = tatami::make_DelayedSubsetBlock<0>(sparse, 7, 13);

    set_sizes(0, dense->ncol());

    create_create_workspaces(true);
    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_expected();
        fill_expected(ref->row(i, expected.data()));

        wipe_output();
        fill_output(sparse_block->row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_block->row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_block->row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetBlockTest, SubsetColumnFullRowAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 0, 1, 2, 3, 4, 5 };
    std::vector<double> buffer_full(dense->ncol());

    auto ref = tatami::make_DelayedSubset<1>(dense, sub);
    auto dense_block = tatami::make_DelayedSubsetBlock<1>(dense, 0, 6);
    auto sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, 0, 6);

    set_sizes(0, sub.size());
    EXPECT_EQ(sub.size(), dense_block->ncol());
    EXPECT_EQ(dense->nrow(), dense_block->nrow());

    create_create_workspaces(true);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(ref->row(i, expected.data()));

        wipe_output();
        fill_output(sparse_block->row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_block->row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_block->row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetBlockTest, SubsetColumnSlicedRowAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 5, 6, 7, 8, 9 };
    std::vector<double> buffer_full(dense->ncol());

    auto ref = tatami::make_DelayedSubset<1>(dense, sub);
    auto dense_block = tatami::make_DelayedSubsetBlock<1>(dense, 5, 10);
    auto sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, 5, 10);

    create_create_workspaces(true);
    size_t LEN = 7;
    first = 0;

    for (size_t i = 0; i < dense->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, sub.size()));

        wipe_expected();
        fill_expected(ref->row(i, expected.data(), first, last));

        wipe_output();
        fill_output(sparse_block->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_block->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_block->row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);

        first += 11;
        first %= sub.size();
    }
}

TEST_F(SubsetBlockTest, SubsetColumnFullColumnAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 4, 5, 6, 7, 8 };
    std::vector<double> buffer_full(dense->ncol());

    auto ref = tatami::make_DelayedSubset<1>(dense, sub);
    auto dense_block = tatami::make_DelayedSubsetBlock<1>(dense, 4, 9);
    auto sparse_block = tatami::make_DelayedSubsetBlock<1>(sparse, 4, 9);

    set_sizes(0, sparse->nrow());

    create_create_workspaces(false);
    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_expected();
        fill_expected(ref->column(i, expected.data()));

        wipe_output();
        fill_output(sparse_block->column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_block->column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_block->sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_block->column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}


