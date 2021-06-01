#include <gtest/gtest.h>

#include "tatami/CompressedSparseMatrix.hpp"
#include "tatami/DenseMatrix.hpp"
#include "tatami/DelayedSubsetOp.hpp"

#include "data.h"
#include "load_sparse.h"
#include "TestCore.h"

#include <vector>
#include <memory>

class SubsetTestCore : public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense;
    std::shared_ptr<tatami::typed_matrix<double, int> > sparse;
    std::shared_ptr<tatami::workspace> work_dense, work_sparse;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = load_matrix_as_sparse_column_matrix(nr, nc, source);
        return;
    }

    void create_create_workspaces(bool row) {
        work_dense = dense->new_workspace(row);
        work_sparse = sparse->new_workspace(row);
        return;
    }
};

class SubsetTest : public SubsetTestCore {
protected:
    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

TEST_F(SubsetTest, SubsetRowFullColumnAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 0, 3, 3, 13, 5, 2, 19, 4, 6, 11, 19, 8 };
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = tatami::DelayedSubsetOp<double, 0>(dense, sub);
    auto sparse_subbed = tatami::DelayedSubsetOp<double, 0>(sparse, sub);

    set_sizes(0, sub.size());
    EXPECT_EQ(sub.size(), dense_subbed.nrow());
    EXPECT_EQ(dense->ncol(), dense_subbed.ncol());

    EXPECT_TRUE(dense->prefer_rows());
    EXPECT_FALSE(sparse->prefer_rows());

    create_create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.column(i, output.data()));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->column(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s : sub) {
            (*eIt) = X[s];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_column(i, outval.data(), outidx.data(), work_sparse));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.column(i, outval.data(), work_dense));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetRowSlicedColumnAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 17, 18, 11, 18, 15, 17, 13, 18, 11, 9, 6, 3, 6, 18, 1 };
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = tatami::DelayedSubsetOp<double, 0>(dense, sub);
    auto sparse_subbed = tatami::DelayedSubsetOp<double, 0>(sparse, sub);

    size_t LEN = 6;
    size_t first = 0;

    create_create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, sub.size()));

        wipe_output();
        fill_output(sparse_subbed.column(i, output.data(), first, last));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->column(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s = first; s < last; ++s) {
            (*eIt) = X[sub[s]];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse));
        EXPECT_EQ(output, expected);

        first += 13;
        first %= sub.size();
    }
}

TEST_F(SubsetTest, SubsetRowFullRowAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 13, 4, 17, 0, 17, 1, 19, 6, 1 };
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = tatami::DelayedSubsetOp<double, 0>(dense, sub);
    auto sparse_subbed = tatami::DelayedSubsetOp<double, 0>(sparse, sub);

    set_sizes(0, dense->ncol());

    create_create_workspaces(true);
    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.row(i, output.data()));

        wipe_expected();
        fill_expected(sparse->row(sub[i], expected.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_row(i, outval.data(), outidx.data(), work_sparse));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.row(i, outval.data(), work_dense));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetColumnFullRowAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 3, 9, 1, 0, 9, 5, 8, 3, 1, 8, 7 };
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = tatami::DelayedSubsetOp<double, 1>(dense, sub);
    auto sparse_subbed = tatami::DelayedSubsetOp<double, 1>(sparse, sub);

    set_sizes(0, sub.size());
    EXPECT_EQ(sub.size(), dense_subbed.ncol());
    EXPECT_EQ(dense->nrow(), dense_subbed.nrow());

    create_create_workspaces(true);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.row(i, output.data()));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->row(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s : sub) {
            (*eIt) = X[s];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_row(i, outval.data(), outidx.data(), work_sparse));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.row(i, outval.data(), work_dense));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetColumnSlicedRowAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 2, 2, 4, 8, 0, 7, 3, 1, 1, 2, 7, 8, 9, 9, 4, 5, 8, 5, 6, 2, 0 };
    size_t LEN = 7;
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = tatami::DelayedSubsetOp<double, 1>(dense, sub);
    auto sparse_subbed = tatami::DelayedSubsetOp<double, 1>(sparse, sub);

    create_create_workspaces(true);
    first = 0;
    for (size_t i = 0; i < dense->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, sub.size()));

        wipe_output();
        fill_output(sparse_subbed.row(i, output.data(), first, last));

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->row(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s = first; s < last && s < sub.size(); ++s) {
            (*eIt) = X[sub[s]];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse));
        EXPECT_EQ(output, expected);

        first += 11;
        first %= sub.size();
    }
}

TEST_F(SubsetTest, SubsetColumnFullColumnAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 7, 8, 0, 5, 1, 4, 1 };
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = tatami::DelayedSubsetOp<double, 1>(dense, sub);
    auto sparse_subbed = tatami::DelayedSubsetOp<double, 1>(sparse, sub);

    set_sizes(0, sparse->nrow());

    create_create_workspaces(false);
    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.column(i, output.data()));

        wipe_expected();
        fill_expected(sparse->column(sub[i], expected.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.sparse_column(i, outval.data(), outidx.data(), work_sparse));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.column(i, outval.data(), work_dense));
        EXPECT_EQ(output, expected);
    }
}


