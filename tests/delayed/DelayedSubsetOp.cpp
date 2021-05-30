#include <gtest/gtest.h>

#include "matrix/CompressedSparseMatrix.hpp"
#include "matrix/DenseMatrix.hpp"
#include "delayed/DelayedSubsetOp.hpp"

#include "../utils/data.h"
#include "../utils/load_sparse.h"
#include "../utils/TestCore.h"

#include <vector>
#include <memory>

class SubsetTestCore : public TestCore {
protected:
    std::shared_ptr<bioc::numeric_matrix> dense;
    std::shared_ptr<bioc::typed_matrix<double, int> > sparse;
    std::unique_ptr<bioc::workspace> work_dense, work_sparse;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<bioc::numeric_matrix>(new bioc::DenseRowMatrix<double>(nr, nc, source));
        sparse = load_matrix_as_sparse_column_matrix(nr, nc, source);
        return;
    }

    void create_workspaces() {
        work_dense.reset(dense->create_workspace());
        work_sparse.reset(sparse->create_workspace());
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
    set_sizes(sub.size());
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = bioc::DelayedSubsetOp<double, 0>(dense, sub);
    auto sparse_subbed = bioc::DelayedSubsetOp<double, 0>(sparse, sub);

    first = 0;
    last = sub.size();
    EXPECT_EQ(sub.size(), dense_subbed.nrow());
    EXPECT_EQ(dense->ncol(), dense_subbed.ncol());

    create_workspaces();
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.get_column(i, output.data()), sub.size());

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->get_column(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s : sub) {
            (*eIt) = X[s];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.get_column(i, output.data()), sub.size());
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_column(i, outval.data(), outidx.data()), sub.size());
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_column(i, outval.data(), outidx.data(), work_sparse.get()), sub.size());
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.get_column(i, outval.data(), work_dense.get()), sub.size());
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetRowSlicedColumnAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 17, 18, 11, 18, 15, 17, 13, 18, 11, 9, 6, 3, 6, 18, 1 };
    size_t LEN = 6;
    set_sizes(LEN);
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = bioc::DelayedSubsetOp<double, 0>(dense, sub);
    auto sparse_subbed = bioc::DelayedSubsetOp<double, 0>(sparse, sub);

    create_workspaces();
    for (size_t i = 0; i < dense->ncol(); ++i) {
        first = i % sub.size();
        last = first + LEN;

        wipe_output();
        fill_output(sparse_subbed.get_column(i, output.data(), first, last), sub.size());

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->get_column(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s = first; s < last && s < sub.size(); ++s) {
            (*eIt) = X[sub[s]];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.get_column(i, output.data(), first, last), sub.size());
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_column(i, outval.data(), outidx.data(), first, last), sub.size());
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()), sub.size());
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetRowFullRowAccess) {
    // Column subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 13, 4, 17, 0, 17, 1, 19, 6, 1 };
    std::vector<double> buffer_full(dense->nrow());

    auto dense_subbed = bioc::DelayedSubsetOp<double, 0>(dense, sub);
    auto sparse_subbed = bioc::DelayedSubsetOp<double, 0>(sparse, sub);

    first = 0;
    size_t NC = sparse->ncol();
    last = NC;
    set_sizes(NC);

    create_workspaces();
    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.get_row(i, output.data()), sub.size());

        wipe_expected();
        fill_expected(sparse->get_row(sub[i], expected.data()), NC);
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.get_row(i, output.data()), NC);
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_row(i, outval.data(), outidx.data()), NC);
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_row(i, outval.data(), outidx.data(), work_sparse.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.get_row(i, outval.data(), work_dense.get()), NC);
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetColumnFullRowAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 3, 9, 1, 0, 9, 5, 8, 3, 1, 8, 7 };
    set_sizes(sub.size());
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = bioc::DelayedSubsetOp<double, 1>(dense, sub);
    auto sparse_subbed = bioc::DelayedSubsetOp<double, 1>(sparse, sub);

    first = 0;
    last = sub.size();
    EXPECT_EQ(sub.size(), dense_subbed.ncol());
    EXPECT_EQ(dense->nrow(), dense_subbed.nrow());

    create_workspaces();
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.get_row(i, output.data()), sub.size());

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->get_row(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s : sub) {
            (*eIt) = X[s];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.get_row(i, output.data()), sub.size());
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_row(i, outval.data(), outidx.data()), sub.size());
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_row(i, outval.data(), outidx.data(), work_sparse.get()), sub.size());
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.get_row(i, outval.data(), work_dense.get()), sub.size());
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetColumnSlicedRowAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 2, 2, 4, 8, 0, 7, 3, 1, 1, 2, 7, 8, 9, 9, 4, 5, 8, 5, 6, 2, 0 };
    size_t LEN = 7;
    set_sizes(LEN);
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = bioc::DelayedSubsetOp<double, 1>(dense, sub);
    auto sparse_subbed = bioc::DelayedSubsetOp<double, 1>(sparse, sub);

    create_workspaces();
    for (size_t i = 0; i < dense->nrow(); ++i) {
        first = i % sub.size();
        last = first + LEN;

        wipe_output();
        fill_output(sparse_subbed.get_row(i, output.data(), first, last), sub.size());

        // Reference extraction, the simple way.
        wipe_expected();
        auto X = sparse->get_row(i, buffer_full.data());
        auto eIt = expected.begin();
        for (auto s = first; s < last && s < sub.size(); ++s) {
            (*eIt) = X[sub[s]];
            ++eIt;
        }
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.get_row(i, output.data(), first, last), sub.size());
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_row(i, outval.data(), outidx.data(), first, last), sub.size());
        EXPECT_EQ(output, expected);        

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()), sub.size());
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SubsetTest, SubsetColumnFullColumnAccess) {
    // Row subsetting by a vector with duplicates, out of order.
    std::vector<size_t> sub = { 7, 8, 0, 5, 1, 4, 1 };
    std::vector<double> buffer_full(dense->ncol());

    auto dense_subbed = bioc::DelayedSubsetOp<double, 1>(dense, sub);
    auto sparse_subbed = bioc::DelayedSubsetOp<double, 1>(sparse, sub);

    first = 0;
    size_t NC = sparse->nrow();
    last = NC;
    set_sizes(NC);

    create_workspaces();
    for (size_t i = 0; i < sub.size(); ++i) {
        wipe_output();
        fill_output(sparse_subbed.get_column(i, output.data()), sub.size());

        wipe_expected();
        fill_expected(sparse->get_column(sub[i], expected.data()), NC);
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_subbed.get_column(i, output.data()), NC);
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_column(i, outval.data(), outidx.data()), NC);
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_subbed.get_sparse_column(i, outval.data(), outidx.data(), work_sparse.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_subbed.get_column(i, outval.data(), work_dense.get()), NC);
        EXPECT_EQ(output, expected);
    }
}


