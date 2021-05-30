#include <gtest/gtest.h>

#include "matrix/CompressedSparseMatrix.hpp"
#include "matrix/DenseMatrix.hpp"
#include "delayed/DelayedIsometricOp.hpp"

#include "../utils/data.h"
#include "../utils/load_sparse.h"
#include "../utils/TestCore.h"

#include <vector>
#include <memory>

class MathTestCore : public TestCore {
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

    void create_workspaces(bool row) {
        work_dense.reset(dense->create_workspace(row));
        work_sparse.reset(sparse->create_workspace(row));
        return;
    }
};

class MathTest : public MathTestCore {
protected:
    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

TEST_F(MathTest, AbsByColumn) {
    bioc::DelayedAbsHelper<double> op;
    auto dense_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(dense, op));
    auto sparse_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op));

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::abs(expected[j]); 
        }

        wipe_output();
        fill_output(sparse_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, SqrtByColumn) {
    bioc::DelayedAbsHelper<double> op0;
    auto dense_mod0 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op0)>(dense, op0));
    auto sparse_mod0 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op0)>(sparse, op0));

    bioc::DelayedSqrtHelper<double> op;
    auto dense_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(dense_mod0, op));
    auto sparse_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(sparse_mod0, op));

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::sqrt(std::abs(expected[j])); 
        }

        wipe_output();
        fill_output(sparse_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, LogByColumn) {
    bioc::DelayedAbsHelper<double> op0;
    auto dense_mod0 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op0)>(dense, op0));
    auto sparse_mod0 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op0)>(sparse, op0));

    bioc::DelayedAddScalarHelper<double> op1(5);
    auto dense_mod1 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op1)>(dense_mod0, op1));
    auto sparse_mod1 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op1)>(sparse_mod0, op1));

    bioc::DelayedLogHelper<double> op;
    auto dense_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(dense_mod1, op));
    auto sparse_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(sparse_mod1, op));

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log(std::abs(expected[j]) + 5);
        }

        wipe_output();
        fill_output(sparse_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }

    // Trying with another base.
    bioc::DelayedLogHelper<double> op2(2);
    auto dense_mod2 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op2)>(dense_mod1, op2));
    auto sparse_mod2 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op2)>(sparse_mod1, op2));

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log(std::abs(expected[j]) + 5)/std::log(2);
        }

        wipe_output();
        fill_output(sparse_mod2->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, Log1pByColumn) {
    bioc::DelayedAbsHelper<double> op0;
    auto dense_mod0 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op0)>(dense, op0));
    auto sparse_mod0 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op0)>(sparse, op0));

    bioc::DelayedLog1pHelper<double> op;
    auto dense_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(dense_mod0, op));
    auto sparse_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(sparse_mod0, op));

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log1p(std::abs(expected[j]));
        }

        wipe_output();
        fill_output(sparse_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }

    // Trying with another base.
    bioc::DelayedLog1pHelper<double> op2(2);
    auto dense_mod2 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op2)>(dense_mod0, op2));
    auto sparse_mod2 = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op2)>(sparse_mod0, op2));

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log1p(std::abs(expected[j]))/std::log(2);
        }

        wipe_output();
        fill_output(sparse_mod2->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, ExpByColumn) {
    bioc::DelayedExpHelper<double> op;
    auto dense_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(dense, op));
    auto sparse_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op));

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_FALSE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::exp(expected[j]); 
        }

        wipe_output();
        fill_output(sparse_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, RoundByColumn) {
    bioc::DelayedRoundHelper<double> op;
    auto dense_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(dense, op));
    auto sparse_mod = std::shared_ptr<bioc::numeric_matrix>(new bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op));

    EXPECT_FALSE(dense_mod->is_sparse());
    EXPECT_TRUE(sparse_mod->is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::round(expected[j]); 
        }

        wipe_output();
        fill_output(sparse_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->get_column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

