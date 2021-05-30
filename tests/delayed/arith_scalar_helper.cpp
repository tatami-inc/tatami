#include <gtest/gtest.h>

#include "matrix/CompressedSparseMatrix.hpp"
#include "matrix/DenseMatrix.hpp"
#include "delayed/DelayedIsometricOp.hpp"

#include "../utils/data.h"
#include "../utils/load_sparse.h"
#include "../utils/TestCore.h"

#include <vector>
#include <memory>

class ArithScalarTestCore : public TestCore {
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

class ArithScalarTest : public ArithScalarTestCore {
protected:
    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

/****************************
 ********* ADDITION *********
 ****************************/

TEST_F(ArithScalarTest, AdditionByColumn) {
    bioc::DelayedAddScalarHelper<double> op(5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    // Works in full.
    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] += 5;
        }

        wipe_output();
        fill_output(sparse_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->nrow()/6, dense->nrow()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] += 5;
        }

        wipe_output();
        fill_output(sparse_mod.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, AdditionByRow) {
    bioc::DelayedAddScalarHelper<double> op(0.5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] += 0.5;
        }

        wipe_output();
        fill_output(sparse_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->ncol()/6, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] += 0.5;
        }

        wipe_output();
        fill_output(sparse_mod.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

/****************************
 ******** SUBTRACTION *******
 ****************************/

TEST_F(ArithScalarTest, SubtractionByColumn) {
    bioc::DelayedSubtractScalarHelper<double, true> op(5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_FALSE(sparse_mod.is_sparse());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    // Works in full.
    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] -= 5;
        }

        wipe_output();
        fill_output(sparse_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->nrow()/6, dense->nrow()/2);
    bioc::DelayedSubtractScalarHelper<double, false> op2(0.7);
    auto dense_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(dense, op2);
    auto sparse_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(sparse, op2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = 0.7 - expected[j];
        }

        wipe_output();
        fill_output(sparse_mod2.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.get_column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, SubtractionByRow) {
    bioc::DelayedSubtractScalarHelper<double, true> op(0.5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] -= 0.5;
        }

        wipe_output();
        fill_output(sparse_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    bioc::DelayedSubtractScalarHelper<double, false> op2(0.1);
    auto dense_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(dense, op2);
    auto sparse_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(sparse, op2);

    set_sizes(dense->ncol()/6, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = 0.1 - expected[j];
        }

        wipe_output();
        fill_output(sparse_mod2.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.get_row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

/*******************************
 ******* MULTIPLICATION ********
 *******************************/

TEST_F(ArithScalarTest, MultiplicationByColumn) {
    bioc::DelayedMultiplyScalarHelper<double> op(5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Works in full.
    set_sizes(0, dense->nrow());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] *= 5;
        }

        wipe_output();
        fill_output(sparse_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->nrow()/6, dense->nrow()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] *= 5;
        }

        wipe_output();
        fill_output(sparse_mod.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, MultiplicationByRow) {
    bioc::DelayedMultiplyScalarHelper<double> op(0.5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] *= 0.5;
        }

        wipe_output();
        fill_output(sparse_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->ncol()/6, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] *= 0.5;
        }

        wipe_output();
        fill_output(sparse_mod.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

/*************************
 ******** DIVISION *******
 *************************/

TEST_F(ArithScalarTest, DivisionByColumn) {
    bioc::DelayedDivideScalarHelper<double, true> op(5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.is_sparse());
    EXPECT_TRUE(sparse->is_sparse());
    EXPECT_TRUE(sparse_mod.is_sparse());

    // Works in full.
    set_sizes(0, dense->nrow());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] /= 5;
        }

        wipe_output();
        fill_output(sparse_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice. Some shenanigans involved to avoid division by zero.
    bioc::DelayedExpHelper<double> op2a;
    auto dense_mod2a = std::shared_ptr<bioc::typed_matrix<double> >(new bioc::DelayedIsometricOp<double, decltype(op2a)>(dense, op2a));
    auto sparse_mod2a = std::shared_ptr<bioc::typed_matrix<double> >(new bioc::DelayedIsometricOp<double, decltype(op2a)>(sparse, op2a));

    bioc::DelayedDivideScalarHelper<double, false> op2(0.7);
    auto dense_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(dense_mod2a, op2);
    auto sparse_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(sparse_mod2a, op2);

    set_sizes(dense->nrow()/10, dense->nrow()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = 0.7 / std::exp(expected[j]);
        }

        wipe_output();
        fill_output(sparse_mod2.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.get_column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.get_column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, DivisionByRow) {
    bioc::DelayedDivideScalarHelper<double, true> op(0.5);
    auto dense_mod = bioc::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = bioc::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] /= 0.5;
        }

        wipe_output();
        fill_output(sparse_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.get_row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.get_sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.get_row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice. Some shenanigans involved to avoid division by zero.
    bioc::DelayedExpHelper<double> op2a;
    auto dense_mod2a = std::shared_ptr<bioc::typed_matrix<double> >(new bioc::DelayedIsometricOp<double, decltype(op2a)>(dense, op2a));
    auto sparse_mod2a = std::shared_ptr<bioc::typed_matrix<double> >(new bioc::DelayedIsometricOp<double, decltype(op2a)>(sparse, op2a));

    bioc::DelayedDivideScalarHelper<double, false> op2(0.1);
    auto dense_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(dense_mod2a, op2);
    auto sparse_mod2 = bioc::DelayedIsometricOp<double, decltype(op2)>(sparse_mod2a, op2);

    set_sizes(dense->ncol()/4, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->get_row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = 0.1 / std::exp(expected[j]);
        }

        wipe_output();
        fill_output(sparse_mod2.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.get_row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.get_row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

