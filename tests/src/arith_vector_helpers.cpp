#include <gtest/gtest.h>

#include "tatami/CompressedSparseMatrix.hpp"
#include "tatami/DenseMatrix.hpp"
#include "tatami/DelayedIsometricOp.hpp"

#include "data.h"
#include "load_sparse.h"
#include "TestCore.h"

#include <vector>
#include <memory>

class ArithVectorTestCore : public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense;
    std::shared_ptr<tatami::typed_matrix<double, int> > sparse;
    std::unique_ptr<tatami::workspace> work_dense, work_sparse;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = load_matrix_as_sparse_column_matrix(nr, nc, source);
        return;
    }

    void create_workspaces(bool row) {
        work_dense.reset(dense->new_workspace(row));
        work_sparse.reset(sparse->new_workspace(row));
        return;
    }

    std::vector<double> create_vector(size_t n, double starter){ 
        std::vector<double> output(n, starter);
        for (size_t i = 1; i < n; ++i) {
            output[i] += i;
        }
        return output;
    }
};

class ArithVectorTest : public ArithVectorTestCore {
protected:
    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

/****************************
 ********* ADDITION *********
 ****************************/

TEST_F(ArithVectorTest, AdditionAlongRows) {
    auto vec = create_vector(dense->nrow(), 0.2);
    tatami::DelayedAddVectorHelper<double> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.sparse());
    EXPECT_FALSE(sparse_mod.sparse());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    EXPECT_TRUE(dense->prefer_rows());
    EXPECT_FALSE(sparse->prefer_rows());

    // Works in full.
    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] += vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->nrow()/6, dense->nrow()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] += vec[first + j];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Row extraction works as well.
    create_workspaces(true);
    set_sizes(0, dense->ncol());

    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] += vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithVectorTest, AdditionAlongColumns) {
    auto vec = create_vector(dense->ncol(), -10.7);
    tatami::DelayedAddVectorHelper<double, 1> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(true);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] += vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->ncol()/6, dense->ncol()/2);

    create_workspaces(true);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] += vec[first + j];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Column extraction works as well.
    create_workspaces(false);
    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] += vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

/****************************
 ******** SUBTRACTION *******
 ****************************/

TEST_F(ArithVectorTest, SubtractionAlongRows) {
    auto vec = create_vector(dense->nrow(), -1.2);
    tatami::DelayedSubtractVectorHelper<double, true> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.sparse());
    EXPECT_FALSE(sparse_mod.sparse());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    // Works in full.
    set_sizes(0, dense->nrow());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] -= vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->nrow()/6, dense->nrow()/2);
    tatami::DelayedSubtractVectorHelper<double, false> op2(vec);
    auto dense_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(dense, op2);
    auto sparse_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(sparse, op2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = vec[first + j] - expected[j];
        }

        wipe_output();
        fill_output(sparse_mod2.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Row extraction works as well.
    create_workspaces(true);
    set_sizes(0, dense->ncol());

    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] -= vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithVectorTest, SubtractionAlongColumns) {
    auto vec = create_vector(dense->ncol(), 24.1);
    tatami::DelayedSubtractVectorHelper<double, true, 1> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] -= vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    tatami::DelayedSubtractVectorHelper<double, false, 1> op2(vec);
    auto dense_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(dense, op2);
    auto sparse_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(sparse, op2);

    set_sizes(dense->ncol()/6, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = vec[first + j] - expected[j];
        }

        wipe_output();
        fill_output(sparse_mod2.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Column extraction works as well.
    create_workspaces(false);
    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] -= vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

/*******************************
 ******* MULTIPLICATION ********
 *******************************/

TEST_F(ArithVectorTest, MultiplicationAlongRows) {
    auto vec = create_vector(dense->nrow(), 0.1);
    tatami::DelayedMultiplyVectorHelper<double> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.sparse());
    EXPECT_TRUE(sparse_mod.sparse());

    // Works in full.
    set_sizes(0, dense->nrow());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] *= vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->nrow()/6, dense->nrow()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] *= vec[first + j];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Row extraction works as well.
    create_workspaces(true);
    set_sizes(0, dense->ncol());

    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] *= vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithVectorTest, MultiplicationAlongColumns) {
    auto vec = create_vector(dense->ncol(), -10.5);
    tatami::DelayedMultiplyVectorHelper<double, 1> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] *= vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice.
    set_sizes(dense->ncol()/6, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] *= vec[j + first];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Column extraction works as well.
    create_workspaces(false);
    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] *= vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

/*************************
 ******** DIVISION *******
 *************************/

TEST_F(ArithVectorTest, DivisionAlongRows) {
    auto vec = create_vector(dense->nrow(), 1.2);
    tatami::DelayedDivideVectorHelper<double, true> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    EXPECT_FALSE(dense_mod.sparse());
    EXPECT_TRUE(sparse->sparse());
    EXPECT_TRUE(sparse_mod.sparse());

    // Works in full.
    set_sizes(0, dense->nrow());
    EXPECT_EQ(dense->nrow(), dense_mod.nrow());
    EXPECT_EQ(dense->ncol(), dense_mod.ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] /= vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice. Some shenanigans involved to avoid division by zero.
    tatami::DelayedExpHelper<double> op2a;
    auto dense_mod2a = std::shared_ptr<tatami::typed_matrix<double> >(new tatami::DelayedIsometricOp<double, decltype(op2a)>(dense, op2a));
    auto sparse_mod2a = std::shared_ptr<tatami::typed_matrix<double> >(new tatami::DelayedIsometricOp<double, decltype(op2a)>(sparse, op2a));

    tatami::DelayedDivideVectorHelper<double, false> op2(vec);
    auto dense_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(dense_mod2a, op2);
    auto sparse_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(sparse_mod2a, op2);

    set_sizes(dense->nrow()/10, dense->nrow()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = vec[j + first] / std::exp(expected[j]);
        }

        wipe_output();
        fill_output(sparse_mod2.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.column(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Row extraction works as well.
    create_workspaces(true);
    set_sizes(0, dense->ncol());

    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] /= vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithVectorTest, DivisionAlongColumns) {
    auto vec = create_vector(dense->ncol(), -10.2);
    tatami::DelayedDivideVectorHelper<double, true, 1> op(vec);
    auto dense_mod = tatami::DelayedIsometricOp<double, decltype(op)>(dense, op);
    auto sparse_mod = tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op);

    // Works in full.
    set_sizes(0, dense->ncol());

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));
        for (size_t j = 0; j < dense->ncol(); ++j) {
            expected[j] /= vec[j];
        }

        wipe_output();
        fill_output(sparse_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.row(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_row(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.row(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }

    // Works in a slice. Some shenanigans involved to avoid division by zero.
    tatami::DelayedExpHelper<double> op2a;
    auto dense_mod2a = std::shared_ptr<tatami::typed_matrix<double> >(new tatami::DelayedIsometricOp<double, decltype(op2a)>(dense, op2a));
    auto sparse_mod2a = std::shared_ptr<tatami::typed_matrix<double> >(new tatami::DelayedIsometricOp<double, decltype(op2a)>(sparse, op2a));

    tatami::DelayedDivideVectorHelper<double, false, 1> op2(vec);
    auto dense_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(dense_mod2a, op2);
    auto sparse_mod2 = tatami::DelayedIsometricOp<double, decltype(op2)>(sparse_mod2a, op2);

    set_sizes(dense->ncol()/4, dense->ncol()/2);

    create_workspaces(false);
    for (size_t i = 0; i < dense->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data(), first, last));
        for (size_t j = 0; j < last - first; ++j) {
            expected[j] = vec[j + first] / std::exp(expected[j]);
        }

        wipe_output();
        fill_output(sparse_mod2.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod2.row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod2.sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2.row(i, outval.data(), first, last, work_dense.get()));
        EXPECT_EQ(output, expected);
    }
    
    // Column extraction works as well.
    create_workspaces(false);
    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] /= vec[i];
        }

        wipe_output();
        fill_output(sparse_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);
        
        // Same result regardless of the backend.
        wipe_output();
        fill_output(dense_mod.column(i, output.data()));
        EXPECT_EQ(output, expected);

        // Works in sparse mode as well.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data()));
        EXPECT_EQ(output, expected);

        // Passes along the workspace.
        wipe_output();
        fill_sparse_output(sparse_mod.sparse_column(i, outval.data(), outidx.data(), work_sparse.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod.column(i, outval.data(), work_dense.get()));
        EXPECT_EQ(output, expected);
    }
}

