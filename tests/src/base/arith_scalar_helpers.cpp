#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

class ArithScalarTestCore : public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense, sparse;
protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column major.
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
    tatami::DelayedAddScalarHelper<double> op(5);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] += 5;
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, SubtractionByColumn) {
    tatami::DelayedSubtractScalarHelper<true, double> op(5);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] -= 5;
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }

    // Trying the other side.
    tatami::DelayedSubtractScalarHelper<false> op2(0.7);
    auto dense_mod2 = tatami::make_DelayedIsometricOp(dense, op2);
    auto sparse_mod2 = tatami::make_DelayedIsometricOp(sparse, op2);

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = 0.7 - expected[j];
        }

        wipe_output();
        fill_output(sparse_mod2->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, MultiplicationByColumn) {
    tatami::DelayedMultiplyScalarHelper<double> op(5);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());
    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] *= 5;
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(ArithScalarTest, DivisionByColumn) {
    tatami::DelayedDivideScalarHelper<true> op(5);
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    // Works in full.
    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] /= 5;
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }

    // Works in the other direction. Some shenanigans involved to avoid division by zero.
    tatami::DelayedExpHelper<double> op2a;
    auto dense_mod2a = tatami::make_DelayedIsometricOp(dense, op2a);
    auto sparse_mod2a = tatami::make_DelayedIsometricOp(sparse, op2a);

    tatami::DelayedDivideScalarHelper<false> op2(0.7);
    auto dense_mod2 = tatami::make_DelayedIsometricOp(dense_mod2a, op2);
    auto sparse_mod2 = tatami::make_DelayedIsometricOp(sparse_mod2a, op2);

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = 0.7 / std::exp(expected[j]);
        }

        wipe_output();
        fill_output(sparse_mod2->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);
    }
}

