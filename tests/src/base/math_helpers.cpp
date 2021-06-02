#include <gtest/gtest.h>

#include "tatami/tatami.h"
#include "../data.h"
#include "../load_sparse.h"
#include "TestCore.h"

#include <vector>
#include <memory>

class MathTestCore : public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense;
    std::shared_ptr<tatami::typed_matrix<double, int> > sparse;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = load_matrix_as_sparse_column_matrix(nr, nc, source);
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
    tatami::DelayedAbsHelper<double> op;
    auto dense_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(dense, op));
    auto sparse_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::abs(expected[j]); 
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, SqrtByColumn) {
    tatami::DelayedAbsHelper<double> op0;
    auto dense_mod0 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op0)>(dense, op0));
    auto sparse_mod0 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op0)>(sparse, op0));

    tatami::DelayedSqrtHelper<double> op;
    auto dense_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(dense_mod0, op));
    auto sparse_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(sparse_mod0, op));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::sqrt(std::abs(expected[j])); 
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, LogByColumn) {
    tatami::DelayedAbsHelper<double> op0;
    auto dense_mod0 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op0)>(dense, op0));
    auto sparse_mod0 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op0)>(sparse, op0));

    tatami::DelayedAddScalarHelper<double> op1(5);
    auto dense_mod1 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op1)>(dense_mod0, op1));
    auto sparse_mod1 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op1)>(sparse_mod0, op1));

    tatami::DelayedLogHelper<double> op;
    auto dense_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(dense_mod1, op));
    auto sparse_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(sparse_mod1, op));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log(std::abs(expected[j]) + 5);
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }

    // Trying with another base.
    tatami::DelayedLogHelper<double> op2(2);
    auto dense_mod2 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op2)>(dense_mod1, op2));
    auto sparse_mod2 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op2)>(sparse_mod1, op2));

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log(std::abs(expected[j]) + 5)/std::log(2);
        }

        wipe_output();
        fill_output(sparse_mod2->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, Log1pByColumn) {
    tatami::DelayedAbsHelper<double> op0;
    auto dense_mod0 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op0)>(dense, op0));
    auto sparse_mod0 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op0)>(sparse, op0));

    tatami::DelayedLog1pHelper<double> op;
    auto dense_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(dense_mod0, op));
    auto sparse_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(sparse_mod0, op));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log1p(std::abs(expected[j]));
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }

    // Trying with another base.
    tatami::DelayedLog1pHelper<double> op2(2);
    auto dense_mod2 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op2)>(dense_mod0, op2));
    auto sparse_mod2 = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op2)>(sparse_mod0, op2));

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::log1p(std::abs(expected[j]))/std::log(2);
        }

        wipe_output();
        fill_output(sparse_mod2->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod2->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, ExpByColumn) {
    tatami::DelayedExpHelper<double> op;
    auto dense_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(dense, op));
    auto sparse_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_FALSE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::exp(expected[j]); 
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(MathTest, RoundByColumn) {
    tatami::DelayedRoundHelper<double> op;
    auto dense_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(dense, op));
    auto sparse_mod = std::shared_ptr<tatami::numeric_matrix>(new tatami::DelayedIsometricOp<double, decltype(op)>(sparse, op));

    EXPECT_FALSE(dense_mod->sparse());
    EXPECT_TRUE(sparse_mod->sparse());
    EXPECT_EQ(dense->nrow(), dense_mod->nrow());
    EXPECT_EQ(dense->ncol(), dense_mod->ncol());

    set_sizes(0, dense->nrow());

    for (size_t i = 0; i < dense->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));
        for (size_t j = 0; j < dense->nrow(); ++j) {
            expected[j] = std::round(expected[j]); 
        }

        wipe_output();
        fill_output(sparse_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(dense_mod->column(i, output.data()));
        EXPECT_EQ(output, expected);
    }
}


