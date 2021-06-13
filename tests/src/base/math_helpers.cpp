#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedIsometricOp.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

class MathTestCore : public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense, sparse;
protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column major.
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
    tatami::DelayedAbsHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

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
    tatami::DelayedAbsHelper op0;
    auto dense_mod0 = tatami::make_DelayedIsometricOp(dense, op0);
    auto sparse_mod0 = tatami::make_DelayedIsometricOp(sparse, op0);

    tatami::DelayedSqrtHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense_mod0, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse_mod0, op);

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
    tatami::DelayedAbsHelper op0;
    auto dense_mod0 = tatami::make_DelayedIsometricOp(dense, op0);
    auto sparse_mod0 = tatami::make_DelayedIsometricOp(sparse, op0);

    tatami::DelayedAddScalarHelper<double> op1(5);
    auto dense_mod1 = tatami::make_DelayedIsometricOp(dense_mod0, op1);
    auto sparse_mod1 = tatami::make_DelayedIsometricOp(sparse_mod0, op1);

    tatami::DelayedLogHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense_mod1, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse_mod1, op);

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
    tatami::DelayedLogHelper op2(2);
    auto dense_mod2 = tatami::make_DelayedIsometricOp(dense_mod1, op2);
    auto sparse_mod2 = tatami::make_DelayedIsometricOp(sparse_mod1, op2);

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
    tatami::DelayedAbsHelper op0;
    auto dense_mod0 = tatami::make_DelayedIsometricOp(dense, op0);
    auto sparse_mod0 = tatami::make_DelayedIsometricOp(sparse, op0);

    tatami::DelayedLog1pHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense_mod0, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse_mod0, op);

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
    tatami::DelayedLog1pHelper op2(2);
    auto dense_mod2 = tatami::make_DelayedIsometricOp(dense_mod0, op2);
    auto sparse_mod2 = tatami::make_DelayedIsometricOp(sparse_mod0, op2);

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
    tatami::DelayedExpHelper op;
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
    tatami::DelayedRoundHelper op;
    auto dense_mod = tatami::make_DelayedIsometricOp(dense, op);
    auto sparse_mod = tatami::make_DelayedIsometricOp(sparse, op);

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


