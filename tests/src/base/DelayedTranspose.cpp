#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "tatami/base/DenseMatrix.hpp"
#include "tatami/base/DelayedTranspose.hpp"
#include "tatami/utils/convert_to_sparse.hpp"

#include "../data/data.h"
#include "TestCore.h"

class TransposeTest: public TestCore {
protected:
    std::shared_ptr<tatami::numeric_matrix> dense;
    std::shared_ptr<tatami::typed_matrix<double, int> > sparse;
    std::shared_ptr<tatami::workspace> work_dense, work_sparse;
protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::shared_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
        sparse = tatami::convert_to_sparse(dense.get(), false); // column-major.
    }

    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

TEST_F(TransposeTest, FullDenseAccess) {
    auto tdense = tatami::make_DelayedTranspose(dense);
    EXPECT_EQ(tdense->ncol(), dense->nrow());
    EXPECT_EQ(tdense->nrow(), dense->ncol());
    EXPECT_EQ(!tdense->prefer_rows(), dense->prefer_rows());

    auto wrk = tdense->new_workspace(false);
    set_sizes(0, tdense->nrow());
    for (size_t i = 0; i < tdense->ncol(); ++i) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data()));

        wipe_output();
        fill_output(tdense->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tdense->column(i, output.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }

    wrk = tdense->new_workspace(true);
    set_sizes(0, tdense->ncol());
    for (size_t i = 0; i < tdense->nrow(); ++i) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data()));

        wipe_output();
        fill_output(tdense->row(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tdense->row(i, expected.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(TransposeTest, SubsetDenseAccess) {
    auto tdense = tatami::make_DelayedTranspose(dense);
    size_t LEN = 6;
    size_t first = 2;

    auto wrk = tdense->new_workspace(false);
    for (size_t i = 0; i < tdense->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, tdense->nrow()));

        wipe_expected();
        fill_expected(dense->row(i, expected.data(), first, last));

        wipe_output();
        fill_output(tdense->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tdense->column(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= tdense->nrow();
    }

    wrk = tdense->new_workspace(true);
    LEN = 7;
    first = 0;
    for (size_t i = 0; i < tdense->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, tdense->ncol()));

        wipe_expected();
        fill_expected(dense->column(i, expected.data(), first, last));

        wipe_output();
        fill_output(tdense->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tdense->row(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= tdense->ncol();
    }
}

TEST_F(TransposeTest, FullSparseAccess) {
    auto tsparse = tatami::make_DelayedTranspose(sparse);
    EXPECT_EQ(tsparse->ncol(), sparse->nrow());
    EXPECT_EQ(tsparse->nrow(), sparse->ncol());
    EXPECT_EQ(!tsparse->prefer_rows(), sparse->prefer_rows());

    auto wrk = tsparse->new_workspace(false);
    set_sizes(0, tsparse->nrow());
    for (size_t i = 0; i < tsparse->ncol(); ++i) {
        wipe_expected();
        fill_expected(sparse->row(i, expected.data()));

        wipe_output();
        fill_output(tsparse->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tsparse->column(i, output.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }

    wrk = tsparse->new_workspace(true);
    set_sizes(0, tsparse->ncol());
    for (size_t i = 0; i < tsparse->nrow(); ++i) {
        wipe_expected();
        fill_expected(sparse->column(i, expected.data()));

        wipe_output();
        fill_output(tsparse->row(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tsparse->row(i, output.data(), wrk.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(TransposeTest, SubsetSparseAccess) {
    auto tsparse = tatami::make_DelayedTranspose(sparse);
    size_t LEN = 3;
    size_t first = 1;

    auto wrk = tsparse->new_workspace(false);
    for (size_t i = 0; i < tsparse->ncol(); ++i) {
        set_sizes(first, std::min(first + LEN, tsparse->nrow()));

        wipe_expected();
        fill_expected(sparse->row(i, expected.data(), first, last));

        wipe_output();
        fill_output(tsparse->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tsparse->column(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= tsparse->nrow();
    }

    LEN = 7;
    first = 0;
    for (size_t i = 0; i < tsparse->nrow(); ++i) {
        set_sizes(first, std::min(first + LEN, tsparse->ncol()));

        wipe_expected();
        fill_expected(sparse->column(i, expected.data(), first, last));

        wipe_output();
        fill_output(tsparse->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(tsparse->row(i, output.data(), first, last, wrk.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= tsparse->ncol();
    }
}
