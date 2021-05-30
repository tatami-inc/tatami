#include <gtest/gtest.h>

#include "matrix/CompressedSparseMatrix.hpp"
#include "matrix/DenseMatrix.hpp"

#include "../utils/load_sparse.h"
#include "../utils/data.h"
#include "../utils/TestCore.h"
#include <vector>
#include <memory>

TEST(CompressedSparseMatrix, ConstructionEmpty) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    bioc::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr);
    EXPECT_TRUE(mat.is_sparse());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
    EXPECT_EQ(mat.type(), bioc::_double);
}

class SparseTestCore : public TestCore {
protected:
    size_t NR, NC;
    std::unique_ptr<bioc::numeric_matrix> dense;
    std::unique_ptr<bioc::typed_matrix<double, int> > sparse_row, sparse_column;
    std::unique_ptr<bioc::workspace> work_dense, work_sparse_row, work_sparse_column;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::unique_ptr<bioc::numeric_matrix>(new bioc::DenseRowMatrix<double>(nr, nc, source));
        sparse_row = load_matrix_as_sparse_row_matrix(nr, nc, source);
        sparse_column = load_matrix_as_sparse_column_matrix(nr, nc, source);

        NR = sparse_column->nrow();
        NC = sparse_column->ncol();
        return;
    }

    void create_workspaces(bool row) {
        work_dense.reset(dense->create_workspace(row));
        work_sparse_column.reset(sparse_column->create_workspace(row));
        work_sparse_row.reset(sparse_row->create_workspace(row));
        return;
    }
};

class SparseTest : public SparseTestCore {
protected:
    void SetUp() {
        assemble(sparse_nrow, sparse_ncol, sparse_matrix);
    }
};

TEST_F(SparseTest, FullColumnAccess) {
    first = 0;
    last = NR;
    EXPECT_EQ(NC, sparse_ncol);
    EXPECT_EQ(NR, sparse_nrow);
    set_sizes(NR);

    // Column access without workspaces.
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->get_column(i, expected.data()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        auto x = sparse_column->get_sparse_column(i, outval.data(), outidx.data());
        fill_sparse_output(x, NR);
        EXPECT_EQ(output, expected);
        EXPECT_FALSE(outval.data()==x.value); // points to internal data.
        EXPECT_FALSE(outidx.data()==x.index);

        wipe_sparse_buffers();
        auto y = sparse_row->get_sparse_column(i, outval.data(), outidx.data());
        fill_sparse_output(y, NR);
        EXPECT_EQ(output, expected);
        EXPECT_TRUE(outval.data()==y.value); // points to buffer.
        EXPECT_TRUE(outidx.data()==y.index);
    }

    // Column access with workspaces.
    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);
    }

    // Column access with workspaces and reverse order.
    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->get_column(NC - i - 1, expected.data(), work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(NC - i - 1, output.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(NC - i - 1, output.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);
    }

    // Column access with workspaces and jump forward.
    create_workspaces(false);
    for (size_t i = 0; i < NC; i+=2) {
        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);
    }

    // Column access with workspaces and jump backward.
    create_workspaces(false);
    for (size_t i = 0; i < NC; i+=3) {
        wipe_expected();
        fill_expected(dense->get_column(NC - i - 1, expected.data(), work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(NC - i - 1, output.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(NC - i - 1, output.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SparseTest, SlicedColumnAccess) {
    first = NR / 5;
    last = NR / 2;
    set_sizes(last - first);

    // Constant slicing, with and without workspaces.
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), first, last), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), first, last), NR);
        EXPECT_EQ(output, expected);
    }

    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), first, last, work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), first, last, work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), first, last, work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);
    }

    // Variable restriction, with and without workspaces.
    size_t LEN = 5;
    set_sizes(LEN);

    first = 0;
    for (size_t i = 0; i < NC; ++i) {
        last = first + LEN;

        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), first, last), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), first, last), NR);
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NR;
    }

    first = 0;
    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        last = first + LEN;

        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), first, last, work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), first, last, work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), first, last, work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NR;
    }

    // Variable restriction, with workspaces and jumps.
    LEN = 10;
    set_sizes(LEN);

    first = 0;
    create_workspaces(false);
    for (size_t i = 0; i < NC; i+=3) {
        last = first + LEN;

        wipe_expected();
        fill_expected(dense->get_column(i, expected.data(), first, last, work_dense.get()), NR);

        wipe_output();
        fill_output(sparse_column->get_column(i, output.data(), first, last, work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_column(i, output.data(), first, last, work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()), NR);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()), NR);
        EXPECT_EQ(output, expected);

        first += 7;
        first %= NR;
    }
}

TEST_F(SparseTest, FullRowAccess) {
    first = 0;
    last = NC;
    set_sizes(NC);

    // Row access without workspaces.
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->get_row(i, expected.data()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(i, output.data()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        auto x = sparse_column->get_sparse_row(i, outval.data(), outidx.data());
        fill_sparse_output(x, NC);
        EXPECT_EQ(output, expected);
        EXPECT_TRUE(outval.data()==x.value); // points to buffer.
        EXPECT_TRUE(outidx.data()==x.index);

        wipe_sparse_buffers();
        auto y = sparse_row->get_sparse_row(i, outval.data(), outidx.data());
        fill_sparse_output(y, NC);
        EXPECT_EQ(output, expected);
        EXPECT_FALSE(outval.data()==y.value); // points to internal data.
        EXPECT_FALSE(outidx.data()==y.index);
    }

    // Row access with workspaces.
    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(i, output.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);
    }

    // Row access with workspaces and reverse order.
    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->get_row(NR - i - 1, expected.data(), work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(NR - i - 1, output.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(NR - i - 1, output.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);
    }

    // Row access with workspaces and jump forward.
    create_workspaces(true);
    for (size_t i = 0; i < NR; i+=2) {
        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output(); 
        fill_output(sparse_row->get_row(i, output.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);
    }

    // Row access with workspaces and jump backward.
    create_workspaces(true);
    for (size_t i = 0; i < NR; i+=3) {
        wipe_expected();
        fill_expected(dense->get_row(NR - i - 1, expected.data(), work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(NR - i - 1, output.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(NR - i - 1, output.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SparseTest, SlicedRowAccess) {
    first = NC / 5;
    last = NC / 2;
    set_sizes(last - first);

    // Constant slicing, with and without workspaces.
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), first, last), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(i, output.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), first, last), NC);
        EXPECT_EQ(output, expected);
    }

    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), first, last, work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), first, last, work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(i, output.data(), first, last, work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);
    }

    // Variable restriction, with and without workspaces.
    size_t LEN = 5;
    set_sizes(LEN);

    first = 0;
    for (size_t i = 0; i < NR; ++i) {
        last = first + LEN;

        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), first, last), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(i, output.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), first, last), NC);
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NC;
    }

    first = 0;
    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        last = first + LEN;

        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), first, last, work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), first, last, work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->get_row(i, output.data(), first, last, work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NC;
    }

    // Variable restriction, with workspaces and jumps.
    LEN = 10;
    set_sizes(LEN);

    first = 0;
    create_workspaces(true);
    for (size_t i = 0; i < NR; i+=3) {
        last = first + LEN;

        wipe_expected();
        fill_expected(dense->get_row(i, expected.data(), first, last, work_dense.get()), NC);

        wipe_output();
        fill_output(sparse_column->get_row(i, output.data(), first, last, work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_output(); 
        fill_output(sparse_row->get_row(i, output.data(), first, last, work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()), NC);
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->get_sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()), NC);
        EXPECT_EQ(output, expected);

        first += 7;
        first %= NC;
    }
}
