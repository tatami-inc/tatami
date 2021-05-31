#include <gtest/gtest.h>

#include "tatami/CompressedSparseMatrix.hpp"
#include "tatami/DenseMatrix.hpp"

#include "load_sparse.h"
#include "data.h"
#include "TestCore.h"
#include <vector>
#include <memory>

TEST(CompressedSparseMatrix, ConstructionEmpty) {
    std::vector<double> values;
    std::vector<int> indices;
    std::vector<size_t> indptr(21);

    tatami::CompressedSparseColumnMatrix<double, int> mat(10, 20, values, indices, indptr);
    EXPECT_TRUE(mat.is_sparse());
    EXPECT_EQ(mat.nrow(), 10);
    EXPECT_EQ(mat.ncol(), 20);
    EXPECT_EQ(mat.type(), tatami::_double);
}

class SparseTestCore : public TestCore {
protected:
    size_t NR, NC;
    std::unique_ptr<tatami::numeric_matrix> dense;
    std::unique_ptr<tatami::typed_matrix<double, int> > sparse_row, sparse_column;
    std::unique_ptr<tatami::workspace> work_dense, work_sparse_row, work_sparse_column;

protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        dense = std::unique_ptr<tatami::numeric_matrix>(new tatami::DenseRowMatrix<double>(nr, nc, source));
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
    set_sizes(0, NR);
    EXPECT_EQ(NC, sparse_ncol);
    EXPECT_EQ(NR, sparse_nrow);

    EXPECT_FALSE(dense->is_sparse());
    EXPECT_TRUE(sparse_column->is_sparse());
    EXPECT_TRUE(sparse_row->is_sparse());

    EXPECT_EQ(sparse_column->preferred_dimension(), 1);
    EXPECT_EQ(sparse_row->preferred_dimension(), 0);

    // Column access without workspaces.
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data()));

        wipe_output();
        fill_output(sparse_column->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        auto x = sparse_column->sparse_column(i, outval.data(), outidx.data());
        fill_sparse_output(x);
        EXPECT_TRUE(x.number < NR);
        EXPECT_EQ(output, expected);
        EXPECT_FALSE(outval.data()==x.value); // points to internal data.
        EXPECT_FALSE(outidx.data()==x.index);

        wipe_sparse_buffers();
        auto y = sparse_row->sparse_column(i, outval.data(), outidx.data());
        fill_sparse_output(y);
        EXPECT_TRUE(y.number < NR);
        EXPECT_EQ(output, expected);
        EXPECT_TRUE(outval.data()==y.value); // points to buffer.
        EXPECT_TRUE(outidx.data()==y.index);
    }

    // Column access with workspaces.
    create_workspaces(false);
    std::vector<size_t> old_offset_row = dynamic_cast<tatami::CompressedSparseRowMatrix<double, int>::compressed_sparse_workspace*>(work_sparse_row.get())->offsets();
    EXPECT_EQ(work_sparse_column.get(), nullptr);

    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    std::vector<size_t> new_offset_row = dynamic_cast<tatami::CompressedSparseRowMatrix<double, int>::compressed_sparse_workspace*>(work_sparse_row.get())->offsets();
    EXPECT_NE(old_offset_row, new_offset_row); // actually has an effect on the offsets.

    // Column access with workspaces and reverse order.
    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->column(NC - i - 1, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(NC - i - 1, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(NC - i - 1, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    // Column access with workspaces and jump forward.
    create_workspaces(false);
    for (size_t i = 0; i < NC; i+=2) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    // Column access with workspaces and jump backward.
    create_workspaces(false);
    for (size_t i = 0; i < NC; i+=3) {
        wipe_expected();
        fill_expected(dense->column(NC - i - 1, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(NC - i - 1, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(NC - i - 1, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(NC - i - 1, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SparseTest, SlicedColumnAccess) {
    first = NR / 5;
    last = NR / 2;
    set_sizes(first, last);

    // Constant slicing, with and without workspaces.
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data(), first, last));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);
    }

    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        wipe_expected();
        fill_expected(dense->column(i, expected.data(), first, last, work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    // Variable restriction, with and without workspaces.
    first = 0;
    size_t LEN = 5;
    for (size_t i = 0; i < NC; ++i) {
        set_sizes(first, std::min(first + LEN, NR));

        wipe_expected();
        fill_expected(dense->column(i, expected.data(), first, last));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NR;
    }

    first = 0;
    create_workspaces(false);
    for (size_t i = 0; i < NC; ++i) {
        set_sizes(first, std::min(first + LEN, NR));

        wipe_expected();
        fill_expected(dense->column(i, expected.data(), first, last, work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NR;
    }

    // Variable restriction, with workspaces and jumps.
    first = 0;
    create_workspaces(false);
    for (size_t i = 0; i < NC; i+=3) {
        set_sizes(first, std::min(first + LEN, NR));

        wipe_expected();
        fill_expected(dense->column(i, expected.data(), first, last, work_dense.get()));

        wipe_output();
        fill_output(sparse_column->column(i, output.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->column(i, output.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_column(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        first += 7;
        first %= NR;
    }
}

TEST_F(SparseTest, FullRowAccess) {
    set_sizes(0, NC);

    // Row access without workspaces.
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data()));

        wipe_output();
        fill_output(sparse_column->row(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(i, output.data()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        auto x = sparse_column->sparse_row(i, outval.data(), outidx.data());
        fill_sparse_output(x);
        EXPECT_TRUE(x.number < NC);
        EXPECT_EQ(output, expected);
        EXPECT_TRUE(outval.data()==x.value); // points to buffer.
        EXPECT_TRUE(outidx.data()==x.index);

        wipe_sparse_buffers();
        auto y = sparse_row->sparse_row(i, outval.data(), outidx.data());
        fill_sparse_output(y);
        EXPECT_TRUE(y.number < NC);
        EXPECT_EQ(output, expected);
        EXPECT_FALSE(outval.data()==y.value); // points to internal data.
        EXPECT_FALSE(outidx.data()==y.index);
    }

    // Row access with workspaces.
    create_workspaces(true);
    std::vector<size_t> old_offset_column = dynamic_cast<tatami::CompressedSparseColumnMatrix<double, int>::compressed_sparse_workspace*>(work_sparse_column.get())->offsets();
    EXPECT_EQ(work_sparse_row.get(), nullptr);

    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(i, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    std::vector<size_t> new_offset_column = dynamic_cast<tatami::CompressedSparseColumnMatrix<double, int>::compressed_sparse_workspace*>(work_sparse_column.get())->offsets();
    EXPECT_NE(old_offset_column, new_offset_column); // actually has an effect on the offsets.

    // Row access with workspaces and reverse order.
    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->row(NR - i - 1, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(NR - i - 1, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(NR - i - 1, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    // Row access with workspaces and jump forward.
    create_workspaces(true);
    for (size_t i = 0; i < NR; i+=2) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output(); 
        fill_output(sparse_row->row(i, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    // Row access with workspaces and jump backward.
    create_workspaces(true);
    for (size_t i = 0; i < NR; i+=3) {
        wipe_expected();
        fill_expected(dense->row(NR - i - 1, expected.data(), work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(NR - i - 1, output.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(NR - i - 1, output.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(NR - i - 1, outval.data(), outidx.data(), work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }
}

TEST_F(SparseTest, SlicedRowAccess) {
    first = NC / 5;
    last = NC / 2;
    set_sizes(first, last);

    // Constant slicing, with and without workspaces.
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data(), first, last));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);
    }

    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        wipe_expected();
        fill_expected(dense->row(i, expected.data(), first, last, work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(i, output.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);
    }

    // Variable restriction, with and without workspaces.
    size_t LEN = 5;
    first = 0;
    for (size_t i = 0; i < NR; ++i) {
        set_sizes(first, std::min(NC, first + LEN));

        wipe_expected();
        fill_expected(dense->row(i, expected.data(), first, last));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(i, output.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), first, last));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NC;
    }

    first = 0;
    create_workspaces(true);
    for (size_t i = 0; i < NR; ++i) {
        set_sizes(first, std::min(NC, first + LEN));

        wipe_expected();
        fill_expected(dense->row(i, expected.data(), first, last, work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output();
        fill_output(sparse_row->row(i, output.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        first += 3;
        first %= NC;
    }

    // Variable restriction, with workspaces and jumps.
    LEN = 10;
    first = 0;
    create_workspaces(true);
    for (size_t i = 0; i < NR; i+=3) {
        set_sizes(first, std::min(NC, first + LEN));

        wipe_expected();
        fill_expected(dense->row(i, expected.data(), first, last, work_dense.get()));

        wipe_output();
        fill_output(sparse_column->row(i, output.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_output(); 
        fill_output(sparse_row->row(i, output.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_column->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_column.get()));
        EXPECT_EQ(output, expected);

        wipe_sparse_buffers();
        fill_sparse_output(sparse_row->sparse_row(i, outval.data(), outidx.data(), first, last, work_sparse_row.get()));
        EXPECT_EQ(output, expected);

        first += 7;
        first %= NC;
    }
}
