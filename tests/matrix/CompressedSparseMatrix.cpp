#include <gtest/gtest.h>

#include "matrix/CompressedSparseMatrix.hpp"
#include "matrix/DenseMatrix.hpp"

#include "data.h"
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

class SparseTestCore : public ::testing::Test {
protected:
    void assemble(size_t nr, size_t nc, const std::vector<double>& source) {
        // Pretend this is row major for the time being.
        dense = std::unique_ptr<bioc::numeric_matrix>(new bioc::DenseRowMatrix<double>(nr, nc, source));

        // Filling the sparse row matrix.
        std::vector<double> values;
        std::vector<int> indices;
        std::vector<size_t> indptr;

        auto sIt = source.begin();
        indptr.resize(nr+1);
        for (size_t i = 0; i < nr; ++i) {
            indptr[i+1] = indptr[i];
            for (size_t j = 0; j < nc; ++j, ++sIt) {
                if (*sIt != 0) {
                    values.push_back(*sIt);
                    indices.push_back(j);
                    ++indptr[i+1];
                }
            }
        }
        sparse_row = std::unique_ptr<bioc::typed_matrix<double, int> >(new bioc::CompressedSparseRowMatrix<double, int>(nr, nc, values, indices, indptr));

        // Filling the sparse column matrix.
        values.clear();
        indices.clear();
        indptr.clear();

        indptr.resize(nc+1);
        for (size_t i = 0; i < nc; ++i) {
            indptr[i+1] = indptr[i];
            auto sIt = source.begin() + i;
            for (size_t j = 0; j < nr; ++j, sIt+=nc) {
                if (*sIt != 0) {
                    values.push_back(*sIt);
                    indices.push_back(j);
                    ++indptr[i+1];
                }
            }
        }
        sparse_column = std::unique_ptr<bioc::typed_matrix<double, int> >(new bioc::CompressedSparseColumnMatrix<double, int>(nr, nc, values, indices, indptr));

        NR = sparse_column->nrow();
        NC = sparse_column->ncol();

        return;
    }

    std::unique_ptr<bioc::numeric_matrix> dense;
    std::unique_ptr<bioc::typed_matrix<double, int> > sparse_row, sparse_column;
    std::unique_ptr<bioc::workspace> work_dense, work_sparse_row, work_sparse_column;

    void create_workspaces() {
        work_dense.reset(dense->create_workspace());
        work_sparse_column.reset(sparse_column->create_workspace());
        work_sparse_row.reset(sparse_row->create_workspace());
        return;
    }

    std::vector<double> expected;
    void fill_expected(const double* ptr, const size_t dim) {
        if (ptr!=expected.data()){ 
            std::copy(ptr, ptr + std::min(dim, last) - first, expected.begin());
        }
        return;
    }

    std::vector<double> output;
    void fill_output(const double* ptr, const size_t dim) {
        if (ptr!=output.data()) {
            std::copy(ptr, ptr + std::min(dim, last) - first , output.begin());
        }
        return;
    }

    void set_sizes(size_t n) {
        output.resize(n);
        outidx.resize(n);
        outval.resize(n);
        expected.resize(n);
        return;
    }

    void wipe_output() {
        std::fill(output.begin(), output.end(), 123);
        return;
    }

    void wipe_expected() {
        std::fill(expected.begin(), expected.end(), 123);
        return;
    }

    std::vector<int> outidx;
    std::vector<double> outval;

    void wipe_sparse_buffers() {
        std::fill(outval.begin(), outval.end(), 123);
        std::fill(outidx.begin(), outidx.end(), 456);
        return;
    }

    void fill_sparse_output(const bioc::sparse_range<double, int>& info, size_t dim) {
        std::fill(output.begin(), output.begin() + std::min(dim, last) - first, 0);
        for (size_t i = 0; i < info.number; ++i) {
            output[info.index[i] - first] = info.value[i];
        }
        return;
    }

    size_t NR, NC;
    size_t first, last;
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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

    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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

    create_workspaces();
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
    create_workspaces();
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
    create_workspaces();
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
