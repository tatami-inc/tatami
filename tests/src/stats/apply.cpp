#include <gtest/gtest.h>
#include <vector>

#include "tatami/stats/apply.hpp"
#include "tatami/base/DenseMatrix.hpp"
#include "tatami/utils/convert_to_dense.hpp"
#include "tatami/utils/convert_to_sparse.hpp"
#include "tatami/stats/variances.hpp"

#include "../data/data.h"

/* This suite of tests just checks that we actually 
 * get different results under each mode. */

struct SillyFactory {
public:
    SillyFactory(double* o, size_t d1) : output(o), dim(d1) {}

private:
    double* output;
    size_t dim;

public:
    struct DenseDirect {
        DenseDirect(double* o) : output(o) {}

        template<typename V>
        void compute(size_t i, const V* ptr, V* buffer) {
            output[i] = 1;
        }
    private:
        double* output;
    };

    DenseDirect dense_direct() {
        return DenseDirect(output);
    }

public:
    struct SparseDirect {
        SparseDirect(double* o) : output(o) {}

        template<typename T, typename IDX>
        void compute(size_t i, const tatami::SparseRange<T, IDX>& range, T* vbuffer = NULL, IDX* ibuffer = NULL) {
            output[i] = 2;
        }
    private:
        double* output;
    };

    SparseDirect sparse_direct() {
        return SparseDirect(output);
    }

public:
    struct DenseRunning {
        DenseRunning(double* o, size_t d) : output(o), dim(d) {}

        template<typename V>
        void add(const V* ptr, V* buffer = NULL) {
            for (size_t d = 0; d < dim; ++d) {
                output[d] = 3;
            }
        }

        void finish() {}
    private:
        double* output;
        size_t dim;
    };

    DenseRunning dense_running() {
        return DenseRunning(output, dim);
    }

    DenseRunning dense_running(size_t start, size_t end) {
        return DenseRunning(output + start, end - start);
    }

public:
    struct SparseRunning {
        SparseRunning(double* o, size_t d) : output(o), dim(d) {}

        template<typename T = double, typename IDX = int>
        void add(const tatami::SparseRange<T, IDX>& range, T* vbuffer = NULL, IDX* ibuffer = NULL) {
            for (size_t d = 0; d < dim; ++d) {
                output[d] = 4;
            }
        }

        void finish() {}
    private:
        double* output;
        size_t dim;
    };

    SparseRunning sparse_running() {
        return SparseRunning(output, dim);
    }

    SparseRunning sparse_running(size_t start, size_t end) {
        return SparseRunning(output + start, end - start);
    }
};

TEST(ApplyCheck, VariableDispatchRows) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_sparse(dense_row.get(), false);

    size_t N = dense_row->nrow();
    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<0>(dense_row.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 1));
    }

    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<0>(sparse_row.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 2));
    }

    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<0>(dense_column.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 3));
    }

    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<0>(sparse_column.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 4));
    }
}

TEST(ApplyCheck, VariableDispatchColumns) {
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double>(sparse_nrow, sparse_ncol, sparse_matrix));
    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_sparse(dense_row.get(), false);

    size_t N = dense_row->ncol();
    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<1>(dense_row.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 3));
    }

    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<1>(sparse_row.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 4));
    }

    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<1>(dense_column.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 1));
    }

    {
        std::vector<double> output(N);
        SillyFactory factory(output.data(), N);
        tatami::apply<1>(sparse_column.get(), factory);
        EXPECT_EQ(output, std::vector<double>(N, 2));
    }
}
