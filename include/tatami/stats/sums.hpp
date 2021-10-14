#ifndef TATAMI_STATS_SUMS_HPP
#define TATAMI_STATS_SUMS_HPP

#include "../base/Matrix.hpp"
#include "apply.hpp"
#include <vector>
#include <numeric>

/**
 * @file sums.hpp
 *
 * Compute row and column sums from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

template<typename O>
struct SumFactory {
public:
    SumFactory(O* o, size_t d1, size_t d2) : output(o), dim(d1), otherdim(d2) {}

private:
    O* output;
    size_t dim, otherdim;

public:
    struct DenseDirect {
        DenseDirect(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename V>
        void compute(size_t i, const V* ptr) {
            output[i] = std::accumulate(ptr, ptr + otherdim, static_cast<O>(0));
        }
    private:
        O* output;
        size_t otherdim;
    };

    DenseDirect dense_direct() {
        return DenseDirect(output, otherdim);
    }

public:
    struct SparseDirect {
        SparseDirect(O* o) : output(o) {}

        template<typename T, typename IDX>
        void compute(size_t i, const SparseRange<T, IDX>& range) {
            output[i] = std::accumulate(range.value, range.value + range.number, static_cast<O>(0));
        }
    private:
        O* output;
    };

    SparseDirect sparse_direct() {
        return SparseDirect(output);
    }

public:
    struct DenseRunning {
        DenseRunning(O* o, size_t d1) : output(o), dim(d1) {}

        template<typename V>
        void add(const V* ptr) {
            for (size_t d = 0; d < dim; ++d) {
                output[d] += ptr[d];
            }
        }

        void finish() {}
    private:
        O* output;
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
        SparseRunning(O* o) : output(o) {}

        template<typename T = double, typename IDX = int>
        void add(const SparseRange<T, IDX>& range) {
            for (size_t j = 0; j < range.number; ++j) {
                output[range.index[j]] += range.value[j];
            }
        }

        void finish() {}
    private:
        O* output;
    };

    SparseRunning sparse_running() {
        return SparseRunning(output);
    }

    SparseRunning sparse_running(size_t start, size_t end) {
        return SparseRunning(output);
    }
};

}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column sums.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> column_sums(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->ncol());
    stats::SumFactory factory(output.data(), p->ncol(), p->nrow());
    apply<1>(p, factory);
    return output;
}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row sums.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> row_sums(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->nrow());
    stats::SumFactory factory(output.data(), p->nrow(), p->ncol());
    apply<0>(p, factory);
    return output;
}

}

#endif
