#ifndef TATAMI_STATS_VARS_HPP
#define TATAMI_STATS_VARS_HPP

#include "../base/Matrix.hpp"
#include "apply.hpp"

#include <vector>
#include <cmath>
#include <numeric>
#include <limits>

/**
 * @file variances.hpp
 *
 * Compute row and column variances from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

namespace variances {

/**
 * @cond
 */
template<typename O>
std::pair<O, O> both_NaN() {
    return std::make_pair(
        std::numeric_limits<O>::quiet_NaN(),
        std::numeric_limits<O>::quiet_NaN()
    );
}

template<typename O> 
O finish_variance_direct(O var, size_t n) {
    if (n > 1) {
        return var / (n - 1);
    } else {
        return std::numeric_limits<O>::quiet_NaN();
    }
}
/**
 * @endcond
 */

/**
 * Compute the mean and variance along a vector.
 * We use the numerically stable calculation involving the sum of squared differences rather than taking the difference of squares.
 *
 * @tparam O Type of the output data.
 * @tparam T Type of the input data.
 *
 * @param[in] ptr Pointer to an array of values of length `n`.
 * @param n Size of the array.
 *
 * @return The sample mean and variance of values in `[ptr, ptr + n)`.
 * If `n == 0`, the mean is set to NaN, and if `n < 2`, the variance is set to NaN.
 */
template<typename O = double, typename T>
std::pair<O, O> compute_direct(const T* ptr, size_t n) {
    if (n < 1) {
        return both_NaN<O>();
    }

    O mean = std::accumulate(ptr, ptr + n, static_cast<O>(0))/n;
    O var = 0;
    for (size_t j = 0; j < n; ++j, ++ptr) {
        var += (*ptr - mean) * (*ptr - mean);
    }

    return std::make_pair(mean, finish_variance_direct(var, n));
}

/**
 * Compute the mean and variance along a sparse vector.
 * This achieves faster processing by only performing summations over non-zero elements.
 *
 * @tparam O Type of the output data.
 * @tparam T Type of the input data.
 * @tparam IDX Type of the indices.
 *
 * @param range A `SparseRange` object specifying the number and values of all non-zero indices.
 * @param n Total length of the vector, including zero values.
 *
 * @return The sample mean and variance of values in the vector.
 * If `n == 0`, the mean is set to NaN, and if `n < 2`, the variance is set to NaN.
 */
template<typename O = double, typename T, typename IDX>
std::pair<O, O> compute_direct(const SparseRange<T, IDX>& range, size_t n) {
    if (n < 1) {
        return both_NaN<O>();
    }

    O mean = std::accumulate(range.value, range.value + range.number, static_cast<O>(0))/n;
    O var = 0;
    for (size_t j = 0; j < range.number; ++j) {
        var += (range.value[j] - mean) * (range.value[j] - mean);
    }
    var += mean * mean * (n - range.number);

    return std::make_pair(mean, finish_variance_direct(var, n));
}

/**
 * Compute a running mean and variance from dense data using Welford's method.
 * This considers a set of "target vectors" that grow by one with each call to this function;
 * the mean and variance of each target vector is computed. 
 *
 * @tparam T Type of the input data.
 * @tparam O Type of the output data.
 *
 * @param[in] ptr Pointer to an array of values of length `n`, with one value per target vector.
 * @param n Size of the array referenced by `ptr`.
 * @param[out] means Pointer to an array containing the running means for each target vector.
 * @param[out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param count Number of times this function has already been called.
 *
 * @return `means` and `vars` are updated with the corresponding elements from `ptr`.
 * `count` is incremented by 1.
 */
template<typename T, typename O>
void compute_running(const T* ptr, size_t n, O* means, O* vars, int& count) { 
    ++count;
    for (size_t i = 0; i < n; ++i, ++means, ++vars, ++ptr) {
        const double delta=*ptr - *means;
        *means += delta/count;
        *vars += delta*(*ptr - *means);
    }
}

/**
 * Compute a running mean and variance from sparse data using Welford's method.
 * This considers a set of "target vectors" that grow by one with each call to this function;
 * the mean and variance of each target vector is computed. 
 *
 * @tparam T Type of the input data.
 * @tparam IDX Type of the indices.
 * @tparam O Type of the output data.
 * @tparam Nz Type fo the non-zero counts.
 *
 * @param range A `SparseRange` object specifying the number, indices and values of all non-zero elements.
 * Each element is assigned to a target vector based on its index.
 * @param[out] means Pointer to an array containing the running means for each target vector.
 * @param[out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param[out] nonzeros Pointer to an array containing the running number of non-zero values for each target vector.
 * @param count Number of times this function has already been called.
 *
 * @return `means` and `vars` are updated with the corresponding elements from `range`.
 * `count` is incremented by 1.
 */
template<typename T, typename IDX, typename O, typename Nz>
void compute_running(const SparseRange<T, IDX>& range, O* means, O* vars, Nz* nonzeros, int& count) {
    ++count;
    for (size_t j = 0; j < range.number; ++j) {
        if (range.value[j]) { // skip observed zeros so 'nonzeros' is what it says.
            auto ri = range.index[j];
            auto& curM = means[ri];
            auto& curS = vars[ri];
            auto& curNZ = nonzeros[ri];
            ++curNZ;

            const auto& curval = range.value[j];
            const double delta = curval - curM;
            curM += delta / curNZ;
            curS += delta * (curval - curM);
        }
    }
}

/**
 * Finish the running variance calculations from `compute_running()`.
 *
 * @tparam O Type of the output data.
 *
 * @param n Number of target vectors.
 * @param[out] means Pointer to an array containing the running means for each target vector.
 * @param[out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param count Number of times the `compute_running()` function was called.
 *
 * @return Elements in `vars` are divided by `count - 1`.
 * For low values of `count`, `means` and/or `vars` may be filled with NaNs.
 */
template<typename O>
void finish_running(size_t n, O* means, O* vars, int count) {
    if (count > 1) {
        for (size_t i = 0; i < n; ++i) {
            vars[i] /= count - 1;
        }
    } else {
        if (count == 0){
            std::fill(means, means + n, std::numeric_limits<double>::quiet_NaN());
        }
        std::fill(vars, vars + n, std::numeric_limits<double>::quiet_NaN());
    }
}

/**
 * Finish the running variance calculations from `compute_running()` with sparse inputs.
 *
 * @tparam O Type of the output data.
 * @tparam Nz Type fo the non-zero counts.
 *
 * @param n Number of target vectors.
 * @param[out] means Pointer to an array containing the running means for each target vector.
 * @param[out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param[out] nonzeros Pointer to an array containing the running number of non-zero values for each target vector.
 * @param count Number of times the `compute_running()` function was called.
 *
 * @return Elements in `vars` are divided by `count - 1`.
 * For low values of `count`, `means` and/or `vars` may be filled with NaNs.
 */
template<typename O, typename Nz>
void finish_running(size_t n, O* means, O* vars, Nz* nonzeros, int count) {
    if (count) {
        for (size_t i = 0; i < n; ++i) {
            const double curNZ = nonzeros[i];
            const double ratio = curNZ / count;
            auto& curM = means[i];
            vars[i] += curM * curM * ratio * (count - curNZ);
            curM *= ratio;
        }
    }
    finish_running(n, means, vars, count);
}

}

template<typename O = double>
struct VarianceFactory {
public:
    VarianceFactory(O* o, size_t d1, size_t d2) : output(o), dim(d1), otherdim(d2) {}

private:
    O* output;
    size_t dim, otherdim;

public:
    struct DenseDirect {
        DenseDirect(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename V>
        void compute(size_t i, const V* ptr) {
            output[i] = variances::compute_direct<O>(ptr, otherdim).second;
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
        SparseDirect(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename T, typename IDX>
        void compute(size_t i, const SparseRange<T, IDX>& range) {
            output[i] = variances::compute_direct<O>(range, otherdim).second;
        }
    private:
        O* output;
        size_t otherdim;
    };

    SparseDirect sparse_direct() {
        return SparseDirect(output, otherdim);
    }

public:
    struct DenseRunning {
        DenseRunning(O* o, size_t d1) : output(o), dim(d1), running_means(dim) {}

        template<typename V>
        void add(const V* ptr) {
            variances::compute_running(ptr, dim, running_means.data(), output, counter);
        }

        void finish() {
            variances::finish_running(dim, running_means.data(), output, counter);
        }
    private:
        O* output;
        size_t dim;
        std::vector<O> running_means;
        int counter = 0;
    };

    DenseRunning dense_running() {
        return DenseRunning(output, dim);
    }

    DenseRunning dense_running(size_t start, size_t end) {
        return DenseRunning(output + start, end - start);
    }

public:
    struct SparseRunning {
        SparseRunning(O* o, size_t dim, size_t s, size_t e) : output(o), start(s), end(e), running_means(dim), running_nzeros(dim) {}

        template<typename T, typename IDX>
        void add(const SparseRange<T, IDX>& range) {
            variances::compute_running(range, running_means.data(), output, running_nzeros.data(), counter);
        }

        void finish() {
            variances::finish_running(end - start, running_means.data() + start, output + start, running_nzeros.data() + start, counter);
        }
    private:
        O* output;
        size_t start, end;
        std::vector<O> running_means;
        std::vector<int> running_nzeros;
        int counter = 0;
    };

    SparseRunning sparse_running() {
        return SparseRunning(output, dim, 0, dim);
    }

    SparseRunning sparse_running(size_t start, size_t end) {
        return SparseRunning(output, dim, start, end);
    }
};

}

/**
 * This uses the usual algorithm for matrices where `tatami::Matrix::prefer_rows()` is false, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value, should be numeric.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column variances.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> column_variances(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->ncol());
    stats::VarianceFactory factory(output.data(), p->ncol(), p->nrow());
    apply<1>(p, factory);
    return output;
}

/**
 * This uses the usual algorithm for matrices where `tatami::Matrix::prefer_rows()` is true, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value, should be numeric.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row variances.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> row_variances(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->nrow());
    stats::VarianceFactory factory(output.data(), p->nrow(), p->ncol());
    apply<0>(p, factory);
    return output;
}

}

#endif
