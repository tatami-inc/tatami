#ifndef TATAMI_STATS_VARS_HPP
#define TATAMI_STATS_VARS_HPP

#include "../base/Matrix.hpp"
#include "utils.hpp"

#include <vector>
#include <cmath>
#include <numeric>
#include <limits>

/**
 * @file variances.hpp
 *
 * @brief Compute row and column variances from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

namespace variances {

/**
 * @cond
 */
template<typename Output_>
std::pair<Output_, Output_> both_NaN() {
    return std::make_pair(
        std::numeric_limits<Output_>::quiet_NaN(),
        std::numeric_limits<Output_>::quiet_NaN()
    );
}

template<typename Output_> 
Output_ finish_variance_direct(Output_ var, size_t n) {
    if (n > 1) {
        return var / (n - 1);
    } else {
        return std::numeric_limits<Output_>::quiet_NaN();
    }
}
/**
 * @endcond
 */

/**
 * Compute the mean and variance along a vector.
 * We use the numerically stable calculation involving the sum of squared differences rather than taking the difference of squares.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 *
 * @param[in] ptr Pointer to an array of values of length `n`.
 * @param n Size of the array.
 *
 * @return The sample mean and variance of values in `[ptr, ptr + n)`.
 * If `n == 0`, the mean is set to NaN, and if `n < 2`, the variance is set to NaN.
 */
template<typename Output_ = double, typename Value_>
std::pair<Output_, Output_> compute_direct(const Value_* ptr, size_t n) {
    if (n < 1) {
        return both_NaN<Output_>();
    }

    Output_ mean = std::accumulate(ptr, ptr + n, static_cast<Output_>(0))/n;
    Output_ var = 0;
    for (size_t j = 0; j < n; ++j, ++ptr) {
        var += (*ptr - mean) * (*ptr - mean);
    }

    return std::make_pair(mean, finish_variance_direct(var, n));
}

/**
 * Compute the mean and variance along a sparse vector.
 * This achieves faster processing by only performing summations over non-zero elements.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the indices.
 *
 * @param range A `SparseRange` object specifying the number and values of all non-zero indices.
 * @param n Total length of the vector, including zero values.
 *
 * @return The sample mean and variance of values in the vector.
 * If `n == 0`, the mean is set to NaN, and if `n < 2`, the variance is set to NaN.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<Output_, Output_> compute_direct(const SparseRange<Value_, Index_>& range, size_t n) {
    if (n < 1) {
        return both_NaN<Output_>();
    }

    Output_ mean = std::accumulate(range.value, range.value + range.number, static_cast<Output_>(0))/n;
    Output_ var = 0;
    for (size_t j = 0; j < range.number; ++j) {
        var += (range.value[j] - mean) * (range.value[j] - mean);
    }
    var += mean * mean * (n - range.number);

    return std::make_pair(mean, finish_variance_direct(var, n));
}

/**
 * Compute running means and variances from dense data using Welford's method.
 * This considers a scenario involving a set of equilength "target" vectors [V1, V2, V3, ..., Vn],
 * where `n` is defined as below and `ptr[i]` contains the `count`-th element for target vector Vi.
 * The aim is to compute the mean and variance for each target vector,
 * by repeatedly calling this function with different `ptr` containing successive elements of all target vectors.
 *
 * @tparam Value_ Type of the input data.
 * @tparam Output_ Type of the output data.
 *
 * @param[in] ptr Pointer to an array of values of length `n`, with one value per target vector.
 * @param n Size of the array referenced by `ptr`.
 * @param[in,out] means Pointer to an array containing the running means for each target vector.
 * @param[in,out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param count Number of times this function has already been called.
 * This is incremented by 1 upon return.
 */
template<typename Value_, typename Output_>
void compute_running(const Value_* ptr, size_t n, Output_* means, Output_* vars, int& count) { 
    ++count;
    for (size_t i = 0; i < n; ++i, ++means, ++vars, ++ptr) {
        const double delta=*ptr - *means;
        *means += delta/count;
        *vars += delta*(*ptr - *means);
    }
}

/**
 * Compute a running mean and variance from sparse data using Welford's method.
 * This does the same as its dense overload.
 *
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the indices.
 * @tparam Output_ Type of the output data.
 * @tparam Nonzero_ Type fo the non-zero counts.
 *
 * @param range A `SparseRange` object specifying the number, indices and values of all non-zero elements.
 * Each element is assigned to a target vector based on its index.
 * @param[in,out] means Pointer to an array containing the running means for each target vector.
 * @param[in,out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param[in,out] nonzeros Pointer to an array containing the running number of non-zero values for each target vector.
 * @param count Number of times this function has already been called.
 * @param skip_zeros Whether non-structural zeros in `range.value` should be skipped.
 * If `false`, the output `nonzeros` instead contains the number of _structural_ non-zero values in each target vector,
 * which may be useful for informing further operations on the compressed sparse matrix structure.
 * Note that this choice has no effect on the computed means or variances, besides some differences due to numeric imprecision.
 * @param subtract Value to subtract from each sparse index before using it to index into `means`, `vars` and `nonzeros`.
 * This is only relevant if `means` and friends do not hold statistics for all target vectors, but just a contiguous block, e.g., during parallelization.
 */
template<typename Value_, typename Index_, typename Output_, typename Nonzero_>
void compute_running(const SparseRange<Value_, Index_>& range, Output_* means, Output_* vars, Nonzero_* nonzeros, int& count, bool skip_zeros = true, Index_ subtract = 0) {
    ++count;
    for (size_t j = 0; j < range.number; ++j) {
        if (!skip_zeros || range.value[j]) { 
            auto ri = range.index[j] - subtract;
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
 * Elements in `vars` are divided by `count - 1`, unless if `count < 2`, in which case it is filled with NaNs.
 * If `count` is zero, `means` will also be filled with NaNs.
 *
 * @tparam Output_ Type of the output data.
 *
 * @param n Number of target vectors.
 * @param[in,out] means Pointer to an array containing the running means for each target vector.
 * @param[in,out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param count Number of times the `compute_running()` function was called.
 *
 */
template<typename Output_>
void finish_running(size_t n, Output_* means, Output_* vars, int count) {
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
 * Elements in `vars` are divided by `count - 1`, unless if `count < 2`, in which case it is filled with NaNs.
 * If `count` is zero, `means` will also be filled with NaNs.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Nonzero_ Type of the non-zero counts.
 *
 * @param n Number of target vectors.
 * @param[in,out] means Pointer to an array containing the running means for each target vector.
 * @param[in,out] vars Pointer to an array containing the running sum of squared differences from the mean for each target vector.
 * @param[in] nonzeros Pointer to an array containing the running number of non-zero values for each target vector.
 * @param count Number of times the `compute_running()` function was called.
 */
template<typename Output_, typename Nonzero_>
void finish_running(size_t n, Output_* means, Output_* vars, const Nonzero_* nonzeros, int count) {
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

/**
 * @cond
 */
template<typename Output_, bool row_, typename Value_, typename Index_>
std::vector<Output_> dimension_variances(const Matrix<Value_, Index_>* p, int threads) {
    auto dim = (row_ ? p->nrow() : p->ncol());
    std::vector<Output_> output(dim);
    auto otherdim = (row_ ? p->ncol() : p->nrow());
    const bool direct = p->prefer_rows() == row_;

    if (p->sparse()) {
        if (direct) {
            Options opt;
            opt.sparse_extract_index = false;
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<row_, true>(p, s, l);
                std::vector<Value_> vbuffer(otherdim);
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), NULL);
                    output[i] = variances::compute_direct<Output_>(out, otherdim).second;
                }
            }, dim, threads);

        } else {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<!row_, true>(p, 0, otherdim, s, l);
                auto len = ext->block_length;
                std::vector<Value_> vbuffer(len);
                std::vector<Index_> ibuffer(len);

                std::vector<Output_> running_means(len);
                std::vector<Index_> running_nzeros(len);
                auto running_vars = output.data() + s;
                int counter = 0;

                for (Index_ i = 0; i < otherdim; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), ibuffer.data());
                    variances::compute_running(out, running_means.data(), running_vars, running_nzeros.data(), counter, /* skip_zeros = */ true, /* subtract = */ s);
                }
                variances::finish_running(len, running_means.data(), running_vars, running_nzeros.data(), counter);
            }, dim, threads);
        }

    } else {
        if (direct) {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<row_, false>(p, s, l);
                std::vector<Value_> buffer(otherdim);
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto out = ext->fetch(i, buffer.data());
                    output[i] = variances::compute_direct<Output_>(out, otherdim).second;
                }
            }, dim, threads);

        } else {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<!row_, false>(p, 0, otherdim, s, l);
                auto len = ext->block_length;
                std::vector<Value_> buffer(len);

                std::vector<Output_> running_means(len);
                auto running_vars = output.data() + s;
                int counter = 0;

                for (Index_ i = 0; i < otherdim; ++i) {
                    auto out = ext->fetch(i, buffer.data());
                    variances::compute_running(out, len, running_means.data(), running_vars, counter);
                }
                variances::finish_running(len, running_means.data(), running_vars, counter);
            }, dim, threads);
        }
    }

    return output;
}
/**
 * @endcond
 */

}

/**
 * This uses the usual algorithm for matrices where `tatami::Matrix::prefer_rows()` is false, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the column variances.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> column_variances(const Matrix<Value_, Index_>* p, int threads = 1) {
    return stats::dimension_variances<Output_, false>(p, threads);
}

/**
 * This uses the usual algorithm for matrices where `tatami::Matrix::prefer_rows()` is true, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the row variances.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> row_variances(const Matrix<Value_, Index_>* p, int threads = 1) {
    return stats::dimension_variances<Output_, true>(p, threads);
}

}

#endif
