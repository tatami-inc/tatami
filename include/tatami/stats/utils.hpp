#ifndef TATAMI_STATS_UTILS_HPP
#define TATAMI_STATS_UTILS_HPP

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "../utils/Oracles.hpp"

#include <cmath>

/**
 * @file utils.hpp
 *
 * @brief Utilities for computing matrix statistics.
 */

namespace tatami {

/**
 * Apply a function to a set of tasks, distributing them to threads via OpenMP if enabled.
 * Callers can specify a custom parallelization scheme by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro, which should accept four arguments:
 * `fun`, `tasks` and `threads` as below, as well as a `worker_size` argument that specifies the size of each task range
 * (the last is provided for convenience only, to avoid the need to re-compute it inside the macro).
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current thread.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `thread`, the thread number executing this task range.
 * - `start`, the start index of the task range.
 * - `end`, the first index past the end of the task range.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * @param tasks Number of tasks.
 * @param threads Number of threads.
 */
template<bool parallel_ = true, class Function_>
void parallelize(Function_ fun, size_t tasks, size_t threads) {
#if defined(_OPENMP) || defined(TATAMI_CUSTOM_PARALLEL)
    if constexpr(parallel_) {
        size_t worker_size = std::ceil(static_cast<double>(tasks) / static_cast<double>(threads));

        if (threads > 1) {
#ifndef TATAMI_CUSTOM_PARALLEL
            #pragma omp parallel for num_threads(threads)
            for (size_t t = 0; t < threads; ++t) {
                size_t start = worker_size * t, end = std::min(dim, start + worker_size);
                if (start < end) {
                    fun(t, start, end);
                }
            }
#else
            TATAMI_CUSTOM_PARALLEL(std::move(fun), tasks, threads, worker_size);
#endif
            return;
        }
    }
#endif

    fun(0, 0, tasks);
    return;
}

/**
 * @cond
 */
template<bool row_, bool sparse_, typename Value_, typename Index_>
auto direct_extractor(const Matrix<Value_, Index_>* mat, bool uses_oracle, Index_ from, Index_ to, const Options& opt) {
    if constexpr(sparse_) {
        auto ext = (row_ ? mat->sparse_row(opt) : mat->sparse_column(opt));
        if (uses_oracle) {
            ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(from, to - from));
        }
        return ext;
    } else {
        auto ext = (row_ ? mat->dense_row(opt) : mat->dense_column(opt));
        if (uses_oracle) {
            ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(from, to - from));
        }
        return ext;
    }
}
/**
 * @endcond
 */

/**
 * @brief Configuration for bidimensional apply.
 *
 * Set up the iteration across a `Matrix` in its preferred iteration dimension.
 * More specifically, we consider each matrix to be a collection of equi-length target vectors, where we wish to compute a statistic for each target vector.
 * We can do so directly by iterating over the target vectors; or we can do so in a "running" manner, where statistics for target vectors can be computed by iterating over the non-target vectors.
 * Use of this class assumes that the desired statistic can be computed in a running manner, which is helpful if the preferred iteration dimension is not the same as that required to iterate directly over the target vectors.
 *
 * @tparam row_ Whether each row is a target vector, e.g., to compute row-wise statistics.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column index.
 */
template<bool row_, typename Value_, typename Index_>
struct BidimensionalApplyConfiguration {
    /**
     * @param mat Pointer to the `Matrix` instance.
     */
    BidimensionalApplyConfiguration(const Matrix<Value_, Index_>* mat) :
        matrix(mat), 
        prefer_rows(matrix->prefer_rows()),
        uses_oracle(matrix->uses_oracle(prefer_rows)),
        target_dim(row_ ? matrix->nrow() : matrix->ncol()),
        other_dim(row_ ? matrix->ncol() : matrix->nrow())
    {}

private:
    const Matrix<Value_, Index_>* matrix;

public:
    /**
     * Whether the matrix prefers iteration over the rows.
     * Same as `Matrix::prefer_rows()` but is just reported here for convenience.
     */
    bool prefer_rows;

private:
    bool uses_oracle; // need to do this little shuffle to avoid problems with initialization order.

public:
    /**
     * Extent of the dimension containing the target vectors, i.e., the number of target vectors.
     */
    Index_ target_dim;

    /**
     * Extent of the dimension containing the non-target vectors, i.e., the "other" dimension, for which we don't want to compute statistics.
     */
    Index_ other_dim;

    /**
     * @param s Index of the first target vector.
     * @param e First index past the last target vector.
     * @param opt Options to pass to the extractor construction methods in `Matrix`.
     *
     * @return An `Extractor` object for iterating over the target dimension in the range `[s, e)`.
     */
    template<bool sparse_>
    auto direct(Index_ s, Index_ e, const Options& opt) const {
        return direct_extractor<row_, sparse_>(matrix, uses_oracle, s, e, opt);
    }

    /**
     * @param s Index of the first target vector.
     * @param e First index past the last target vector.
     * @param opt Options to pass to the extractor construction methods in `Matrix`.
     *
     * @return An `Extractor` object for iterating over the non-target dimensions,
     * extracting the range `[s, e)` from each target vector.
     */
    template<bool sparse_>
    auto running(Index_ s, Index_ e, const Options& opt) const {
        Index_ len = e - s;
        if constexpr(sparse_) {
            auto ext = (row_ ? matrix->sparse_column(s, len, opt) : matrix->sparse_row(s, len, opt));
            if (uses_oracle) {
                ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(0, other_dim));
            }
            return ext;
        } else {
            auto ext = (row_ ? matrix->dense_column(s, len, opt) : matrix->dense_row(s, len, opt));
            if (uses_oracle) {
                ext->set_oracle(std::make_unique<ConsecutiveOracle<Index_> >(0, other_dim));
            }
            return ext;
        }
    }
};

/**
 * @brief Configuration for direct apply.
 *
 * Set up the iteration across a `Matrix` in the specified iteration dimension.
 * More specifically, we consider each matrix to be a collection of equi-length target vectors, where we wish to compute a statistic for each target vector.
 * We do so directly by iterating over the target vectors, ignoring any preference in `Matrix::prefer_rows()`.
 * This is useful when a statistic cannot be computed in a running manner.
 *
 * @tparam row_ Whether each row is a target vector, e.g., to compute row-wise statistics.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column index.
 */
template<bool row_, typename Value_, typename Index_>
struct DirectApplyConfiguration {
    DirectApplyConfiguration(const Matrix<Value_, Index_>* mat) :
        matrix(mat), 
        uses_oracle(matrix->uses_oracle(row_)),
        target_dim(row_ ? matrix->nrow() : matrix->ncol()),
        other_dim(row_ ? matrix->ncol() : matrix->nrow())
    {}

private:
    const Matrix<Value_, Index_>* matrix;
    bool uses_oracle;

public:
    /**
     * Extent of the dimension containing the target vectors, i.e., the number of target vectors.
     */
    Index_ target_dim;

    /**
     * Extent of the dimension containing the non-target vectors, i.e., the "other" dimension, for which we don't want to compute statistics.
     */
    Index_ other_dim;

    /**
     * @param s Index of the first target vector.
     * @param e First index past the last target vector.
     * @param opt Options to pass to the extractor construction methods in `Matrix`.
     *
     * @return An `Extractor` object for iterating over the target dimension in the range `[s, e)`.
     */
    template<bool sparse_>
    auto extractor(Index_ s, Index_ e, const Options& opt) const {
        return direct_extractor<row_, sparse_>(matrix, uses_oracle, s, e, opt);
    }
};

}

#endif
