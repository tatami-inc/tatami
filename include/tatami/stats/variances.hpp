#ifndef TATAMI_VARS_HPP
#define TATAMI_VARS_HPP

#include "../base/typed_matrix.hpp"
#include "apply.hpp"

#include <vector>
#include <cmath>
#include <numeric>
#include <limits>

/**
 * @file variances.hpp
 *
 * Compute row and column variances from a `tatami::typed_matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @brief Helper to compute the variance along each dimension.
 */
struct VarianceHelper {
public:
    /**
     * This statistic can be computed from sparse inputs.
     */
    static const bool supports_sparse = true;

    /**
     * This statistic can be computed in a running manner.
     */
    static const bool supports_running = true;

private:
    static std::pair<double, double> both_NaN() {
        return std::make_pair(
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        );
    }

    static double finish_variance(double var, size_t n) {
        if (n > 1) {
            return var / (n - 1);
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

public:
    /**
     * Compute the variance along a vector.
     * We use the numerically stable calculation involvingof the sum of squared differences rather than taking the difference of squares.
     *
     * @tparam T Type of the input data.
     *
     * @param ptr Pointer to an array of values of length `n`.
     * @param n Size of the array.
     * @param buffer Unused, provided here for consistency only.
     *
     * @return The sample variance of values in `[ptr, ptr + n)`.
     * If `n <= 1`, an NaN value is returned.
     */
    template<typename T = double>
    static double compute(const T* ptr, size_t n, T* buffer = NULL) {
        return compute_with_mean(ptr, n, buffer).second;
    }

    /**
     * Compute the mean and variance along a vector.
     * This avoids a redundant calculation of the mean as it is already computed as part of the variance calculation.
     *
     * @tparam T Type of the input data.
     *
     * @param ptr Pointer to an array of values of length `n`.
     * @param n Size of the array, must be greater than 1.
     * @param buffer Unused, provided here for consistency only.
     *
     * @return The sample mean and variance of values in `[ptr, ptr + n)`.
     * If `n == 0`, the mean is set to NaN, and if `n < 2`, the variance is set to NaN.
     */
    template<typename T = double>
    static std::pair<double, double> compute_with_mean(const T* ptr, size_t n, T* buffer = NULL) {
        if (n < 1) {
            return both_NaN();
        }

        double mean = std::accumulate(ptr, ptr + n, 0.0)/n;
        double var = 0;
        for (size_t j = 0; j < n; ++j, ++ptr) {
            var += (*ptr - mean) * (*ptr - mean);
        }

        return std::make_pair(mean, finish_variance(var, n));
    }

public:
    /**
     * Compute the variance along a sparse vector.
     * This achieves faster processing by only performing summations over non-zero elements.
     *
     * @tparam T Type of the input data.
     * @tparam IDX Type of the indices.
     *
     * @param range A `sparse_range` object specifying the number and values of all non-zero indices.
     * @param n Total length of the vector, including zero values.
     * @param vbuffer,ibuffer Unused, provided here for consistency only.
     *
     * @return The sample variance of values in the vector.
     */
    template<typename T = double, typename IDX = int>
    static double compute(const sparse_range<T, IDX>& range, size_t n, T* vbuffer = NULL, IDX* ibuffer = NULL) {
        return compute_with_mean(range, n, vbuffer, ibuffer).second;
    }

    /**
     * Compute the mean and variance along a sparse vector.
     * Again, this avoids a redundant calculation of the mean.
     *
     * @tparam T Type of the input data.
     * @tparam IDX Type of the indices.
     *
     * @param range A `sparse_range` object specifying the number and values of all non-zero indices.
     * @param n Total length of the vector, including zero values.
     * @param vbuffer,ibuffer Unused, provided here for consistency only.
     *
     * @return The sample mean and variance of values in the vector.
     */
    template<typename T = double, typename IDX = int>
    static std::pair<double, double> compute_with_mean(const sparse_range<T, IDX>& range, size_t n, T* vbuffer = NULL, IDX* ibuffer = NULL) {
        if (n < 1) {
            return both_NaN();
        }

        double mean = std::accumulate(range.value, range.value + range.number, 0.0)/n;
        double var = 0;
        for (size_t j = 0; j < range.number; ++j) {
            var += (range.value[j] - mean) * (range.value[j] - mean);
        }
        var += mean * mean * (n - range.number);

        return std::make_pair(mean, finish_variance(var, n));
    }

private:
    struct Common {
    public:
        Common(size_t n) : running_var(n), running_mean(n) {}

        void finish0() {
            if (dim > 1) {
                for (auto& s : running_var) {
                    s /= dim - 1;
                }
            } else {
                if (dim == 0){
                    std::fill(running_mean.begin(), running_mean.end(), std::numeric_limits<double>::quiet_NaN());
                }
                std::fill(running_var.begin(), running_var.end(), std::numeric_limits<double>::quiet_NaN());
            }
            return;
        }
    protected:
        std::vector<double> running_var;
        std::vector<double> running_mean;
        size_t dim = 0;
    };

public:
    /**
     * @brief Helper to compute the running variance from dense inputs.
     */
    struct Dense : private Common {
        /**
         * @param n Number of parallel vectors for which to compute running statistics.
         */
        Dense(size_t n) : Common(n) {}

        /**
         * Add another vector to the running variance calculations.
         * This uses Welford's algorithm to compute the running variance of each parallel vector.
         * Each entry in the `ptr` array contains the latest values of the set of parallel vectors.
         *
         * @tparam T Type of the input data.
         *
         * @param ptr Pointer to an array of length equal to the number of parallel vectors.
         * @param buffer Ignored.
         */
        template<typename T>
        void add(const T* ptr, T* buffer = NULL) {
            auto mIt = running_mean.begin();
            ++dim;

            for (auto sIt = running_var.begin(); sIt < running_var.end(); ++sIt, ++ptr, ++mIt) {
                const double delta=*ptr - *mIt;
                *mIt += delta/dim;
                *sIt += delta*(*ptr - *mIt);
            }
            return;
        }

        /**
         * Finish the running calculations.
         */
        void finish() {
            finish0();
            return;
        }

        /**
         * Obtain the sample mean for each parallel vector. 
         * This will contain NaN values if `add()` was not called at least once.
         */
        const std::vector<double>& means() {
            return running_mean;
        }

        /**
         * Obtain the sample variance for each parallel vector. 
         * This will contain NaN values if `add()` was not called more than once.
         */
        const std::vector<double>& statistics() {
            return running_var;
        }
    };

    /**
     * @brief Helper to compute the running variance from sparse inputs.
     */
    struct Sparse : private Common {
        /**
         * @param n Number of parallel vectors for which to compute running statistics.
         */
        Sparse(size_t n) : Common(n), running_nnzero(n) {}

        /**
         * Add another sparse vector to the running variance calculations.
         * Again, we use Welford's methods, but only across the non-zero elements.
         *
         * @tparam T Type of the input data.
         * @tparam IDX Type of the indices.
         *
         * @param range A `sparse_range` object identifying the non-zero elements in the sparse vector.
         * @param vbuffer,ibuffer Ignored.
         */
        template<typename T, typename IDX>
        void add(sparse_range<T, IDX> range, T* vbuffer = NULL, IDX* ibuffer = NULL) {
            auto sIt = running_var.begin();
            auto mIt = running_mean.begin();
            auto nzIt = running_nnzero.begin();
            ++dim;

            for (size_t j = 0; j < range.number; ++j, ++range.index, ++range.value) {
                auto ri = *range.index;
                auto& curM = *(mIt + ri);
                auto& curS = *(sIt + ri);
                auto& curNZ = *(nzIt + ri);
                ++curNZ;

                const auto& curval = *range.value;
                const double delta = curval - curM;
                curM += delta / curNZ;
                curS += delta * (curval - curM);
            }

            return;
        }

        /**
         * Finish the running calculations.
         */
        void finish() {
            if (dim) {
                auto mIt = running_mean.begin();
                auto nzIt = running_nnzero.begin();
                for (auto sIt = running_var.begin(); sIt != running_var.end(); ++mIt, ++sIt, ++nzIt) {
                    const double curNZ = *nzIt;
                    const double ratio = curNZ / dim;
                    auto& curM = *mIt;
                    *sIt += curM * curM * ratio * (dim - curNZ);
                    curM *= ratio;
                }
            }

            finish0();
            return;
        }

        /**
         * Obtain the sample mean for each parallel vector. 
         * This will contain NaN values if `add()` was not called at least once.
         */
        const std::vector<double>& means() {
            return running_mean;
        }

        /**
         * Obtain the sample variance for each parallel vector. 
         * This will contain NaN values if `add()` was not called more than once.
         */
        const std::vector<double>& statistics() {
            return running_var;
        }
    private:
        std::vector<int> running_nnzero;
    };
};

}

/**
 * This uses the usual algorithm for matrices where `tatami::matrix::prefer_rows()` is false, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam T Type of the matrix value, should be numeric.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column variances.
 */
template<typename T, typename IDX>
inline std::vector<T> column_variances(const typed_matrix<T, IDX>* p) {
    return apply<1, T, IDX, stats::VarianceHelper>(p);
}

/**
 * This uses the usual algorithm for matrices where `tatami::matrix::prefer_rows()` is true, otherwise it uses Welford's algorithm.
 * As a result, the computed variances will be slightly different (within numerical precision) for row- and column-major matrices.
 *
 * @tparam T Type of the matrix value, should be numeric.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row variances.
 */
template<typename T, typename IDX>
inline std::vector<T> row_variances(const typed_matrix<T, IDX>* p) {
    return apply<0, T, IDX, stats::VarianceHelper>(p);
}

}

#endif
