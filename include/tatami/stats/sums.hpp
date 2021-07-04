#ifndef TATAMI_STATS_SUMS_HPP
#define TATAMI_STATS_SUMS_HPP

#include "../base/typed_matrix.hpp"
#include "apply.hpp"
#include <vector>
#include <numeric>

/**
 * @file sums.hpp
 *
 * Compute row and column sums from a `tatami::typed_matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @brief Helper to compute the sum along each dimension.
 */
struct SumHelper {
public:
    /**
     * Type of the computed statistic.
     */
    typedef double value;

    /**
     * This statistic can be computed from sparse inputs.
     */
    static const bool supports_sparse = true;

    /**
     * This statistic can be computed in a running manner.
     */
    static const bool supports_running = true;

public:
    /**
     * Compute the mean along a vector.
     *
     * @tparam T Type of the input data.
     *
     * @param ptr Pointer to an array of values of length `n`.
     * @param n Size of the array.
     * @param buffer Unused, provided here for consistency only.
     *
     * @return The sample mean of values in `[ptr, ptr + n)`.
     */
    template<typename T = double>
    static double compute(const T* ptr, size_t n, T* buffer = NULL) {
        return std::accumulate(ptr, ptr + n, 0.0);
    }

    /**
     * Compute the mean along a sparse vector.
     * This achieves faster processing by only performing summations over non-zero elements.
     *
     * @tparam T Type of the input data.
     * @tparam IDX Type of the indices.
     *
     * @param range A `sparse_range` object specifying the number and values of all non-zero indices.
     * @param n Total length of the vector, including zero values.
     * @param vbuffer,ibuffer Unused, provided here for consistency only.
     *
     * @return The sample mean of values in the vector.
     */
    template<typename T = double, typename IDX = int>
    static double compute(const sparse_range<T, IDX>& range, size_t n, T* vbuffer = NULL, IDX* ibuffer = NULL) {
        return std::accumulate(range.value, range.value + range.number, 0.0);
    }

public:
    /**
     * @brief Helper to compute the running sum from dense inputs.
     */
    struct Dense {
        /**
         * @param n Number of parallel vectors for which to compute running statistics.
         */
        Dense(size_t n) : store(n) {}

        /**
         * Add another vector to the running sum calculations.
         * Each entry in the `ptr` array contains the latest values of the set of parallel vectors.
         *
         * @tparam T Type of the input data.
         *
         * @param ptr Pointer to an array of length equal to the number of parallel vectors.
         * @param buffer Ignored.
         */
        template<typename T = double>
        void add(const T* ptr, T* buffer = NULL) {
            for (auto sIt = store.begin(); sIt != store.end(); ++sIt, ++ptr) {
                *sIt += *ptr;
            }
            return;
        }

        /**
         * Finish the running calculation.
         */
        void finish() {}

        /**
         * Obtain the sum of values for each parallel vector. 
         */
        const std::vector<double>& statistics() {
            return store;
        }

    private:
        std::vector<double> store;
    };

    /**
     * @brief Helper to compute the running sum from sparse inputs.
     */
    struct Sparse {
        /**
         * @param n Number of parallel vectors for which to compute running statistics.
         */
        Sparse(size_t n) : store(n) {}

        /**
         * Add another sparse vector to the running sum calculations.
         *
         * @tparam T Type of the input data.
         * @tparam IDX Type of the indices.
         *
         * @param range A `sparse_range` object identifying the non-zero elements in the sparse vector.
         * @param vbuffer,ibuffer Ignored.
         */
        template<typename T = double, typename IDX = int>
        void add(sparse_range<T, IDX> range, T* vbuffer = NULL, IDX* ibuffer = NULL) {
            for (size_t j = 0; j < range.number; ++j, ++range.index, ++range.value) {
                store[*range.index] += *range.value;
            }
            return;
        }

        /**
         * Finish the running calculation.
         */
        void finish() {}

        /**
         * Obtain the sum of values for each parallel vector. 
         */
        const std::vector<double>& statistics() {
            return store;
        }
    private:
        std::vector<double> store;
    };
};

}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column sums.
 */
template<typename T, typename IDX>
inline std::vector<T> column_sums(const typed_matrix<T, IDX>* p) {
    return apply<1, T, IDX, stats::SumHelper>(p);
}

/**
 * @tparam T Type of the matrix value, should be summable.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::typed_matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row sums.
 */
template<typename T, typename IDX>
inline std::vector<T> row_sums(const typed_matrix<T, IDX>* p) {
    return apply<0, T, IDX, stats::SumHelper>(p);
}

}

#endif
