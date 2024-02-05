#ifndef TATAMI_ORACLE_HPP
#define TATAMI_ORACLE_HPP

/**
 * @file Oracle.hpp
 *
 * @brief Oracle for data access.
 */

namespace tatami {

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future access requests.
 *
 * This allows `Matrix` implementations to pre-fetch data for future requests to `DenseOracleAwareExtractor::fetch()` or `SparseOracleAwareExtractor::fetch()`.
 */
template<typename Index_>
struct Oracle {
    /**
     * @cond
     */
    virtual ~Oracle() = default;
    /**
     * @endcond
     */

    /**
     * @return Total number of predictions.
     */
    virtual size_t total() const = 0;

    /**
     * @param i Which prediction to return.
     * @return The index of the `i`-th prediction.
     */
    virtual Index_ get(size_t i) const = 0; 
};

}

#endif
