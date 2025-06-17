#ifndef TATAMI_ORACLE_HPP
#define TATAMI_ORACLE_HPP

#include <cstddef>

/**
 * @file Oracle.hpp
 *
 * @brief Oracle for data access.
 */

namespace tatami {

/**
 * Unsigned integer type for the prediction index.
 * This is used to obtain the `i`-th prediction from the `Oracle`.
 */
typedef std::size_t PredictionIndex;

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future access requests on the target dimension.
 *
 * This allows `Matrix` implementations to pre-fetch data for future requests to `OracularDenseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
 * Check out `ConsecutiveOracle` and `FixedVectorOracle` for some examples of concrete subclasses.
 */
template<typename Index_>
class Oracle {
public:
    /**
     * @cond
     */
    Oracle() = default;
    Oracle(const Oracle&) = default;
    Oracle& operator=(const Oracle&) = default;
    Oracle(Oracle&&) = default;
    Oracle& operator=(Oracle&&) = default;
    virtual ~Oracle() = default;
    /**
     * @endcond
     */

    /**
     * @return Total number of predictions.
     */
    virtual PredictionIndex total() const = 0;

    /**
     * @param i Which prediction to return.
     * @return The `i`-th prediction, to be interpreted as an index on the target dimension.
     */
    virtual Index_ get(PredictionIndex i) const = 0;
};

}

#endif
