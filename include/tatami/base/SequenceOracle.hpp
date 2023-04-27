#ifndef TATAMI_SEQUENCE_ORACLE_HPP
#define TATAMI_SEQUENCE_ORACLE_HPP

#include <vector>
#include <memory>
#include <numeric>

/**
 * @file SequenceOracle.hpp
 *
 * @brief Predict future accesses during iteration.
 */

namespace tatami {

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future access requests.
 *
 * This allows `Matrix` implementations to pre-fetch data for future requests to `DenseExtractor::fetch()` or `SparseExtractor::fetch()`.
 */
template<typename Index_>
struct SequenceOracle {
    /**
     * @cond
     */
    virtual ~SequenceOracle() = default;
    /**
     * @endcond
     */

    /**
     * Predict the indices to be accessed in future `fetch()` calls.
     *
     * @param[out] predicted Pointer to an array in which to store the predicted indices of future elements to be accessed by `fetch()`.
     * @param number Maximum number of indices to predict.
     *
     * @return Number of indices that were predicted.
     * This is guaranteed to be no greater than `number`.
     */
    virtual size_t predict(Index_* predicted, size_t number) = 0;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a known sequence.
 */
template<typename Index_>
struct FixedOracle : public SequenceOracle<Index_> {
    /**
     * @param r Pointer to a constant array of indices.
     * @param n Length of the array for `r`.
     */
    FixedOracle(const Index_* r, size_t n) : reference(r), length(n) {}

    size_t predict(Index_* predicted, size_t number) {
        auto upto = counter + number;
        auto out = reference + counter;
        if (upto >= length) {
            number = length - counter;
            upto = length;
        }
        counter = upto;
        std::copy(out, out + number, predicted);
        return number;
    }
private:
    const Index_* reference;

    size_t length;

    size_t counter = 0;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses of a consecutive sequence.
 */
template<typename Index_>
struct ConsecutiveOracle : public SequenceOracle<Index_> {
    /**
     * @param s Start index of the consecutive sequence.
     * @param e One past the end of the sequence.
     */
    ConsecutiveOracle(Index_ s, Index_ e) : end(e), counter(s) {}

    size_t predict(Index_* buffer, size_t number) {
        auto upto = counter + number;
        if (upto >= end) {
            number = end - counter;
            upto = end;
        }
        std::iota(buffer, buffer + number, counter);
        counter = upto;
        return number;
    }
private:
    size_t end;

    size_t counter = 0;
};

}

#endif
