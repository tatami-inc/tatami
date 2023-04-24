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
protected:
    /**
     * @cond
     */
    SequenceOracle() = default;
    virtual ~SequenceOracle() = default;
    /**
     * @endcond
     */

public:
    /**
     * Predict the indices to be accessed in future `fetch()` calls.
     *
     * @param n Maximum number of indices to predict.
     *
     * @return Pointer to an array of indices of the future elements to be accessed by `fetch()`.
     * The length of this array is also returned and is guaranteed to be less than `n`.
     */
    virtual std::pair<const Index_*, size_t> predict(size_t n) = 0;
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

    std::pair<const Index_*, size_t> predict(size_t n) {
        auto upto = counter + n;
        auto out = reference + counter;
        if (upto >= length) {
            n = length - counter;
            upto = length;
        }
        counter = upto;
        return std::make_pair(out, n);
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

    std::pair<const Index_*, size_t> predict(size_t n) {
        auto upto = counter + n;
        if (upto >= end) {
            n = end - counter;
            upto = end;
        }
        buffer.resize(n);
        std::iota(buffer.begin(), buffer.end(), counter);
        counter = upto;
        return std::make_pair(buffer.data(), buffer.size());
    }
private:
    size_t end;

    size_t counter = 0;

    std::vector<Index_> buffer;
};


}

#endif
