#ifndef TATAMI_OPTIONS_HPP
#define TATAMI_OPTIONS_HPP

#include <vector>
#include <memory>

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
 * This allows implementations to pre-fetch data for future requests to `DenseExtractor::fetch()` or `SparseExtractor::fetch()`.
 */
template<typename Index_>
struct SequenceOracle {
    /**
     * Predict the indices to be accessed in future `fetch()` calls.
     *
     * @param n Maximum number of indices to predict.
     *
     * @return Pointer to an array of indices of the future elements to be accessed by `fetch()`.
     * The length of this array is also returned and is guaranteed to be less than `n`.
     */
    std::pair<const Index_*, size_t> predict(size_t n) = 0;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a known sequence.
 *
 * Once the sequence is fully iterated, this instance will wrap around to the start of the sequence.
 */
template<typename Index_>
struct FixedOracle : public SequenceOracle {
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
            counter = 0;
        } else {
            counter = upto;
        }
        return std::make_pair(out, n);
    }
private:
    const Index_* reference;

    size_t length;

    size_t counter = 0;
};

}

#endif
