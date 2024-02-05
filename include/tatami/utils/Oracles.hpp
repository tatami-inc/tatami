#ifndef TATAMI_ORACLES_HPP
#define TATAMI_ORACLES_HPP

#include "../base/Options.hpp"
#include <numeric>

/**
 * @file Oracles.hpp
 *
 * @brief Predict future accesses during iteration.
 */

namespace tatami {

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a known sequence.
 */
template<typename Index_>
struct FixedOracle : public Oracle<Index_> {
    /**
     * @param r Pointer to a constant array of indices.
     * The underlying array should be valid for the lifetime of this `FixedOracle` instance.
     * @param n Length of the array at `r`.
     */
    FixedOracle(const Index_* r, size_t n) : reference(r), length(n) {}

    size_t total() const {
        return length;
    }

    Index_ get(size_t i) const {
        return reference[i];
    }

private:
    const Index_* reference;
    size_t length;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses of a consecutive sequence.
 */
template<typename Index_>
struct ConsecutiveOracle : public Oracle<Index_> {
    /**
     * @param s Start index of the consecutive sequence.
     * @param l Length of the sequence.
     */
    ConsecutiveOracle(Index_ s, Index_ l) : offset(s), length(l) {}

    size_t total() const {
        return length;
    }

    Index_ get(size_t i) const {
        return offset + i;
    }

private:
    Index_ offset, length;
};

}

#endif
