#ifndef TATAMI_FIXED_ORACLE_HPP
#define TATAMI_FIXED_ORACLE_HPP

#include "../base/Oracle.hpp"
#include <numeric>

/**
 * @file FixedOracle.hpp
 *
 * @brief Iterate across a fixed sequence of elements on the target dimension.
 */

namespace tatami {

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a view on a fixed sequence.
 */
template<typename Index_>
struct FixedViewOracle : public Oracle<Index_> {
    /**
     * @param r Pointer to a constant array of indices on the target dimension.
     * The underlying array should be valid for the lifetime of this `FixedOracle` instance.
     * @param n Length of the array at `r`.
     */
    FixedViewOracle(const Index_* r, size_t n) : reference(r), length(n) {}

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
 * @brief Predict future accesses from a vector containing a fixed sequence.
 */
template<typename Index_>
struct FixedVectorOracle : public Oracle<Index_> {
    /**
     * @param v Vector containing a fixed sequence of indices on the target dimension.
     */
    FixedVectorOracle(std::vector<Index_> v) : sequence(std::move(v)) {}

    size_t total() const {
        return sequence.size();
    }

    Index_ get(size_t i) const {
        return sequence[i];
    }

private:
    std::vector<Index_> sequence;
};

}

#endif
