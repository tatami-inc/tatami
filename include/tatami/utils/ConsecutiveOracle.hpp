#ifndef TATAMI_CONSECUTIVE_ORACLE_HPP
#define TATAMI_CONSECUTIVE_ORACLE_HPP

#include "../base/Oracle.hpp"
#include <numeric>

/**
 * @file ConsecutiveOracle.hpp
 *
 * @brief Iterate across consecutive elements of a dimension.
 */

namespace tatami {

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
