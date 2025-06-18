#ifndef TATAMI_CONSECUTIVE_ORACLE_HPP
#define TATAMI_CONSECUTIVE_ORACLE_HPP

#include "../base/Oracle.hpp"

#include "sanisizer/sanisizer.hpp"

/**
 * @file ConsecutiveOracle.hpp
 *
 * @brief Iterate across consecutive elements of the target dimension.
 */

namespace tatami {

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses along a consecutive sequence. 
 */
template<typename Index_>
class ConsecutiveOracle final : public Oracle<Index_> {
public:
    /**
     * @param start Start index of the consecutive sequence on the target dimension.
     * @param length Length of the sequence.
     */
    ConsecutiveOracle(Index_ start, Index_ length) : my_offset(start), my_length(sanisizer::cast<PredictionIndex>(length)) {}

    PredictionIndex total() const {
        return my_length;
    }

    Index_ get(PredictionIndex i) const {
        return my_offset + i;
    }

private:
    Index_ my_offset;
    PredictionIndex my_length;
};

}

#endif
