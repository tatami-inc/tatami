#ifndef TATAMI_FIXED_ORACLE_HPP
#define TATAMI_FIXED_ORACLE_HPP

#include "../base/Oracle.hpp"

#include <vector>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

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
class FixedViewOracle final : public Oracle<Index_> {
public:
    /**
     * @param ptr Pointer to a constant array of indices on the target dimension.
     * The underlying array should be valid for the lifetime of this `FixedViewOracle` instance.
     * @param number Length of the array at `ptr`.
     */
    FixedViewOracle(const Index_* ptr, PredictionIndex number) : my_reference(ptr), my_length(number) {
        sanisizer::can_cast<std::size_t>(number); // make sure that array is addressable by all [0, number).
    }

    PredictionIndex total() const {
        return my_length;
    }

    Index_ get(PredictionIndex i) const {
        return my_reference[i];
    }

private:
    const Index_* my_reference;
    PredictionIndex my_length;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a vector containing a fixed sequence.
 */
template<typename Index_>
class FixedVectorOracle final : public Oracle<Index_> {
public:
    /**
     * @param vector Vector containing a fixed sequence of indices on the target dimension.
     */
    FixedVectorOracle(std::vector<Index_> vector) : my_sequence(std::move(vector)) {
        sanisizer::can_cast<PredictionIndex>(my_sequence.size()); // make sure that total() will return a sensible value.
    }

    PredictionIndex total() const {
        return my_sequence.size();
    }

    Index_ get(PredictionIndex i) const {
        return my_sequence[i];
    }

private:
    std::vector<Index_> my_sequence;
};

}

#endif
