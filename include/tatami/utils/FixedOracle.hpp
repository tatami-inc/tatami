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
class FixedViewOracle : public Oracle<Index_> {
public:
    /**
     * @param ptr Pointer to a constant array of indices on the target dimension.
     * The underlying array should be valid for the lifetime of this `FixedViewOracle` instance.
     * @param number Length of the array at `ptr`.
     */
    FixedViewOracle(const Index_* ptr, size_t number) : my_reference(ptr), my_length(number) {}

    size_t total() const {
        return my_length;
    }

    Index_ get(size_t i) const {
        return my_reference[i];
    }

private:
    const Index_* my_reference;
    size_t my_length;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a vector containing a fixed sequence.
 */
template<typename Index_>
class FixedVectorOracle : public Oracle<Index_> {
public:
    /**
     * @param vector Vector containing a fixed sequence of indices on the target dimension.
     */
    FixedVectorOracle(std::vector<Index_> vector) : my_sequence(std::move(vector)) {}

    size_t total() const {
        return my_sequence.size();
    }

    Index_ get(size_t i) const {
        return my_sequence[i];
    }

private:
    std::vector<Index_> my_sequence;
};

}

#endif
