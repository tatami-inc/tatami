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
 * @tparam Pointer_ Pointer type to the array of indices.
 *
 * @brief Predict future accesses from a view on a fixed sequence.
 */
template<typename Index_, typename Pointer_ = const Index_*>
class FixedViewOracle final : public Oracle<Index_> {
public:
    /**
     * @param ptr Pointer to a constant array of indices on the target dimension.
     * For non-smart pointers, the underlying array should be valid for the lifetime of this `FixedViewOracle` instance.
     * @param number Length of the array at `ptr`.
     */
    FixedViewOracle(Pointer_ ptr, const PredictionIndex number) : my_reference(ptr), my_length(number) {
        sanisizer::can_cast<std::size_t>(number); // make sure that array is addressable by all [0, number).
    }

    PredictionIndex total() const {
        return my_length;
    }

    Index_ get(PredictionIndex i) const {
        return my_reference[i];
    }

private:
    Pointer_ my_reference;
    PredictionIndex my_length;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam Pointer_ Container of indices.
 * This should support `size()` and access by `[]`.
 *
 * @brief Predict future accesses from a vector containing a fixed sequence.
 */
template<typename Index_, class Container_ = std::vector<Index_> >
class FixedVectorOracle final : public Oracle<Index_> {
public:
    /**
     * @param sequence Fixed sequence of indices on the target dimension.
     */
    FixedVectorOracle(Container_ sequence) : my_sequence(std::move(sequence)) {
        sanisizer::can_cast<PredictionIndex>(my_sequence.size()); // make sure that total() will return a sensible value.
    }

    PredictionIndex total() const {
        return my_sequence.size();
    }

    Index_ get(PredictionIndex i) const {
        return my_sequence[i];
    }

private:
    Container_ my_sequence;
};

}

#endif
