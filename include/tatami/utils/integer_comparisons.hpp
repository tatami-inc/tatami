#ifndef TATAMI_SAFE_NON_NEGATIVE_EQUAL_HPP
#define TATAMI_SAFE_NON_NEGATIVE_EQUAL_HPP

#include <type_traits>

/**
 * @file integer_comparisons.hpp
 * @brief Sign-aware integer comparisons.
 */

namespace tatami {

/**
 * This is effectively a copy of C++20's `std::cmp_equal()`, but we aren't guaranteed access to C++20 when compiling **tatami**.
 *
 * @tparam Left_ Integer type.
 * @tparam Right_ Another integer type.
 * @param l Integer value.
 * @param r Integer value.
 * @return Whether the two integer values are non-negative and equal.
 */
template<typename Left_, typename Right_>
bool safe_non_negative_equal(Left_ l, Right_ r) {
    return l >= 0 && r >= 0 && static_cast<typename std::make_unsigned<Left_>::type>(l) == static_cast<typename std::make_unsigned<Right_>::type>(r);
}

}

#endif
