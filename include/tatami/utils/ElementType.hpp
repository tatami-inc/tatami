#ifndef TATAMI_ELEMENTTYPE_HPP
#define TATAMI_ELEMENTTYPE_HPP

/**
 * @file ElementType.hpp
 * @brief Get type of elements in an array.
 */

namespace tatami {

/**
 * @tparam Array_ Some array of values that are accessed with `[`.
 *
 * Extract the type of array elements, after stripping away references and `const`-ness.
 */
template<class Array_>
using ElementType = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Array_>()[0])>::type>::type;

}

#endif
