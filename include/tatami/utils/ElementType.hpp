#ifndef TATAMI_ELEMENTTYPE_HPP
#define TATAMI_ELEMENTTYPE_HPP

#include "copy.hpp"

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
using ElementType = I<decltype(std::declval<Array_>()[0])>;

}

#endif
