#ifndef TATAMI_HAS_DATA_HPP
#define TATAMI_HAS_DATA_HPP

#include <utility>
#include <type_traits>

/**
 * @file has_data.hpp
 * @brief Compile-time checks for the `data()` method.
 */

namespace tatami {

/**
 * @brief Compile time check for the `data()` method.
 * @tparam Type_ Expected type of the data.
 * @tparam Container_ Class to check for `data()`.
 *
 * This returns `false` by default.
 */
template<typename Type_, class Container_, typename = int>
struct has_data {
    /**
     * Compile-time constant indicating whether `data()` exists.
     */
    static const bool value = false;
};

/**
 * @brief Compile time check for the `data()` method.
 * @tparam Type_ Expected type of the data.
 * @tparam Container_ Class to check for `data()`.
 *
 * This only returns `true` if `Container_` has a `data()` method that returns a (possibly `const`) pointer to an array of `Type_`'s.
 */
template<typename Type_, class Container_>
struct has_data<Type_, Container_, decltype((void) std::declval<Container_>().data(), 0)> { 
    /**
     * @cond
     */
    typedef decltype(std::declval<Container_>().data()) DataResult;
    /**
     * @endcond
     */

    /**
     * Compile-time constant indicating whether `data()` exists.
     */
    static const bool value = std::is_same<Type_*, DataResult>::value || std::is_same<const Type_*, DataResult>::value;
};

}

#endif
