#ifndef TATAMI_HAS_DATA_HPP
#define TATAMI_HAS_DATA_HPP

#include <type_traits>

namespace tatami {

/**
 * @brief Compile time check for the `data()` method.
 * @tparam T Expected type of the data.
 * @tparam V Class to check for `data()`.
 *
 * This returns `false` by default.
 */
template<typename T, class V, typename = int>
struct has_data {
    static const bool value = false;
};

/**
 * @brief Compile time check for the `data()` method.
 * @tparam T Expected type of the data.
 * @tparam V Class to check for `data()`.
 *
 * This only returns `true` if `V` has a `data()` method that returns a pointer to an array of `T`'s.
 */
template<typename T, class V>
struct has_data<T, V, decltype((void) std::declval<V>().data(), 0)> { 
    static const bool value = std::is_same<T*, decltype(std::declval<V>().data())>::value;
};

}

#endif
