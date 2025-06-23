#ifndef TATAMI_INDEX_TO_CONTAINER_HPP
#define TATAMI_INDEX_TO_CONTAINER_HPP

#include <type_traits>
#include <limits>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

namespace tatami {

/**
 * @cond
 */
template<typename Left_, typename Right_> // provided only for back-compatibility.
bool safe_non_negative_equal(Left_ l, Right_ r) {
    return l >= 0 && r >= 0 && sanisizer::is_equal(l, r);
}
/**
 * @endcond
 */

/**
 * @tparam Container_ Container with a `size()` method.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param x Size of the container as an `Index_`, typically the dimension extents.
 * It is assumed that `x` can be represented by a `std::size_t`, even if `Index_` is a larger type.
 *
 * @return `x` as its input type.
 * An error is raised if casting to the type of `Container_::size()` would result in overflow.
 */
template<typename Container_, typename Index_>
Index_ can_cast_Index_to_container_size(Index_ x) {
    // If the Index_ type is larger than size_t, we cast it to size_t to provide can_cast() with some opportunities for compile-time optimization.
    // This is because we know that all uses of Index_ must fit into a size_t in order for the Extractor::fetch calls to address the user-supplied arrays.
    // By casting away larger types, we allow can_cast() to eliminate the run-time overflow checks completely.
    typedef typename std::conditional<std::numeric_limits<std::size_t>::max() < std::numeric_limits<Index_>::max(), std::size_t, Index_>::type Intermediate;
    sanisizer::can_cast<decltype(std::declval<Container_>().size())>(static_cast<Intermediate>(x));
    return x;
}

/**
 * @tparam Container_ Container with a `size()` method.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param x Size of the container as an `Index_`, typically the dimension extents.
 * It is assumed that `x` can be represented by a `std::size_t`, even if `Index_` is a larger type.
 *
 * @return `x` as the type of `Container_::size()`.
 * An error is raised if the cast would result in overflow.
 */
template<typename Container_, typename Index_>
decltype(std::declval<Container_>().size()) cast_Index_to_container_size(Index_ x) {
    return can_cast_Index_to_container_size<Container_>(x);
}

/**
 * @tparam Container_ Container with a `size()` method.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param x Size of the container as an `Index_`, typically the dimension extents.
 * It is assumed that `x` can be represented by a `std::size_t`, even if `Index_` is a larger type.
 *
 * @return Instance of a `Container_` of size `x`.
 * An error is raised if `x` is too large.
 */
template<typename Container_, typename Index_>
Container_ create_container_of_Index_size(Index_ x) {
    // Same logic as described above.
    typedef typename std::conditional<std::numeric_limits<std::size_t>::max() < std::numeric_limits<Index_>::max(), std::size_t, Index_>::type Intermediate;
    return sanisizer::create<Container_>(static_cast<Intermediate>(x));
}

/**
 * @tparam Container_ Container with a `size()` method.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param container_ Instance of a container.
 * On output, this is resized to size `x`.
 * An error is raised if `x` is too large.
 * @param x Size of the container as an `Index_`, typically the dimension extents.
 * It is assumed that `x` can be represented by a `std::size_t`, even if `Index_` is a larger type.
 */
template<typename Container_, typename Index_>
void resize_container_to_Index_size(Container_& container, Index_ x) {
    container.resize(cast_Index_to_container_size<Container_>(x));
}

}

#endif
