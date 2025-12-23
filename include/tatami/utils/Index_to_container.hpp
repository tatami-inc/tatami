#ifndef TATAMI_INDEX_TO_CONTAINER_HPP
#define TATAMI_INDEX_TO_CONTAINER_HPP

#include <type_traits>
#include <limits>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

#include "copy.hpp"

/**
 * @file Index_to_container.hpp
 * @brief Convert index type to container size.
 */

namespace tatami {

/**
 * @cond
 */
template<typename Left_, typename Right_> // provided only for back-compatibility.
bool safe_non_negative_equal(const Left_ l, const Right_ r) {
    return l >= 0 && r >= 0 && sanisizer::is_equal(l, r);
}

// By the time we get any kind of Index_ in a tatami context, it should represent a dimension extent, so we can assume that 'x' is non-negative. 
// We also know that all uses of Index_ must fit into a size_t in order for the Extractor::fetch calls to address the user-supplied arrays.
// By casting away larger types, we allow sanisizer functions to eliminate the run-time overflow checks completely.
template<typename Index_> 
auto attest_for_Index(Index_ x) {
    return sanisizer::attest_gez(sanisizer::attest_max_by_type<std::size_t>(x));
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
Index_ can_cast_Index_to_container_size(const Index_ x) {
    sanisizer::can_cast<I<decltype(std::declval<Container_>().size())> >(attest_for_Index(x));
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
I<decltype(std::declval<Container_>().size())> cast_Index_to_container_size(const Index_ x) {
    return can_cast_Index_to_container_size<Container_>(x);
}

/**
 * @tparam Container_ Container with a `size()` method and a constructor that accepts the size as the first argument.
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam Args_ Further arguments to pass to the `Container_`'s constructor.
 *
 * @param x Size of the container as an `Index_`, typically the dimension extents.
 * It is assumed that `x` can be represented by a `std::size_t`, even if `Index_` is a larger type.
 * @param args Further arguments to pass to the `Container_` constructor, after `x`.
 *
 * @return Instance of a `Container_` of size `x`.
 * An error is raised if `x` is too large.
 */
template<typename Container_, typename Index_, typename ... Args_>
Container_ create_container_of_Index_size(const Index_ x, Args_&& ... args) {
    return sanisizer::create<Container_>(attest_for_Index(x), std::forward<Args_>(args)...);
}

/**
 * @tparam Container_ Container with a `size()` method and a `resize()` method that accepts the new length as the first argument.
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam Args_ Further arguments to pass to `Container_::resize()`.
 *
 * @param container Instance of a container.
 * On output, this is resized to size `x`.
 * An error is raised if `x` is too large.
 * @param x Size of the container as an `Index_`, typically the dimension extents.
 * It is assumed that `x` can be represented by a `std::size_t`, even if `Index_` is a larger type.
 * @param args Further arguments to pass to the `resize()` method, after `x`.
 */
template<typename Container_, typename Index_, typename ... Args_>
void resize_container_to_Index_size(Container_& container, const Index_ x, Args_&& ... args) {
    container.resize(cast_Index_to_container_size<Container_>(x), std::forward<Args_>(args)...);
}

}

#endif
