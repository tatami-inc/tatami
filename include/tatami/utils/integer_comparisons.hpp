#ifndef TATAMI_SAFE_NON_NEGATIVE_EQUAL_HPP
#define TATAMI_SAFE_NON_NEGATIVE_EQUAL_HPP

#include <type_traits>
#include <limits>

#include "sanisizer/sanisizer.hpp"

namespace tatami {

/**
 * @cond
 */ 
template<typename Left_, typename Right_> // provided only for back-compatibility.
bool safe_non_negative_equal(Left_ l, Right_ r) {
    return l >= 0 && r >= 0 && sanisizer::is_equal(l, r);
}

template<typename Container_, typename Index_>
Index_ can_cast_Index_to_container_size(Index_ x) {
    // If the Index_ type is larger than size_t, we cast it to size_t to provide can_cast() with some opportunities for compile-time optimization.
    // This is because we know that all uses of Index_ must fit into a size_t in order for the Extractor::fetch calls to address the user-supplied arrays.
    // By casting away larger types, we allow can_cast() to eliminate the run-time overflow checks completely.
    typedef typename std::conditional<std::numeric_limits<std::size_t>::max() < std::numeric_limits<Index_>::max(), std::size_t, Index_>::type Intermediate;
    sanisizer::can_cast<decltype(std::declval<Container_>().size())>(static_cast<Intermediate>(x));
    return x;
}

template<typename Container_, typename Index_>
decltype(std::declval<Container_>().size()) cast_Index_to_container_size(Index_ x) {
    return can_cast_Index_to_container_size<Container_>(x);
}

template<typename Container_, typename Index_>
Container_ create_container_of_Index_size(Index_ x) {
    // Same logic as described above.
    typedef typename std::conditional<std::numeric_limits<std::size_t>::max() < std::numeric_limits<Index_>::max(), std::size_t, Index_>::type Intermediate;
    return sanisizer::create<Container_>(static_cast<Intermediate>(x));
}

template<typename Container_, typename Index_>
void resize_container_to_Index_size(Container_& container, Index_ x) {
    container.resize(cast_Index_to_container_size<Container_>(x));
}


/**
 * @endcond
 */ 

}

#endif
