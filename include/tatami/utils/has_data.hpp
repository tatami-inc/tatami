#ifndef TATAMI_HAS_DATA_HPP
#define TATAMI_HAS_DATA_HPP

#include <utility>
#include <type_traits>

namespace tatami {

template<typename Type_, class Container_, typename = int>
struct has_data {
    static const bool value = false;
};

template<typename Type_, class Container_>
struct has_data<Type_, Container_, decltype((void) std::declval<Container_>().data(), 0)> { 
    typedef decltype(std::declval<Container_>().data()) DataResult;

    static const bool value = std::is_same<Type_*, DataResult>::value || std::is_same<const Type_*, DataResult>::value;
};

}

#endif
