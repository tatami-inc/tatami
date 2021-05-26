#ifndef IS_CONTIGUOUS_H
#define IS_CONTIGUOUS_H

#include <type_traits>
#include <utility>

namespace bioc {

// Default to false.
template<typename T, class V, typename = int>
struct has_data {
    static const bool value = false;
};

// Specialization is only run if it _has_ a data method.
template<typename T, class V>
struct has_data<T, V, decltype((void) V().data(), 0)> { 
    static const bool value = std::is_same<T*, decltype(V().data())>::value;
};

}

#endif
