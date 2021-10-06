#ifndef TATAMI_APPLY_CONFIG_HPP
#define TATAMI_APPLY_CONFIG_HPP

#include <utility>

namespace tatami {

namespace stats {

template<class V, typename = int>
struct has_sparse_direct {
    static constexpr bool value = false;
};

template<class V>
struct has_sparse_direct<V, decltype((void) std::declval<V>().sparse_direct(), 0)> {
    static constexpr bool value = true;
};

template<class V, typename = int>
struct has_sparse_running {
    static constexpr bool value = false;
};

template<class V>
struct has_sparse_running<V, decltype((void) std::declval<V>().sparse_running(), 0)> {
    static constexpr bool value = true;
};

template<class V, typename = int>
struct has_sparse_running_parallel {
    static constexpr bool value = false; 
};

template<class V>
struct has_sparse_running_parallel<V, decltype((void) std::declval<V>().sparse_running(0, 0), 0)> {
    static constexpr bool value = true;
};

template<class V, typename = int>
struct has_dense_running {
    static constexpr bool value = false;
};

template<class V>
struct has_dense_running<V, decltype((void) std::declval<V>().dense_running(), 0)> {
    static constexpr bool value = true;
};

template<class V, typename = int>
struct has_dense_running_parallel {
    static constexpr bool value = false;
};

template<class V>
struct has_dense_running_parallel<V, decltype((void) std::declval<V>().dense_running(0, 0), 0)> {
    static constexpr bool value = true;
};

}

}

#endif 
