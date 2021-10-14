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

/******************/

template<class V, typename = int>
struct has_sparse_running {
    static constexpr bool value = false;
};

template<class V>
struct has_sparse_running<V, decltype((void) std::declval<V>().sparse_running(), 0)> {
    static constexpr bool value = true;
};

/******************/

template<class V, typename = int>
struct has_sparse_running_parallel {
    static constexpr bool value = false; 
};

template<class V>
struct has_sparse_running_parallel<V, decltype((void) std::declval<V>().sparse_running(0, 0), 0)> {
    static constexpr bool value = true;
};

/******************/

template<class V, typename = int>
struct has_dense_running {
    static constexpr bool value = false;
};

template<class V>
struct has_dense_running<V, decltype((void) std::declval<V>().dense_running(), 0)> {
    static constexpr bool value = true;
};

/******************/

template<class V, typename = int>
struct has_dense_running_parallel {
    static constexpr bool value = false;
};

template<class V>
struct has_dense_running_parallel<V, decltype((void) std::declval<V>().dense_running(0, 0), 0)> {
    static constexpr bool value = true;
};

/******************/

template<class V, typename T, typename = int>
struct has_nonconst_dense_compute {
    static constexpr bool value = false;
};

template<class V, typename T>
struct has_nonconst_dense_compute<V, T, decltype((void) std::declval<V>().compute_copy(0, std::declval<T*>()), 0)> { // note the const-ness.
    static constexpr bool value = true;
};

/******************/

template<class V, typename T, typename IDX, typename = int>
struct has_nonconst_sparse_compute {
    static constexpr bool value = false;
};

template<class V, typename T, typename IDX>
struct has_nonconst_sparse_compute<V, T, IDX, decltype((void) std::declval<V>().compute_copy(0, 0, std::declval<T*>(), std::declval<IDX*>()), 0)> {
    static constexpr bool value = true;
};

/******************/

template<class V, typename = int>
struct nonconst_sparse_compute_copy_mode {
    static constexpr SparseCopyMode value = SPARSE_COPY_BOTH;
};

template<class V>
struct nonconst_sparse_compute_copy_mode<V, decltype((void) V::copy_mode, 0)> {
    static constexpr SparseCopyMode value = V::copy_mode;
};

}

}

#endif 
