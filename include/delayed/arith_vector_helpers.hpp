#ifndef DELAYED_ARITH_SCALAR_HELPER_H
#define DELAYED_ARITH_SCALAR_HELPER_H

#include <vector>

namespace bioc {

template<typename T, int MARGIN = 0, class V = std::vector<T> >
struct DelayedAddVectorHelper {
    DelayedAddVectorHelper(V&& v) : vec(v) {}
    DelayedAddVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            return val + vec[i]; 
        } else {
            return val + vec[c]; 
        }
    }
    static const bool sparse = false; 
private:
    const V vec;
};

template<typename T, bool RIGHT, int MARGIN = 0, class V = std::vector<T> >
struct DelayedSubtractVectorHelper {
    DelayedSubtractVectorHelper(V&& v) : vec(v) {}
    DelayedSubtractVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            if constexpr(RIGHT) {
                return val - vec[i];
            } else {
                return vec[i] - val;
            }
        } else {
            if constexpr(RIGHT) {
                return val - vec[c];
            } else {
                return vec[c] - val;
            }
        }
    }
    static const bool sparse = false; 
private:
    const V vec;
};

template<typename T, int MARGIN = 0, class V = std::vector<T> >
struct DelayedMultiplyVectorHelper {
    DelayedMultiplyVectorHelper(V&& v) : vec(v) {}
    DelayedMultiplyVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            return val * vec[i]; 
        } else {
            return val * vec[c]; 
        }
    }
    static const bool sparse = true;
private:
    const V vec;
};

template<typename T, bool RIGHT, int MARGIN = 0, class V = std::vector<T> >
struct DelayedDivideVectorHelper {
    DelayedDivideVectorHelper(V&& v) : vec(v) {}
    DelayedDivideVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, T val) const { 
        if constexpr(MARGIN==0) {
            if constexpr(RIGHT) {
                return val / vec[i];
            } else {
                return vec[i] / val;
            }
        } else {
            if constexpr(RIGHT) {
                return val / vec[c];
            } else {
                return vec[c] / val;
            }
        }
    }
    static const bool sparse = true;
private:
    const V vec;
};

}

#endif
