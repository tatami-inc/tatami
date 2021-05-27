#ifndef DELAYED_ARITH_SCALAR_HELPER_H
#define DELAYED_ARITH_SCALAR_HELPER_H

#include "DelayedIsometricOp.hpp"

namespace bioc {

template<typename T, class V, int MARGIN = 1, typename X = T>
struct DelayedAddVectorHelper {
    DelayedAddVectorHelper(V&& v) : vec(v) {}
    DelayedAddVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, X val) const { 
        if constexpr(MARGIN==1) {
            return val + vec[i]; 
        } else {
            return val + vec[c]; 
        }
    }
    bool sparse() const { return false; }
private:
    const V vec;
};

template<typename T, bool RIGHT, class V, int MARGIN = 1, typename X = T>
struct DelayedSubtractVectorHelper {
    DelayedSubtractVectorHelper(V&& v) : vec(v) {}
    DelayedSubtractVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, X val) const { 
        if constexpr(MARGIN==1) {
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
    bool sparse() const { return false; }
private:
    const V vec;
};

template<typename T, class V, int MARGIN = 1, typename X = T>
struct DelayedMultiplyVectorHelper {
    DelayedMultiplyVectorHelper(V&& v) : vec(v) {}
    DelayedMultiplyVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, X val) const { 
        if constexpr(MARGIN==1) {
            return val * vec[i]; 
        } else {
            return val * vec[c]; 
        }
    }
    bool sparse() const { return true; }
private:
    const V vec;
};

template<typename T, bool RIGHT, class V, int MARGIN = 1, typename X = T>
struct DelayedDivideVectorHelper {
    DelayedDivideVectorHelper(V&& v) : vec(v) {}
    DelayedDivideVectorHelper(const V& v) : vec(v) {}
    T operator()(size_t i, size_t c, X val) const { 
        if constexpr(MARGIN==1) {
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
    bool sparse() const { return true; }
private:
    const V vec;
};

}

#endif
