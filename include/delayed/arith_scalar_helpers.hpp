#ifndef DELAYED_ARITH_SCALAR_HELPER_H
#define DELAYED_ARITH_SCALAR_HELPER_H

#include "DelayedIsometricOp.hpp"

namespace bioc {

template<typename T, typename X = T, typename S = T>
struct DelayedAddScalarHelper {
    DelayedAddScalarHelper(S s) : scalar(s) {}
    T operator()(size_t i, size_t c, X val) const { 
        return val + scalar; 
    }
    bool sparse() const { return scalar == 0; } 
private:
    const S scalar;
};

template<typename T, typename X = T, typename S = T>
struct DelayedMultiplyScalarHelper { 
    DelayedMultiplyScalarHelper(S s) : scalar(s) {}
    T operator()(size_t i, size_t c, X val) const { 
        return val * scalar; 
    }
    bool sparse() const { return true; }
private:
    const S scalar;
};

template<typename T, bool RIGHT, typename X = T, typename S = T>
struct DelayedSubtractScalarHelper {
    DelayedSubtractScalarHelper(S s) : scalar(s) {}
    T operator()(size_t i, size_t c, X val) const { 
        if constexpr(RIGHT) {
            return val - scalar; 
        } else {
            return scalar - val;
        }
    }
    bool sparse() const { return scalar == 0; } 
private:
    const S scalar;
};

template<typename T, bool RIGHT, typename X = T, typename S = T>
struct DelayedDivideScalarHelper { 
    DelayedDivideScalarHelper(S s) : scalar(s) {}
    T operator()(size_t i, size_t c, X val) const { 
        if constexpr(RIGHT) {
            return val / scalar; 
        } else {
            return scalar / val;
        }
    }
    bool sparse() const { return true; } 
private:
    const S scalar;
};

}

#endif
