#ifndef TATAMI_ARITH_SCALAR_HELPERS_H
#define TATAMI_ARITH_SCALAR_HELPERS_H

namespace tatami {

template<typename T>
struct DelayedAddScalarHelper {
    DelayedAddScalarHelper(T s) : scalar(s) {}
    T operator()(size_t i, size_t c, T val) const { 
        return val + scalar; 
    }
    static const bool sparse = false;
private:
    const T scalar;
};

template<typename T>
struct DelayedMultiplyScalarHelper { 
    DelayedMultiplyScalarHelper(T s) : scalar(s) {}
    T operator()(size_t i, size_t c, T val) const { 
        return val * scalar; 
    }
    static const bool sparse = true;
private:
    const T scalar;
};

template<typename T, bool RIGHT>
struct DelayedSubtractScalarHelper {
    DelayedSubtractScalarHelper(T s) : scalar(s) {}
    T operator()(size_t i, size_t c, T val) const { 
        if constexpr(RIGHT) {
            return val - scalar; 
        } else {
            return scalar - val;
        }
    }
    static const bool sparse = false;
private:
    const T scalar;
};

template<typename T, bool RIGHT>
struct DelayedDivideScalarHelper { 
    DelayedDivideScalarHelper(T s) : scalar(s) {}
    T operator()(size_t i, size_t c, T val) const { 
        if constexpr(RIGHT) {
            return val / scalar; 
        } else {
            return scalar / val;
        }
    }
    static const bool sparse = true;
private:
    const T scalar;
};

}

#endif
