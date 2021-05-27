#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include "DelayedIsometricOp.hpp"
#include <cmath>

namespace bioc {

template<typename T, typename X = T>
struct DelayedLogHelper {
    DelayedLogHelper() : log_base(1) {}
    DelayedLogHelper(double base) : log_base(std::log(base)) {}
    T operator()(size_t i, size_t c, X val) const {
        return std::log(val)/log_base;
    }
    bool sparse() const { return false; }
private:
    const double log_base;
};

template<typename T, typename X = T>
struct DelayedSqrtHelper {
    T operator()(size_t i, size_t c, X val) const {
        return std::sqrt(val);
    }
    bool sparse() const { return true; }
};

template<typename T, typename X = T>
struct DelayedLog1pHelper {
    DelayedLogHelper() : log_base(1) {}
    DelayedLogHelper(double base) : log_base(std::log(base)) {}
    T operator()(size_t i, size_t c, X val) const {
        return std::log1p(val)/log_base;
    }
    bool sparse() const { return true; }
private:
    const double log_base;
};

template<typename T, typename X = T>
struct DelayedRoundHelper {
    T operator()(size_t i, size_t c, X val) const {
        return std::round(val);
    }
    bool sparse() const { return true; }
};

template<typename T, typename X = T>
struct DelayedExpHelper {
    T operator()(size_t i, size_t c, X val) const {
        return std::exp(val);
    }
    bool sparse() const { return false; }
};

}

#endif
