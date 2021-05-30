#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include <cmath>

namespace bioc {

template<typename T>
struct DelayedAbsHelper {
    DelayedAbsHelper() {};
    T operator()(size_t i, size_t c, T val) const {
        return std::abs(val);
    }
    static const bool sparse = true;
};

template<typename T>
struct DelayedLogHelper {
    DelayedLogHelper() : log_base(1) {}
    DelayedLogHelper(double base) : log_base(std::log(base)) {}
    T operator()(size_t i, size_t c, T val) const {
        return std::log(val)/log_base;
    }
    static const bool sparse = false;
private:
    const double log_base;
};

template<typename T>
struct DelayedSqrtHelper {
    T operator()(size_t i, size_t c, T val) const {
        return std::sqrt(val);
    }
    static const bool sparse = true;
};

template<typename T>
struct DelayedLog1pHelper {
    DelayedLog1pHelper() : log_base(1) {}
    DelayedLog1pHelper(double base) : log_base(std::log(base)) {}
    T operator()(size_t i, size_t c, T val) const {
        return std::log1p(val)/log_base;
    }
    static const bool sparse = true;
private:
    const double log_base;
};

template<typename T>
struct DelayedRoundHelper {
    T operator()(size_t i, size_t c, T val) const {
        return std::round(val);
    }
    static const bool sparse = true;
};

template<typename T>
struct DelayedExpHelper {
    T operator()(size_t i, size_t c, T val) const {
        return std::exp(val);
    }
    static const bool sparse = false;
};

}

#endif
