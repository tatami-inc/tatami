#ifndef TATAMI_TEST_ISOMETRIC_UTILS_HPP
#define TATAMI_TEST_ISOMETRIC_UTILS_HPP

#include "tatami/base/Matrix.hpp"
#include "tatami_test/tatami_test.hpp"

inline double careful_division(double left, double right) {
    if (right == 0) {
        if (left > 0) {
            return std::numeric_limits<double>::infinity();
        } else if (left < 0) {
            return -std::numeric_limits<double>::infinity();
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    } else {
        return left / right;
    }
}

inline double careful_modulo(double left, double right) {
    if (right == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    } else {
        auto out = std::fmod(left, right);
        return out + (left/right < 0 && out != 0 ? right : 0);
    }
}

template<typename Value_, typename Index_>
void quick_test_all(const tatami::Matrix<Value_, Index_>& mat, const tatami::Matrix<Value_, Index_>& ref) {
    tatami_test::TestAccessOptions opts;

    opts.use_row = true;
    tatami_test::test_full_access(mat, ref, opts);
    tatami_test::test_block_access(mat, ref, 0.25, 0.6, opts);
    tatami_test::test_indexed_access(mat, ref, 0.4, 0.2, opts);

    opts.use_row = false;
    tatami_test::test_full_access(mat, ref, opts);
    tatami_test::test_block_access(mat, ref, 0.4, 0.5, opts);
    tatami_test::test_indexed_access(mat, ref, 0.25, 0.1, opts);
}

#endif
