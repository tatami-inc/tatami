#ifndef TATAMI_TEST_ISOMETRIC_UTILS_HPP
#define TATAMI_TEST_ISOMETRIC_UTILS_HPP

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

template<class Matrix1, class Matrix2>
void quick_test_all(const Matrix1* mat, const Matrix2* ref, bool has_nan = false) {
    tatami_test::TestAccessParameters params;
    params.has_nan = has_nan;
    int nrow = mat->nrow();
    int ncol = mat->ncol();

    params.use_row = true;
    tatami_test::test_full_access(params, mat, ref);
    tatami_test::test_block_access(params, mat, ref, ncol * 0.25, ncol * 0.9);
    tatami_test::test_indexed_access(params, mat, ref, ncol * 0.4, 11);

    params.use_row = false;
    tatami_test::test_full_access(params, mat, ref);
    tatami_test::test_block_access(params, mat, ref, nrow * 0.4, nrow * 0.9);
    tatami_test::test_indexed_access(params, mat, ref, nrow * 0.25, 10);
}

#endif
