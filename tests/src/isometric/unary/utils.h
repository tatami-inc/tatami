#ifndef TATAMI_TEST_ISOMETRIC_UNARY_UTILS_HPP
#define TATAMI_TEST_ISOMETRIC_UNARY_UTILS_HPP

template<typename T>
T careful_division(T left, T right) {
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

template<class Matrix1, class Matrix2>
void test_nan_access(const Matrix1* ref, const Matrix2* mat) {
    auto refext = ref->dense_row();
    auto matext = mat->dense_row();
    int nrow = ref->nrow();
    EXPECT_EQ(nrow, mat->nrow());

    for (int r = 0; r < nrow; ++r) {
        auto refout = refext->fetch(r);
        auto matout = matext->fetch(r);

        constexpr double placeholder = 1234567890;
        for (auto& x : refout) {
            if (std::isnan(x)) {
                x = placeholder;
            }
        }

        for (auto& x : matout) {
            if (std::isnan(x)) {
                x = placeholder;
            }
        }
        EXPECT_EQ(refout, matout);
    }
}

template<class Matrix1, class Matrix2>
void quick_test_all(const Matrix1* mat, const Matrix2* ref) {
    // Full access.
    tatami_test::test_simple_column_access(mat, ref, true, 1);
    tatami_test::test_simple_row_access(mat, ref, true, 1);

    int nrow = mat->nrow();
    int ncol = mat->ncol();

    // Block access.
    tatami_test::test_sliced_column_access(mat, ref, true, 1, nrow * 0.25, nrow * 0.9);
    tatami_test::test_sliced_row_access(mat, ref, true, 1, ncol * 0.4, ncol * 0.9);

    // Indexed access.
    tatami_test::test_indexed_column_access(mat, ref, true, 1, nrow * 0.25, 10);
    tatami_test::test_indexed_row_access(mat, ref, true, 1, ncol * 0.4, 11);
}

#endif
