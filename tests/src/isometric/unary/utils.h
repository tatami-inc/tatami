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

#endif
