#ifndef TEST_ROW_ACCESS_H
#define TEST_ROW_ACCESS_H
#include "utils.h"

template<class Matrix, class Matrix2>
void test_simple_row_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, size_t jump = 1) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto wrk = ptr->new_workspace(true);
    for (size_t i = 0; i < NR; i += jump) {
        size_t c = (forward ? i : NR - i - 1);

        auto expected = ref->row(i);
        {
            auto observed = ptr->row(i);
            EXPECT_EQ(expected, observed);

            // Now with a workspace.
            auto observedW = ptr->row(i, wrk.get());
            EXPECT_EQ(expected, observedW);
        }

        {
            auto observed = ptr->sparse_row(i);
            EXPECT_EQ(expected, expand(observed, NC));

            // Now with a workspace.
            auto observedW = ptr->sparse_row(i, wrk.get());
            EXPECT_EQ(expected, expand(observedW, NC));
        }
    }
}

template<class Matrix, class Matrix2>
void test_sliced_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t first, size_t len, size_t shift) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto wrk = ptr->new_workspace(true);
    for (size_t i = 0; i < NR; i += jump, first += shift) {
        size_t c = (forward ? i : NR - i - 1);
        auto interval = wrap_intervals(first, first + len, NC);
        size_t start = interval.first, end = interval.second;
        
        auto expected = ref->row(i, start, end);
        {
            auto observed = ptr->row(i, start, end);
            EXPECT_EQ(expected, observed);

            // Now with a workspace.
            auto observedW = ptr->row(i, start, end, wrk.get());
            EXPECT_EQ(expected, observedW);
        }

        {
            auto observed = ptr->sparse_row(i, start, end);
            EXPECT_EQ(expected, expand(observed, start, end));

            // Now with a workspace.
            auto observedW = ptr->sparse_row(i, start, end, wrk.get());
            EXPECT_EQ(expected, expand(observedW, start, end));
        }
    }
}

#endif
