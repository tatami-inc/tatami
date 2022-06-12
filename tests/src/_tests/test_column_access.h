#ifndef TEST_COLUMN_ACCESS_H
#define TEST_COLUMN_ACCESS_H
#include "utils.h"

template<class Matrix, class Matrix2>
void test_simple_column_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, size_t jump = 1) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto wrk = ptr->new_workspace(false);
    auto wrk_bi = ptr->new_workspace(false);

    for (size_t i = 0; i < NC; i += jump) {
        size_t c = (forward ? i : NC - i - 1);

        auto expected = ref->column(c);
        {
            auto observed = ptr->column(c);
            EXPECT_EQ(expected, observed);

            // Now with a workspace.
            auto observedW = ptr->column(c, wrk.get());
            EXPECT_EQ(expected, observedW);
        }

        {
            auto observed = ptr->sparse_column(c);
            EXPECT_EQ(expected, expand(observed, NR));

            // Now with a workspace.
            auto observedW = ptr->sparse_column(c, wrk.get());
            EXPECT_EQ(expected, expand(observedW, NR));
        }

        // Check workspace caching when access is bidirectional,
        // i.e., not purely increasing or decreasing.
        if (jump > 1 && i) {
            auto observed = ptr->column(c, wrk_bi.get());
            EXPECT_EQ(expected, observed);

            auto subc = (forward ? c-1 : c+1);
            auto observedm1 = ptr->column(subc, wrk_bi.get());
            auto expectedm1 = ref->column(subc);
            EXPECT_EQ(expectedm1, observedm1);
        }
    }
}

template<class Matrix, class Matrix2>
void test_sliced_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t first, size_t len, size_t shift) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto wrk = ptr->new_workspace(false);
    auto wrk_bi = ptr->new_workspace(false);

    for (size_t i = 0; i < NC; i += jump, first += shift) {
        size_t c = (forward ? i : NC - i - 1);
        auto interval = wrap_intervals(first, first + len, NR);
        size_t start = interval.first, end = interval.second;

        auto expected = ref->column(c, start, end);
        {
            auto observed = ptr->column(c, start, end);
            EXPECT_EQ(expected, observed);

            // Now with a workspace.
            auto observedW = ptr->column(c, start, end, wrk.get());
            EXPECT_EQ(expected, observedW);
        }

        {
            auto observed = ptr->sparse_column(c, start, end);
            EXPECT_EQ(expected, expand(observed, start, end));

            // Now with a workspace.
            auto observedW = ptr->sparse_column(c, start, end, wrk.get());
            EXPECT_EQ(expected, expand(observedW, start, end));
        }

        // Check workspace caching when access is bidirectional,
        // i.e., not purely increasing or decreasing.
        if (jump > 1 && i) {
            auto observed = ptr->column(c, start, end, wrk_bi.get());
            EXPECT_EQ(expected, observed);

            auto subc = (forward ? c-1 : c+1);
            auto observedm1 = ptr->column(subc, start, end, wrk_bi.get());
            auto expectedm1 = ref->column(subc, start, end);
            EXPECT_EQ(expectedm1, observedm1);
        }
    }
}

#endif
