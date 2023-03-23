#ifndef TEST_ROW_ACCESS_H
#define TEST_ROW_ACCESS_H
#include "utils.h"

template<class Matrix, class Matrix2>
void test_simple_row_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, size_t jump = 1) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto rwork = ref->new_row_workspace();
    auto pwork = ptr->new_row_workspace();
    auto swork = ptr->new_row_workspace();
    auto pwork_bi = ptr->new_row_workspace();
    auto rwork_bi = ref->new_row_workspace();

    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);

        auto expected = ref->row(r, rwork.get());
        EXPECT_EQ(expected.size(), NC);
        {
            auto observed = ptr->row(r, pwork.get());
            EXPECT_EQ(expected, observed);
        }
        {
            auto observed = ptr->sparse_row(r, swork.get());
            EXPECT_EQ(expected, expand(observed, NC));
        }

        // Check workspace caching when access is bidirectional,
        // i.e., not purely increasing or decreasing.
        if (jump > 1 && i) {
            auto observed = ptr->row(r, pwork_bi.get());
            EXPECT_EQ(expected, observed);

            auto subr = (forward ? r-1 : r+1);
            auto observedm1 = ptr->row(subr, pwork_bi.get());
            auto expectedm1 = ref->row(subr, rwork_bi.get());
            EXPECT_EQ(expectedm1, observedm1);
        }
    }
}

template<class Matrix, class Matrix2>
void test_sliced_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t end) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto rwork = ref->new_row_workspace();
    auto pwork = ptr->new_row_workspace(start, end - start);
    auto swork = ptr->new_row_workspace(start, end - start);
    auto pwork_bi = ptr->new_row_workspace(start, end - start);
    auto rwork_bi = ref->new_row_workspace(start, end - start);

    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);

        auto raw_expected = ref->row(r, rwork.get());
        std::vector<typename Matrix::data_type> expected(raw_expected.begin() + start, raw_expected.begin() + end);
        {
            auto observed = ptr->row(r, pwork.get());
            EXPECT_EQ(expected, observed);
        }
        {
            auto observed = ptr->sparse_row(r, swork.get());
            EXPECT_EQ(expected, expand(observed, start, end));
        }

        // Check workspace caching when access is bidirectional,
        // i.e., not purely increasing or decreasing.
        if (jump > 1 && i) {
            auto observed = ptr->row(r, pwork_bi.get());
            EXPECT_EQ(expected, observed);

            auto subr = (forward ? r-1 : r+1);
            auto observedm1 = ptr->row(subr, pwork_bi.get());
            auto expectedm1 = ptr->row(subr, rwork_bi.get());
            EXPECT_EQ(expectedm1, observedm1);
        }
    }
}

template<class Matrix, class Matrix2>
void test_indexed_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t step) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    std::vector<typename Matrix::index_type> indices;
    {
        int counter = start;
        while (counter < NC) {
            indices.push_back(counter);
            counter += step;
        }
    }

    auto rwork = ref->new_row_workspace();
    auto pwork = ptr->new_row_workspace(indices.size(), indices.data());
    auto swork = ptr->new_row_workspace(indices.size(), indices.data());
    auto pwork_bi = ptr->new_row_workspace(indices.size(), indices.data());
    auto rwork_bi = ref->new_row_workspace(indices.size(), indices.data());

    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);

        auto raw_expected = ref->row(r, rwork.get());
        std::vector<typename Matrix::data_type> expected;
        expected.reserve(indices.size());
        for (auto idx : indices) {
            expected.push_back(raw_expected[idx]);
        }

        {
            auto observed = ptr->row(r, pwork.get());
            EXPECT_EQ(expected, observed);
        }
        {
            auto observed = ptr->sparse_row(r, swork.get());
            auto full = expand(observed, NC);
            std::vector<double> sub;
            sub.reserve(indices.size());
            for (auto idx : indices) {
                sub.push_back(full[idx]);
            }
            EXPECT_EQ(expected, sub);
        }

        // Check workspace caching when access is bidirectional,
        // i.e., not purely increasing or decreasing.
        if (jump > 1 && i) {
            auto observed = ptr->row(r, pwork_bi.get());
            EXPECT_EQ(expected, observed);

            auto subr = (forward ? r-1 : r+1);
            auto observedm1 = ptr->row(subr, pwork_bi.get());
            auto expectedm1 = ptr->row(subr, rwork_bi.get());
            EXPECT_EQ(expectedm1, observedm1);
        }
    }
}

#endif
