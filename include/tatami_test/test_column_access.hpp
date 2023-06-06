#ifndef TATAMI_TEST_COLUMN_ACCESS_H
#define TATAMI_TEST_COLUMN_ACCESS_H

#include "utils.hpp"
#include "test_access_base.hpp"
#include <type_traits>

namespace tatami_test {

template<bool has_nan_ = false, bool less_sparse_ = true, class Matrix_, class Matrix2_>
void test_simple_column_access(const Matrix_* ptr, const Matrix2_* ref, bool forward, int jump) {
    int NR = ref->nrow();

    auto rwork = ref->dense_column();
    test_access_base<false, has_nan_, less_sparse_>(
        ptr, 
        ref, 
        forward, 
        jump, 
        [&](int c) -> auto { 
            auto expected = rwork->fetch(c);
            EXPECT_EQ(expected.size(), NR);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, NR);
        }
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column();
    EXPECT_EQ(pwork->full_length, NR);

    auto swork = ptr->sparse_column();
    EXPECT_EQ(swork->full_length, NR);
}

template<bool has_nan_ = false, bool less_sparse_ = true, class Matrix_, class Matrix2_>
void test_sliced_column_access(const Matrix_* ptr, const Matrix2_* ref, bool forward, int jump, int start, int end) {
    auto rwork = ref->dense_column();
    test_access_base<false, has_nan_, less_sparse_>(
        ptr, 
        ref, 
        forward, 
        jump, 
        [&](int c) -> auto { 
            auto raw_expected = rwork->fetch(c);
            return std::vector<typename Matrix_::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& range) -> auto {
            return expand(range, start, end);
        },
        start,
        end - start
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column(start, end - start);
    EXPECT_EQ(pwork->block_start, start);
    EXPECT_EQ(pwork->block_length, end - start);

    auto swork = ptr->sparse_column(start, end - start);
    EXPECT_EQ(swork->block_start, start);
    EXPECT_EQ(swork->block_length, end - start);
}

template<bool has_nan_ = false, bool less_sparse_ = true, class Matrix_, class Matrix2_>
void test_indexed_column_access(const Matrix_* ptr, const Matrix2_* ref, bool forward, int jump, int start, int step) {
    int NR = ref->nrow();
    std::vector<typename Matrix_::index_type> indices;
    {
        int counter = start;
        while (counter < NR) {
            indices.push_back(counter);
            counter += step;
        }
    }

    // First trying with the index pointers.
    auto rwork = ref->dense_column();
    test_access_base<false, has_nan_, less_sparse_>(
        ptr, 
        ref, 
        forward, 
        jump, 
        [&](int c) -> auto { 
            auto raw_expected = rwork->fetch(c);
            std::vector<typename Matrix_::value_type> expected;
            expected.reserve(indices.size());
            for (auto idx : indices) {
                expected.push_back(raw_expected[idx]);
            }
            return expected;
        }, 
        [&](const auto& range) -> auto {
            auto full = expand(range, NR);
            std::vector<double> sub;
            sub.reserve(indices.size());
            for (auto idx : indices) {
                sub.push_back(full[idx]);
            }
            return sub;
        },
        indices
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column(indices);
    EXPECT_EQ(std::vector<int>(pwork->index_start(), pwork->index_start() + pwork->index_length), indices);

    auto swork = ptr->sparse_column(indices);
    EXPECT_EQ(std::vector<int>(swork->index_start(), swork->index_start() + swork->index_length), indices);
}

}

#endif
