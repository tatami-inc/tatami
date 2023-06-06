#ifndef TATAMI_TEST_ROW_ACCESS_HPP
#define TATAMI_TEST_ROW_ACCESS_HPP

#include "utils.hpp"
#include "test_access_base.hpp"
#include <type_traits>

namespace tatami_test {

template<bool has_nan_ = false, bool less_sparse_ = true, class Matrix_, class Matrix2_>
void test_simple_row_access(const Matrix_* ptr, const Matrix2_* ref, bool forward = true, int jump = 1) {
    int NC = ref->ncol();

    auto rwork = ref->dense_row();
    test_access_base<true, has_nan_, less_sparse_>(
        ptr, 
        ref, 
        forward, 
        jump, 
        [&](int r) -> auto { 
            auto expected = rwork->fetch(r);
            EXPECT_EQ(expected.size(), NC);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, NC);
        }
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_row();
    EXPECT_EQ(pwork->full_length, NC);

    auto swork = ptr->sparse_row();
    EXPECT_EQ(swork->full_length, NC);
}

template<bool has_nan_ = false, bool less_sparse_ = true, class Matrix_, class Matrix2_>
void test_sliced_row_access(const Matrix_* ptr, const Matrix2_* ref, bool forward, int jump, int start, int end) {
    auto rwork = ref->dense_row();
    test_access_base<true, has_nan_, less_sparse_>(
        ptr, 
        ref, 
        forward, 
        jump, 
        [&](int r) -> auto { 
            auto raw_expected = rwork->fetch(r);
            return std::vector<typename Matrix_::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& range) -> auto {
            return expand(range, start, end);
        },
        start,
        end - start
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_row(start, end - start);
    EXPECT_EQ(pwork->block_start, start);
    EXPECT_EQ(pwork->block_length, end - start);

    auto swork = ptr->sparse_row(start, end - start);
    EXPECT_EQ(swork->block_start, start);
    EXPECT_EQ(swork->block_length, end - start);
}

template<bool has_nan_ = false, bool less_sparse_ = true, class Matrix_, class Matrix2_>
void test_indexed_row_access(const Matrix_* ptr, const Matrix2_* ref, bool forward, int jump, int start, int step) {
    int NC = ptr->ncol();
    std::vector<typename Matrix_::index_type> indices;
    {
        int counter = start;
        while (counter < NC) {
            indices.push_back(counter);
            counter += step;
        }
    }

    auto rwork = ref->dense_row();
    test_access_base<true, has_nan_, less_sparse_>(
        ptr, 
        ref, 
        forward, 
        jump, 
        [&](int r) -> auto { 
            auto raw_expected = rwork->fetch(r);
            std::vector<typename Matrix_::value_type> expected;
            expected.reserve(indices.size());
            for (auto idx : indices) {
                expected.push_back(raw_expected[idx]);
            }
            return expected;
        },
        [&](const auto& range) -> auto {
            auto full = expand(range, NC);
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
    auto pwork = ptr->dense_row(indices);
    EXPECT_EQ(std::vector<int>(pwork->index_start(), pwork->index_start() + pwork->index_length), indices);

    auto swork = ptr->sparse_row(indices);
    EXPECT_EQ(std::vector<int>(swork->index_start(), swork->index_start() + swork->index_length), indices);
}

}

#endif
