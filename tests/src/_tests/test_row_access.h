#ifndef TEST_ROW_ACCESS_H
#define TEST_ROW_ACCESS_H
#include "utils.h"

template<class Matrix, class Matrix2, class Function1, class Function2, typename ...Args>
void test_simple_row_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, int jump, Function1 expector, Function2 sparse_expand, Args... args) {
    int NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    int NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto pwork = ptr->dense_row(args...);
    auto swork = ptr->sparse_row(args...);

    tatami::Options<int> opt;
    opt.sparse.ordered_index = false;
    auto swork_uns = ptr->sparse_row(args..., opt);
    opt.sparse.ordered_index = true;

    opt.sparse.extract_index = false;
    auto swork_v = ptr->sparse_row(args..., opt);
    opt.sparse.extract_value = false;
    auto swork_n = ptr->sparse_row(args..., opt);
    opt.sparse.extract_index = true;
    auto swork_i = ptr->sparse_row(args..., opt);
    opt.sparse.extract_value = true;

    opt.access.cache_for_reuse = true;
    auto cwork = ptr->dense_row(args..., opt); 

    for (int i = 0; i < NR; i += jump) {
        int r = (forward ? i : NR - i - 1);

        auto expected = expector(r);
        {
            auto observed = pwork->fetch(r);
            EXPECT_EQ(expected, observed);
        }

        // Various flavors of sparse retrieval.
        {
            auto observed = swork->fetch(r);
            EXPECT_EQ(expected, sparse_expand(observed));
            EXPECT_TRUE(is_increasing(observed.index));

            std::vector<int> indices(expected.size()); // using the dense expected size as a proxy for the extraction length in block/indexed cases.
            auto observed_i = swork_i->fetch(r, NULL, indices.data());
            EXPECT_TRUE(observed_i.value == NULL);
            EXPECT_EQ(observed.index, std::vector<int>(observed_i.index, observed_i.index + observed_i.number));

            std::vector<double> values(expected.size());
            auto observed_v = swork_v->fetch(r, values.data(), NULL);
            EXPECT_TRUE(observed_v.index == NULL);
            EXPECT_EQ(observed.value, std::vector<double>(observed_v.value, observed_v.value + observed_v.number));

            auto observed_n = swork_n->fetch(r, NULL, NULL);
            EXPECT_TRUE(observed_n.value == NULL);
            EXPECT_TRUE(observed_n.index == NULL);
            EXPECT_EQ(observed.value.size(), observed_n.number);

            auto observed_uns = swork_uns->fetch(r);
            EXPECT_EQ(expected, sparse_expand(observed_uns));
        } 

        // Checks for caching.
        {
            auto observed = cwork->fetch(r);
            EXPECT_EQ(expected, observed);
        }
    }

    // Check that the caching gives the same results as uncached pass.
    for (int i = 0; i < NR; i += jump) {
        int r = (forward ? i : NR - i - 1);
        auto expected = pwork->fetch(r);
        auto observed = cwork->fetch(r);
        EXPECT_EQ(expected, observed);
    }
}

template<class Matrix, class Matrix2>
void test_simple_row_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, int jump = 1) {
    int NC = ref->ncol();

    auto rwork = ref->dense_row();
    test_simple_row_access_base(ptr, ref, forward, jump, 
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

template<class Matrix, class Matrix2>
void test_sliced_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, int jump, int start, int end) {
    auto rwork = ref->dense_row();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](int r) -> auto { 
            auto raw_expected = rwork->fetch(r);
            return std::vector<typename Matrix::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
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

template<class Matrix, class Matrix2>
void test_indexed_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, int jump, int start, int step) {
    int NC = ptr->ncol();
    std::vector<typename Matrix::index_type> indices;
    {
        int counter = start;
        while (counter < NC) {
            indices.push_back(counter);
            counter += step;
        }
    }

    auto rwork = ref->dense_row();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](int r) -> auto { 
            auto raw_expected = rwork->fetch(r);
            std::vector<typename Matrix::value_type> expected;
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

#endif
