#ifndef TEST_COLUMN_ACCESS_H
#define TEST_COLUMN_ACCESS_H
#include "utils.h"
#include <type_traits>

template<class Matrix, class Matrix2, class Function1, class Function2, typename ...Args>
void test_simple_column_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, int jump, Function1 expector, Function2 sparse_expand, Args... args) {
    int NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    int NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto pwork = ptr->dense_column(args...);
    auto swork = ptr->sparse_column(args...);

    tatami::Options<int> opt;
    opt.sparse.ordered_index = false;
    auto swork_uns = ptr->sparse_column(args..., opt);
    opt.sparse.ordered_index = true;

    opt.sparse.extract_index = false;
    auto swork_v = ptr->sparse_column(args..., opt);
    opt.sparse.extract_value = false;
    auto swork_n = ptr->sparse_column(args..., opt);
    opt.sparse.extract_index = true;
    auto swork_i = ptr->sparse_column(args..., opt);
    opt.sparse.extract_value = true;

    opt.access.cache_for_reuse = true;
    auto cwork = ptr->dense_column(args..., opt); 

    for (int i = 0; i < NC; i += jump) {
        int c = (forward ? i : NC - i - 1);

        auto expected = expector(c);
        {
            auto observed = pwork->fetch(c);
            EXPECT_EQ(expected, observed);
        }

        // Various flavors of sparse retrieval.
        {
            auto observed = swork->fetch(c);
            EXPECT_EQ(expected, sparse_expand(observed));
            EXPECT_TRUE(is_increasing(observed.index));

            std::vector<int> indices(expected.size()); // using the dense expected size as a proxy for the extraction length in block/indexed cases.
            auto observed_i = swork_i->fetch(c, NULL, indices.data());
            EXPECT_TRUE(observed_i.value == NULL);
            EXPECT_EQ(observed.index, std::vector<int>(observed_i.index, observed_i.index + observed_i.number));

            std::vector<double> values(expected.size());
            auto observed_v = swork_v->fetch(c, values.data(), NULL);
            EXPECT_TRUE(observed_v.index == NULL);
            EXPECT_EQ(observed.value, std::vector<double>(observed_v.value, observed_v.value + observed_v.number));

            auto observed_n = swork_n->fetch(c, NULL, NULL);
            EXPECT_TRUE(observed_n.value == NULL);
            EXPECT_TRUE(observed_n.index == NULL);
            EXPECT_EQ(observed.value.size(), observed_n.number);

            auto observed_uns = swork_uns->fetch(c);
            EXPECT_EQ(expected, sparse_expand(observed_uns));
        } 

        // Checks for caching.
        {
            auto observed = cwork->fetch(c);
            EXPECT_EQ(expected, observed);
        }
    }

    // Check that the caching gives the same results as uncached pass.
    for (int i = 0; i < NC; i += jump) {
        int c = (forward ? i : NC - i - 1);
        auto expected = pwork->fetch(c);
        auto observed = cwork->fetch(c);
        EXPECT_EQ(expected, observed);
    }
}

template<class Matrix, class Matrix2>
void test_simple_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, int jump) {
    int NR = ref->nrow();

    auto rwork = ref->dense_column();
    test_simple_column_access_base(ptr, ref, forward, jump, 
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

template<class Matrix, class Matrix2>
void test_sliced_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, int jump, int start, int end) {
    auto rwork = ref->dense_column();
    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](int c) -> auto { 
            auto raw_expected = rwork->fetch(c);
            return std::vector<typename Matrix::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
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

template<class Matrix, class Matrix2>
void test_indexed_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, int jump, int start, int step) {
    int NR = ref->nrow();
    std::vector<typename Matrix::index_type> indices;
    {
        int counter = start;
        while (counter < NR) {
            indices.push_back(counter);
            counter += step;
        }
    }

    // First trying with the index pointers.
    auto rwork = ref->dense_column();
    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](int c) -> auto { 
            auto raw_expected = rwork->fetch(c);
            std::vector<typename Matrix::value_type> expected;
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
        indices.data(),
        indices.size()
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column(indices.data(), indices.size());
    EXPECT_EQ(std::vector<int>(pwork->index_start(), pwork->index_start() + pwork->index_length), indices);

    auto swork = ptr->sparse_column(indices.data(), indices.size());
    EXPECT_EQ(std::vector<int>(swork->index_start(), swork->index_start() + swork->index_length), indices);
}

#endif
