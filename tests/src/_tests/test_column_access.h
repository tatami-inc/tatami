#ifndef TEST_COLUMN_ACCESS_H
#define TEST_COLUMN_ACCESS_H
#include "utils.h"
#include <type_traits>

using tatami::DimensionLimit;
using tatami::ExtractOptions;

template<class Matrix, class Matrix2, class Function1, class Function2>
void test_simple_column_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, Function1 expector, Function2 sparse_expand, const DimensionLimit<int>& rlimits) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    ExtractOptions opt;
    if (!forward) {
        opt.access_order = tatami::AccessOrder::RANDOM;
    }

    DimensionLimit<int> climits;
    auto pwork = ptr->dense_column(rlimits, climits, opt);
    auto swork = ptr->sparse_column(rlimits, climits, opt);

    opt.sparse_ordered_index = false;
    auto swork_uns = ptr->sparse_column(rlimits, climits, opt);
    opt.sparse_ordered_index = true;

    opt.sparse_extract_index = false;
    auto swork_v = ptr->sparse_column(rlimits, climits, opt);
    opt.sparse_extract_value = false;
    auto swork_n = ptr->sparse_column(rlimits, climits, opt);
    opt.sparse_extract_index = true;
    auto swork_i = ptr->sparse_column(rlimits, climits, opt);
    opt.sparse_extract_value = true;

    opt.cache_for_reuse = true;
    auto cwork = ptr->dense_column(rlimits, climits, opt); 

    for (size_t i = 0; i < NC; i += jump) {
        size_t c = (forward ? i : NC - i - 1);

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
    for (size_t i = 0; i < NC; i += jump) {
        size_t c = (forward ? i : NC - i - 1);
        auto expected = pwork->fetch(c);
        auto observed = cwork->fetch(c);
        EXPECT_EQ(expected, observed);
    }
}

template<class Matrix, class Matrix2>
void test_simple_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump) {
    size_t NR = ref->nrow();

    auto rwork = ref->dense_column();
    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](size_t c) -> auto { 
            auto expected = rwork->fetch(c);
            EXPECT_EQ(expected.size(), NR);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, NR);
        },
        DimensionLimit<int>()
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column();
    EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::NONE);
    EXPECT_EQ(pwork->extracted_length, NR);

    auto swork = ptr->sparse_column();
    EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::NONE);
    EXPECT_EQ(swork->extracted_length, NR);
}

template<class Matrix, class Matrix2>
void test_sliced_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t end) {
    DimensionLimit<int> rsub;
    rsub.type = tatami::DimensionLimitType::BLOCK;
    rsub.block_start = start;
    rsub.block_length = end - start;

    auto rwork = ref->dense_column();
    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](size_t c) -> auto { 
            auto raw_expected = rwork->fetch(c);
            return std::vector<typename Matrix::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& range) -> auto {
            return expand(range, start, end);
        },
        rsub
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column(rsub, DimensionLimit<int>(), ExtractOptions());
    EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::BLOCK);
    EXPECT_EQ(pwork->extracted_block, start);
    EXPECT_EQ(pwork->extracted_length, end - start);

    auto swork = ptr->sparse_column(rsub, DimensionLimit<int>(), ExtractOptions());
    EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::BLOCK);
    EXPECT_EQ(swork->extracted_block, start);
    EXPECT_EQ(swork->extracted_length, end - start);
}

template<class Matrix, class Matrix2>
void test_indexed_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t step) {
    size_t NR = ref->nrow();
    std::vector<int> indices;
    {
        int counter = start;
        while (counter < NR) {
            indices.push_back(counter);
            counter += step;
        }
    }

    // First trying with the index pointers.
    {
        DimensionLimit<int> rsub;
        rsub.type = tatami::DimensionLimitType::INDEX;
        rsub.index_length = indices.size();
        rsub.index_start = indices.data();

        auto rwork = ref->dense_column();
        test_simple_column_access_base(ptr, ref, forward, jump, 
            [&](size_t c) -> auto { 
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
            rsub
        );

        // Checking that properties are correctly passed down.
        auto pwork = ptr->dense_column(rsub, DimensionLimit<int>(), ExtractOptions());
        EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(pwork->extracted_index(), indices.data());
        EXPECT_EQ(pwork->extracted_length, indices.size());

        auto swork = ptr->sparse_column(rsub, DimensionLimit<int>(), ExtractOptions());
        EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(swork->extracted_index(), indices.data());
        EXPECT_EQ(swork->extracted_length, indices.size());
    }

    // Now trying with full indices.
    {
        DimensionLimit<int> rsub;
        rsub.type = tatami::DimensionLimitType::INDEX;
        rsub.indices = indices;

        auto rwork = ref->dense_column();
        test_simple_column_access_base(ptr, ref, forward, jump, 
            [&](size_t c) -> auto { 
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
            rsub
        );

        // Checking that properties are correctly passed down.
        auto pwork = ptr->dense_column(rsub, DimensionLimit<int>(), ExtractOptions());
        EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(pwork->extracted_length, indices.size());
        EXPECT_EQ(std::vector<int>(pwork->extracted_index(), pwork->extracted_index() + pwork->extracted_length), indices);

        auto swork = ptr->sparse_column(rsub, DimensionLimit<int>(), ExtractOptions());
        EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(swork->extracted_length, indices.size());
        EXPECT_EQ(std::vector<int>(swork->extracted_index(), swork->extracted_index() + swork->extracted_length), indices);
    }}

#endif
