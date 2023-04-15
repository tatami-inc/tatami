#ifndef TEST_ROW_ACCESS_H
#define TEST_ROW_ACCESS_H
#include "utils.h"

using tatami::DimensionLimit;
using tatami::ExtractOptions;

template<class Matrix, class Matrix2, class Function1, class Function2, typename ...Args>
void test_simple_row_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, Function1 expector, Function2 sparse_expand, const DimensionLimit<int>& climits) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    ExtractOptions opt;
    if (!forward) {
        opt.access_order = tatami::AccessOrder::RANDOM;
    }

    DimensionLimit<int> rlimits;
    auto pwork = ptr->dense_row(rlimits, climits, opt);
    auto swork = ptr->sparse_row(rlimits, climits, opt);

    opt.sparse_ordered_index = false;
    auto swork_uns = ptr->sparse_row(rlimits, climits, opt);
    opt.sparse_ordered_index = true;

    opt.sparse_extract_index = false;
    auto swork_v = ptr->sparse_row(rlimits, climits, opt);
    opt.sparse_extract_value = false;
    auto swork_n = ptr->sparse_row(rlimits, climits, opt);
    opt.sparse_extract_index = true;
    auto swork_i = ptr->sparse_row(rlimits, climits, opt);
    opt.sparse_extract_value = true;

    opt.cache_for_reuse = true;
    auto cwork = ptr->dense_row(rlimits, climits, opt); 

    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);

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
    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);
        auto expected = pwork->fetch(r);
        auto observed = cwork->fetch(r);
        EXPECT_EQ(expected, observed);
    }
}

template<class Matrix, class Matrix2>
void test_simple_row_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, size_t jump = 1) {
    size_t NC = ref->ncol();

    auto rwork = ref->dense_row();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](size_t r) -> auto { 
            auto expected = rwork->fetch(r);
            EXPECT_EQ(expected.size(), NC);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, NC);
        },
        DimensionLimit<int>()
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_row();
    EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::NONE);
    EXPECT_EQ(pwork->extracted_length, NC);

    auto swork = ptr->sparse_row();
    EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::NONE);
    EXPECT_EQ(swork->extracted_length, NC);
}

template<class Matrix, class Matrix2>
void test_sliced_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t end) {
    DimensionLimit<int> csub;
    csub.type = tatami::DimensionLimitType::BLOCK;
    csub.block_start = start;
    csub.block_length = end - start;

    auto rwork = ref->dense_row();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](size_t r) -> auto { 
            auto raw_expected = rwork->fetch(r);
            return std::vector<typename Matrix::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& range) -> auto {
            return expand(range, start, end);
        },
        csub
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_row(DimensionLimit<int>(), csub, ExtractOptions());
    EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::BLOCK);
    EXPECT_EQ(pwork->extracted_block, start);
    EXPECT_EQ(pwork->extracted_length, end - start);

    auto swork = ptr->sparse_row(DimensionLimit<int>(), csub, ExtractOptions());
    EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::BLOCK);
    EXPECT_EQ(swork->extracted_block, start);
    EXPECT_EQ(swork->extracted_length, end - start);
}

template<class Matrix, class Matrix2>
void test_indexed_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t step) {
    size_t NC = ptr->ncol();
    std::vector<typename Matrix::index_type> indices;
    {
        int counter = start;
        while (counter < NC) {
            indices.push_back(counter);
            counter += step;
        }
    }

    // First trying with the index pointers.
    {
        DimensionLimit<int> csub;
        csub.type = tatami::DimensionLimitType::INDEX;
        csub.index_length = indices.size();
        csub.index_start = indices.data();

        auto rwork = ref->dense_row();
        test_simple_row_access_base(ptr, ref, forward, jump, 
            [&](size_t r) -> auto { 
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
            csub
        );

        // Checking that properties are correctly passed down.
        auto pwork = ptr->dense_row(DimensionLimit<int>(), csub, ExtractOptions());
        EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(pwork->extracted_index(), indices.data());
        EXPECT_EQ(pwork->extracted_length, indices.size());

        auto swork = ptr->sparse_row(DimensionLimit<int>(), csub, ExtractOptions());
        EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(swork->extracted_index(), indices.data());
        EXPECT_EQ(swork->extracted_length, indices.size());
    }

    // First trying with the index pointers.
    {
        DimensionLimit<int> csub;
        csub.type = tatami::DimensionLimitType::INDEX;
        csub.indices = indices;

        auto rwork = ref->dense_row();
        test_simple_row_access_base(ptr, ref, forward, jump, 
            [&](size_t r) -> auto { 
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
            csub
        );

        // Checking that properties are correctly passed down.
        auto pwork = ptr->dense_row(DimensionLimit<int>(), csub, ExtractOptions());
        EXPECT_EQ(pwork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(std::vector<int>(pwork->extracted_index(), pwork->extracted_index() + pwork->extracted_length), indices);
        EXPECT_EQ(pwork->extracted_length, indices.size());

        auto swork = ptr->sparse_row(DimensionLimit<int>(), csub, ExtractOptions());
        EXPECT_EQ(swork->extracted_limit, tatami::DimensionLimitType::INDEX);
        EXPECT_EQ(std::vector<int>(swork->extracted_index(), swork->extracted_index() + swork->extracted_length), indices);
        EXPECT_EQ(swork->extracted_length, indices.size());
    }
}

#endif
