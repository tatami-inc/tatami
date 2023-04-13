#ifndef TEST_COLUMN_ACCESS_H
#define TEST_COLUMN_ACCESS_H
#include "utils.h"
#include <type_traits>

template<class Matrix, class Matrix2, class Function1, class Function2, typename ...Args>
void test_simple_column_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, Function1 expector, Function2 sparse_expand, Args... args) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto pwork = ptr->dense_column_workspace(args...);
    auto swork = ptr->sparse_column_workspace(args...);

    tatami::WorkspaceOptions opt;
    opt.sparse_ordered_index = false;
    auto swork_uns = ptr->sparse_column_workspace(args..., opt);
    opt.sparse_ordered_index = true;

    opt.sparse_extract_mode = tatami::SparseExtractMode::INDEX;
    auto swork_i = ptr->sparse_column_workspace(args..., opt);

    opt.sparse_extract_mode = tatami::SparseExtractMode::VALUE;
    auto swork_v = ptr->sparse_column_workspace(args..., opt);

    opt.sparse_extract_mode = tatami::SparseExtractMode::NONE;
    auto swork_n = ptr->sparse_column_workspace(args..., opt);
    opt.sparse_extract_mode = tatami::SparseExtractMode::BOTH;

    auto pwork_bi = ptr->dense_column_workspace(args...);
    auto rwork_bi = ref->dense_column_workspace(args...);

    opt.cache_for_reuse = true;
    auto cwork = ptr->dense_column_workspace(args..., opt); 

    for (size_t i = 0; i < NC; i += jump) {
        size_t c = (forward ? i : NC - i - 1);

        auto expected = expector(c);

        {
            auto observed = ptr->column(c, pwork.get());
            EXPECT_EQ(expected, observed);
        }

        {
            auto observed = ptr->column(c, swork.get());
            EXPECT_EQ(expected, sparse_expand(observed));
            EXPECT_TRUE(is_increasing(observed.index));

            std::vector<int> indices(expected.size()); // using the dense expected size as a proxy for the extraction length in block/indexed cases.
            auto observed_i = ptr->column(c, NULL, indices.data(), swork_i.get());
            EXPECT_TRUE(observed_i.value == NULL);
            EXPECT_EQ(observed.index, std::vector<int>(observed_i.index, observed_i.index + observed_i.number));

            std::vector<double> values(expected.size());
            auto observed_v = ptr->column(c, values.data(), NULL, swork_v.get());
            EXPECT_TRUE(observed_v.index == NULL);
            EXPECT_EQ(observed.value, std::vector<double>(observed_v.value, observed_v.value + observed_v.number));

            auto observed_n = ptr->column(c, NULL, NULL, swork_n.get());
            EXPECT_TRUE(observed_n.value == NULL);
            EXPECT_TRUE(observed_n.index == NULL);
            EXPECT_EQ(observed.value.size(), observed_n.number);

            auto observed_n2 = ptr->column(c, swork_n.get()); // just another request for some coverage of the Matrix::copy_over function.
            EXPECT_EQ(observed.value.size(), observed_n2.value.size());
        }

        {
            auto observed = ptr->column(c, swork_uns.get());
            EXPECT_EQ(expected, sparse_expand(observed));
        } 

        // Check for proper workspace behavior when access is bidirectional,
        // i.e., not purely increasing or decreasing.
        if (jump > 1 && i) {
            auto observed = ptr->column(c, pwork_bi.get());
            EXPECT_EQ(expected, observed);

            auto subc = (forward ? c-1 : c+1);
            auto observedm1 = ptr->column(subc, pwork_bi.get());
            auto expectedm1 = ref->column(subc, rwork_bi.get());
            EXPECT_EQ(expectedm1, observedm1);
        }

        {
            auto observed = ptr->column(c, cwork.get());
            EXPECT_EQ(expected, observed);
        }
    }

    // Check that the caching gives the same results as uncached pass.
    for (size_t i = 0; i < NC; i += jump) {
        size_t c = (forward ? i : NC - i - 1);
        auto expected = ptr->column(c, pwork.get());
        auto observed = ptr->column(c, cwork.get());
        EXPECT_EQ(expected, observed);
    }
}

template<class Matrix, class Matrix2>
void test_simple_column_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, size_t jump = 1) {
    auto rwork = ref->dense_column_workspace();
    size_t NR = ref->nrow();

    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](size_t c) -> auto { 
            auto expected = ref->column(c, rwork.get());
            EXPECT_EQ(expected.size(), NR);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, NR);
        }
    );
}

template<class Matrix, class Matrix2>
void test_sliced_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t end) {
    auto rwork = ref->dense_column_workspace();
    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](size_t c) -> auto { 
            auto raw_expected = ref->column(c, rwork.get());
            return std::vector<typename Matrix::data_type>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& range) -> auto {
            return expand(range, start, end);
        },
        start, 
        end - start
    );

    // Checking that properties are correctly passed down.
    auto pwork = ptr->dense_column_workspace(start, end - start);
    EXPECT_EQ(pwork->start, start);
    EXPECT_EQ(pwork->length, end - start);

    auto swork = ptr->sparse_column_workspace(start, end - start);
    EXPECT_EQ(pwork->start, start);
    EXPECT_EQ(pwork->length, end - start);
}

template<class Matrix, class Matrix2>
void test_indexed_column_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t step) {
    size_t NR = ptr->nrow();
    std::vector<typename Matrix::index_type> indices;
    {
        int counter = start;
        while (counter < NR) {
            indices.push_back(counter);
            counter += step;
        }
    }

    auto rwork = ref->dense_column_workspace();
    test_simple_column_access_base(ptr, ref, forward, jump, 
        [&](size_t c) -> auto { 
            auto raw_expected = ref->column(c, rwork.get());
            std::vector<typename Matrix::data_type> expected;
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
    auto pwork = ptr->dense_column_workspace(indices);
    EXPECT_EQ(pwork->indices(), indices);
    EXPECT_EQ(pwork->length, indices.size());

    auto swork = ptr->sparse_column_workspace(indices);
    EXPECT_EQ(pwork->indices(), indices);
    EXPECT_EQ(pwork->length, indices.size());
}

#endif
