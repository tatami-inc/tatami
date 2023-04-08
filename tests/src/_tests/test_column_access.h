#ifndef TEST_COLUMN_ACCESS_H
#define TEST_COLUMN_ACCESS_H
#include "utils.h"

template<class Matrix, class Matrix2, class Function1, class Function2, typename ...Args>
void test_simple_column_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, Function1 expector, Function2 sparse_expand, Args... args) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto pwork = ptr->dense_column_workspace(args...);
    auto swork = ptr->sparse_column_workspace(args...);

    tatami::WorkspaceOptions opt;
    opt.sorted = false;
    auto swork_uns = ptr->sparse_column_workspace(args..., opt);
    opt.sorted = true;

    opt.mode = tatami::SparseExtractMode::INDEX;
    auto swork_i = ptr->sparse_column_workspace(args..., opt);

    opt.mode = tatami::SparseExtractMode::VALUE;
    auto swork_v = ptr->sparse_column_workspace(args..., opt);

    opt.mode = tatami::SparseExtractMode::NONE;
    auto swork_n = ptr->sparse_column_workspace(args..., opt);
    opt.mode = tatami::SparseExtractMode::BOTH;

    auto pwork_bi = ptr->dense_column_workspace(args...);
    auto rwork_bi = ref->dense_column_workspace(args...);

    opt.cache = true;
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

            auto observed_i = ptr->column(c, swork_i.get());
            EXPECT_EQ(observed_i.index, observed.index);

            auto observed_v = ptr->column(c, swork_v.get());
            EXPECT_EQ(observed_v.value, observed_v.value);

            auto observed_n = ptr->column(c, NULL, NULL, swork_n.get());
            EXPECT_EQ(observed.value.size(), observed_n.number);
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
}

#endif
