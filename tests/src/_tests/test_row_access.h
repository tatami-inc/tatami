#ifndef TEST_ROW_ACCESS_H
#define TEST_ROW_ACCESS_H
#include "utils.h"

template<class Matrix, class Matrix2, class Function1, class Function2, typename ...Args>
void test_simple_row_access_base(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, Function1 expector, Function2 sparse_expand, Args... args) {
    size_t NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    size_t NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto pwork = ptr->dense_row_workspace(args...);
    auto swork = ptr->sparse_row_workspace(args...);

    tatami::WorkspaceOptions opt;
    opt.sorted = false;
    auto swork_uns = ptr->sparse_row_workspace(args..., opt);
    opt.sorted = true;

    opt.mode = tatami::SparseExtractMode::INDEX;
    auto swork_i = ptr->sparse_row_workspace(args..., opt);

    opt.mode = tatami::SparseExtractMode::VALUE;
    auto swork_v = ptr->sparse_row_workspace(args..., opt);

    opt.mode = tatami::SparseExtractMode::NONE;
    auto swork_n = ptr->sparse_row_workspace(args..., opt);
    opt.mode = tatami::SparseExtractMode::BOTH;

    auto pwork_bi = ptr->dense_row_workspace(args...);
    auto rwork_bi = ref->dense_row_workspace(args...);

    opt.cache = true;
    auto cwork = ptr->dense_row_workspace(args..., opt); 

    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);

        auto expected = expector(r);
        {
            auto observed = ptr->row(r, pwork.get());
            EXPECT_EQ(expected, observed);
        }

        {
            auto observed = ptr->row(r, swork.get());
            EXPECT_EQ(expected, sparse_expand(observed));
            EXPECT_TRUE(is_increasing(observed.index));

            std::vector<int> indices(expected.size()); // using the dense expected size as a proxy for the extraction length in block/indexed cases.
            auto observed_i = ptr->row(r, NULL, indices.data(), swork_i.get());
            EXPECT_TRUE(observed_i.value == NULL);
            EXPECT_EQ(observed.index, std::vector<int>(observed_i.index, observed_i.index + observed_i.number));

            std::vector<double> values(expected.size());
            auto observed_v = ptr->row(r, values.data(), NULL, swork_v.get());
            EXPECT_TRUE(observed_v.index == NULL);
            EXPECT_EQ(observed.value, std::vector<double>(observed_v.value, observed_v.value + observed_v.number));

            auto observed_n = ptr->row(r, NULL, NULL, swork_n.get());
            EXPECT_TRUE(observed_n.value == NULL);
            EXPECT_TRUE(observed_n.index == NULL);
            EXPECT_EQ(observed.value.size(), observed_n.number);

            auto observed_n2 = ptr->row(r, swork_n.get()); // just another request for some coverage of the Matrix::copy_over function.
            EXPECT_EQ(observed.value.size(), observed_n2.value.size());
        }

        {
            auto observed = ptr->row(r, swork_uns.get());
            EXPECT_EQ(expected, sparse_expand(observed));
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

        {
            auto observed = ptr->row(r, cwork.get());
            EXPECT_EQ(expected, observed);
        }
    }

    // Check that the caching gives the same results as uncached pass.
    for (size_t i = 0; i < NR; i += jump) {
        size_t r = (forward ? i : NR - i - 1);
        auto expected = ptr->row(r, pwork.get());
        auto observed = ptr->row(r, cwork.get());
        EXPECT_EQ(expected, observed);
    }
}

template<class Matrix, class Matrix2>
void test_simple_row_access(const Matrix* ptr, const Matrix2* ref, bool forward = true, size_t jump = 1) {
    auto rwork = ref->dense_row_workspace();
    size_t NC = ref->ncol();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](size_t r) -> auto { 
            auto expected = ref->row(r, rwork.get());
            EXPECT_EQ(expected.size(), NC);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, NC);
        }
    );
}

template<class Matrix, class Matrix2>
void test_sliced_row_access(const Matrix* ptr, const Matrix2* ref, bool forward, size_t jump, size_t start, size_t end) {
    auto rwork = ref->dense_row_workspace();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](size_t r) -> auto { 
            auto raw_expected = ref->row(r, rwork.get());
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

    auto rwork = ref->dense_row_workspace();
    test_simple_row_access_base(ptr, ref, forward, jump, 
        [&](size_t r) -> auto { 
            auto raw_expected = ref->row(r, rwork.get());
            std::vector<typename Matrix::data_type> expected;
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
}

#endif
