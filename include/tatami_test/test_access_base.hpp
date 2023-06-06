#ifndef TATAMI_TEST_ACCESS_BASE_HPP
#define TATAMI_TEST_ACCESS_BASE_HPP

#include "utils.hpp"
#include <type_traits>
#include "../tatami/base/utils.hpp"

namespace tatami_test {

template<bool use_row_, bool has_nan_, bool less_sparse_, class Matrix_, class Matrix2_, class Function1_, class Function2_, typename ...Args_>
void test_access_base(const Matrix_* ptr, const Matrix2_* ref, bool forward, int jump, Function1_ expector, Function2_ sparse_expand, Args_... args) {
    int NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    int NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto pwork = tatami::new_extractor<use_row_, false>(ptr, args...);
    auto swork = tatami::new_extractor<use_row_, true>(ptr, args...);

    tatami::Options opt;
    opt.sparse_ordered_index = false;
    auto swork_uns = tatami::new_extractor<use_row_, true>(ptr, args..., opt);
    opt.sparse_ordered_index = true;

    opt.sparse_extract_index = false;
    auto swork_v = tatami::new_extractor<use_row_, true>(ptr, args..., opt);
    opt.sparse_extract_value = false;
    auto swork_n = tatami::new_extractor<use_row_, true>(ptr, args..., opt);
    opt.sparse_extract_index = true;
    auto swork_i = tatami::new_extractor<use_row_, true>(ptr, args..., opt);
    opt.sparse_extract_value = true;

    opt.cache_for_reuse = true;
    auto cwork = tatami::new_extractor<use_row_, false>(ptr, args...);

    size_t sparse_extracted = 0, dense_extracted = 0;
    int limit = (use_row_ ? NR : NC);

    for (int i = 0; i < limit; i += jump) {
        int c = (forward ? i : limit - i - 1);

        auto expected = expector(c);
        sanitize_nan<has_nan_>(expected);
        dense_extracted += expected.size();

        {
            auto observed = pwork->fetch(c);
            sanitize_nan<has_nan_>(observed);
            EXPECT_EQ(expected, observed);
        }

        // Various flavors of sparse retrieval.
        {
            auto observed = swork->fetch(c);
            sanitize_nan<has_nan_>(observed.value);
            EXPECT_EQ(expected, sparse_expand(observed));
            EXPECT_TRUE(is_increasing(observed.index));

            sparse_extracted += observed.number;

            std::vector<int> indices(expected.size()); // using the dense expected size as a proxy for the extraction length in block/indexed cases.
            auto observed_i = swork_i->fetch(c, NULL, indices.data());
            EXPECT_TRUE(observed_i.value == NULL);
            EXPECT_EQ(observed.index, std::vector<int>(observed_i.index, observed_i.index + observed_i.number));

            std::vector<double> vbuffer(expected.size());
            auto observed_v = swork_v->fetch(c, vbuffer.data(), NULL);
            EXPECT_TRUE(observed_v.index == NULL);
            std::vector<double> values_only(observed_v.value, observed_v.value + observed_v.number);
            sanitize_nan<has_nan_>(values_only);
            EXPECT_EQ(observed.value, values_only);

            auto observed_n = swork_n->fetch(c, NULL, NULL);
            EXPECT_TRUE(observed_n.value == NULL);
            EXPECT_TRUE(observed_n.index == NULL);
            EXPECT_EQ(observed.value.size(), observed_n.number);

            auto observed_uns = swork_uns->fetch(c);
            sanitize_nan<has_nan_>(observed_uns.value);
            EXPECT_EQ(expected, sparse_expand(observed_uns));
        } 

        // Checks for caching.
        {
            auto observed = cwork->fetch(c);
            sanitize_nan<has_nan_>(observed);
            EXPECT_EQ(expected, observed);
        }
    }

    // Check that the caching gives the same results as uncached pass.
    for (int i = 0; i < limit; i += jump) {
        int c = (forward ? i : limit - i - 1);
        auto expected = pwork->fetch(c);
        auto observed = cwork->fetch(c);

        sanitize_nan<has_nan_>(expected);
        sanitize_nan<has_nan_>(observed);
        EXPECT_EQ(expected, observed);
    }

    // Check that sparse extraction isn't just the same as the dense extractor. 
    if constexpr(less_sparse_) {
        if (ptr->sparse()) {
            EXPECT_TRUE(dense_extracted > sparse_extracted);
        }
    }
}

}

#endif
