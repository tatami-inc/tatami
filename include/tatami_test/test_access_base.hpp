#ifndef TATAMI_TEST_ACCESS_BASE_HPP
#define TATAMI_TEST_ACCESS_BASE_HPP

#include "utils.hpp"
#include <type_traits>
#include "../tatami/base/utils.hpp"

namespace tatami_test {

enum TestAccessOrder { FORWARD, REVERSE, RANDOM };

struct TestAccessParameters {
    // Whether to expect NaNs (and sanitize them in the comparison).
    bool has_nan;

    // Whether the input matrix is actually sparse and should have fewer values from the sparse extractors. 
    bool less_sparse; 

    // Ordering of rows/columns to test.
    Order order;

    // Minimum jump between rows/columns.
    int jump;
};

template<bool use_row_, class TestMatrix_, class RefMatrix_, class DenseExtract_, class SparseExpand_, typename ...Args_>
void test_access_base_no_oracle(const TestAccessParameters& params, const TestMatrix_* ptr, const RefMatrix_* ref, DenseExtract_ expector, SparseExpand_ sparse_expand, Args_... args) {
    int NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    int NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());
    int limit = (use_row_ ? NR : NC);

    // Setting up the sequence over which to iterate.
    typename typename TestMatrix_::index_type Index_;
    std::shared_ptr<tatami::Oracle<Index_> > oracle;
    std::vector<Index_> sequence;
    if constexpr(use_oracle_) {
        if (jump == 1 && order == FORWARD) {
            oracle.reset(new tatami::ConsecutiveOracle<Index_>(0, limit);
        } else {
            if (order == REVERSE) {
                for (int i = limit; i > 0; i -= jump) {
                    sequence.push_back(i - 1);
                }
            } else {
                for (int i = 0; i < limit; i += jump) {
                    sequence.push_back(i);
                }
                if (order == RANDOM) {
                    uint64_t seed = static_cast<uint64_t>(NR) * NC + limit * static_cast<uint64_t>(order);
                    std::mt19937_64 rng(seed);
                    std::random_shuffle(sequence.begin(), sequence.end(), rng);
                }
            }
            oracle.reset(new tatami::FixedOracle<Index_>(sequence.data(), sequence.size()));
        }
    }

    // Constructing the various extractors.
    auto pwork = [&]() {
        if constexpr(use_oracle_) {
            return tatami::new_extractor<use_row_, false>(ptr, oracle, args...);
        } else {
            return tatami::new_extractor<use_row_, false>(ptr, args...);
        }
    }();
    auto swork = [&]() {
        if constexpr(use_oracle_) {
            return tatami::new_extractor<use_row_, true>(ptr, oracle, args...);
        } else {
            return tatami::new_extractor<use_row_, true>(ptr, args...);
        }
    }();

    tatami::Options opt;
    opt.sparse_ordered_index = false;
    auto swork_uns = [&]() {
        if constexpr(use_oracle_) {
            return tatami::new_extractor<use_row_, true>(ptr, oracle, args..., opt);
        } else {
            return tatami::new_extractor<use_row_, true>(ptr, args..., opt);
        }
    }();
    opt.sparse_ordered_index = true;

    opt.sparse_extract_index = false;
    auto swork_v = [&]() {
        if constexpr(use_oracle_) {
            return tatami::new_extractor<use_row_, true>(ptr, oracle, args..., opt);
        } else {
            return tatami::new_extractor<use_row_, true>(ptr, args..., opt);
        }
    }();

    opt.sparse_extract_value = false;
    auto swork_n = [&]() {
        if constexpr(use_oracle_) {
            return tatami::new_extractor<use_row_, true>(ptr, oracle, args..., opt);
        } else {
            return tatami::new_extractor<use_row_, true>(ptr, args..., opt);
        }
    }();

    opt.sparse_extract_index = true;
    auto swork_i = [&]() {
        if constexpr(use_oracle_) {
            return tatami::new_extractor<use_row_, true>(ptr, oracle, args..., opt);
        } else {
            return tatami::new_extractor<use_row_, true>(ptr, args..., opt);
        }
    }();
    opt.sparse_extract_value = true;

    // Looping over rows/columns and checking extraction against the reference.
    size_t sparse_extracted = 0, dense_extracted = 0;

    auto get_index = [&](int i0) -> Index_ {
        Index_ i;
        switch (order) {
            case FORWARD: i = i0; break;
            case REVERSE: i = limit - i0 - 1; break;
            case RANDOM: i = sequence[i0/jump]; break;
        }
        return i;
    };

    for (int i0 = 0; i0 < limit; i0 += jump) {
        auto i = get_index(i0);
        auto expected = expector(i);
        sanitize_nan(expected, has_nan);
        dense_extracted += expected.size();

        // Checking dense retrieval first.
        {
            auto observed = [&]() {
                if constexpr(use_oracle_) {
                    Index_ j = limit;
                    auto output = pwork->fetch(j);
                    EXPECT_EQ(i, j);
                    return output;
                } else {
                    return pwork->fetch(i);
                }
            }();
            sanitize_nan(observed, has_nan);
            ASSERT_EQ(expected, observed);
        }

        // Various flavors of sparse retrieval.
        {
            auto observed = [&]() {
                if constexpr(use_oracle_) {
                    Index_ j = limit;
                    auto output = swork->fetch(j);
                    EXPECT_EQ(i, j);
                    return output;
                } else {
                    return swork->fetch(i);
                }
            }();
            sanitize_nan(observed.value, has_nan);
            ASSERT_EQ(expected, sparse_expand(observed));
            ASSERT_TRUE(is_increasing(observed.index));

            sparse_extracted += observed.number;

            std::vector<int> indices(expected.size()); // using the dense expected size as a proxy for the extraction length in block/indexed cases.
            auto observed_i = [&]() {
                if constexpr(use_oracle_) {
                    Index_ j = limit;
                    auto output = swork_i->fetch(j, NULL, indices.data());
                    EXPECT_EQ(i, j);
                    return output;
                } else {
                    return swork_i->fetch(i, NULL, indices.data());
                }
            }();
            ASSERT_TRUE(observed_i.value == NULL);
            ASSERT_EQ(observed.index, std::vector<int>(observed_i.index, observed_i.index + observed_i.number));

            std::vector<double> vbuffer(expected.size());
            auto observed_v = [&]() {
                if constexpr(use_oracle_) {
                    Index_ j = limit;
                    auto output = swork_v->fetch(j, vbuffer.data(), NULL);
                    EXPECT_EQ(i, j);
                    return output;
                } else {
                    return swork_v->fetch(i, vbuffer.data(), NULL);
                }
            }();
            ASSERT_TRUE(observed_v.index == NULL);
            std::vector<double> values_only(observed_v.value, observed_v.value + observed_v.number);
            sanitize_nan(values_only, has_nan);
            ASSERT_EQ(observed.value, values_only);

            auto observed_n = [&]() {
                if constexpr(use_oracle_) {
                    Index_ j = limit;
                    auto output = swork_n->fetch(j, NULL, NULL);
                    EXPECT_EQ(i, j);
                    return output;
                } else {
                    return swork_n->fetch(i, NULL, NULL);
                }
            }();
            ASSERT_TRUE(observed_n.value == NULL);
            ASSERT_TRUE(observed_n.index == NULL);
            ASSERT_EQ(observed.value.size(), observed_n.number);

            auto observed_uns = swork_uns->fetch(c);
            sanitize_nan(observed_uns.value, has_nan);
            ASSERT_EQ(expected, sparse_expand(observed_uns));
        } 
    }

    // Check that sparse extraction isn't just the same as the dense extractor. 
    if constexpr(less_sparse_) {
        if (ptr->sparse()) {
            ASSERT_TRUE(dense_extracted > sparse_extracted);
        }
    }
}

}

#endif
