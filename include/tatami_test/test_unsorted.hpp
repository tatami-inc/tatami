#ifndef TATAMI_TEST_UNSORTED_HPP
#define TATAMI_TEST_UNSORTED_HPP

#include "fetch.hpp"
#include "test_access.hpp"

#include "../tatami/utils/new_extractor.hpp"
#include "../tatami/utils/ConsecutiveOracle.hpp"
#include "../tatami/utils/FixedOracle.hpp"

#include <vector>
#include <limits>
#include <random>
#include <cmath>
#include <memory>

namespace tatami_test {

namespace internal {

template<bool use_oracle_, typename Value_, typename Index_, typename ...Args_>
void test_unsorted_access_base(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr, Index_ extent, Args_... args) {
    auto NR = ptr->nrow();
    auto NC = ptr->ncol();

    auto sequence = simulate_sequence(params, NR, NC);
    auto oracle = create_oracle<use_oracle_>(params, sequence);

    auto swork = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args...);

    tatami::Options opt;
    opt.sparse_ordered_index = false;
    auto swork_uns = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    opt.sparse_extract_index = false;
    auto swork_uns_v = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    opt.sparse_extract_value = false;
    auto swork_uns_n = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    opt.sparse_extract_index = true;
    auto swork_uns_i = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    // Looping over rows/columns and checking extraction for various unsorted combinations.
    for (auto i : sequence) {
        auto observed = [&]() {
            if constexpr(use_oracle_) {
                return fetch(swork.get(), extent);
            } else {
                return fetch(swork.get(), i, extent);
            }
        }();

        {
            auto observed_uns = [&]() {
                if constexpr(use_oracle_) {
                    return fetch(swork_uns.get(), extent);
                } else {
                    return fetch(swork_uns.get(), i, extent);
                }
            }();

            // Poor man's zip + unzip.
            std::vector<std::pair<Index_, Value_> > collected;
            collected.reserve(observed.value.size());
            for (Index_ i = 0, end = observed_uns.value.size(); i < end; ++i) {
                collected.emplace_back(observed_uns.index[i], observed_uns.value[i]);
            }
            std::sort(collected.begin(), collected.end());

            std::vector<Value_> sorted_v;
            std::vector<Index_> sorted_i;
            sorted_v.reserve(collected.size());
            sorted_i.reserve(collected.size());
            for (const auto& p : collected) {
                sorted_i.push_back(p.first);
                sorted_v.push_back(p.second);
            }

            ASSERT_EQ(sorted_i, observed.index);
            ASSERT_EQ(sorted_v, observed.value);
        }

        {
            std::vector<int> indices(extent);
            auto observed_i = [&]() {
                if constexpr(use_oracle_) {
                    return swork_uns_i->fetch(NULL, indices.data());
                } else {
                    return swork_uns_i->fetch(i, NULL, indices.data());
                }
            }();
            ASSERT_TRUE(observed_i.value == NULL);

            tatami::copy_n(observed_i.index, observed_i.number, indices.data());
            indices.resize(observed_i.number);
            std::sort(indices.begin(), indices.end());
            ASSERT_EQ(observed.index, indices);
        }

        {
            std::vector<double> vbuffer(extent);
            auto observed_v = [&]() {
                if constexpr(use_oracle_) {
                    return swork_uns_v->fetch(vbuffer.data(), NULL);
                } else {
                    return swork_uns_v->fetch(i, vbuffer.data(), NULL);
                }
            }();
            ASSERT_TRUE(observed_v.index == NULL);

            tatami::copy_n(observed_v.value, observed_v.number, vbuffer.data());
            vbuffer.resize(observed_v.number);
            sanitize_nan(vbuffer, params.has_nan);
            std::sort(vbuffer.begin(), vbuffer.end());

            auto vexpected = observed.value;
            std::sort(vexpected.begin(), vexpected.end());
            ASSERT_EQ(vexpected, vbuffer);
        }

        {
            auto observed_n = [&]() {
                if constexpr(use_oracle_) {
                    return swork_uns_n->fetch(NULL, NULL);
                } else {
                    return swork_uns_n->fetch(i, NULL, NULL);
                }
            }();
            ASSERT_TRUE(observed_n.value == NULL);
            ASSERT_TRUE(observed_n.index == NULL);
            ASSERT_EQ(observed.value.size(), observed_n.number);
        }
    }
}

template<bool use_oracle_, typename Value_, typename Index_>
void test_unsorted_full_access(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr) {
    internal::test_unsorted_access_base<use_oracle_>(params, ptr, (params.use_row ? ptr->ncol() : ptr->nrow()));
}

template<bool use_oracle_, typename Value_, typename Index_>
void test_unsorted_block_access(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr, int start, int end) {
    auto len = end - start;
    internal::test_unsorted_access_base<use_oracle_>(params, ptr, len, start, len);
}

template<bool use_oracle_, typename Value_, typename Index_>
void test_unsorted_indexed_access(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr, int start, int step) {
    std::vector<Index_> indices;
    int counter = start;
    int nsecondary = (params.use_row ? ptr->ncol() : ptr->nrow());
    while (counter < nsecondary) {
        indices.push_back(counter);
        counter += step;
    }
    internal::test_unsorted_access_base<use_oracle_>(params, ptr, static_cast<Index_>(indices.size()), indices);
}

}

/********************************************************
 ********************************************************/

// Finally, some user-visible functions that convert compile-time parameters to run-time
// parameters, to make it easier to define parametrized tests across oracle/row status.
template<typename Value_, typename Index_>
void test_unsorted_full_access(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr) {
    if (params.use_oracle) {
        internal::test_unsorted_full_access<true>(params, ptr);
    } else {
        internal::test_unsorted_full_access<false>(params, ptr);
    }
}

template<typename Value_, typename Index_>
void test_unsorted_block_access(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr, int start, int end) {
    if (params.use_oracle) {
        internal::test_unsorted_block_access<true>(params, ptr, start, end);
    } else {
        internal::test_unsorted_block_access<false>(params, ptr, start, end);
    }
}

template<typename Value_, typename Index_>
void test_unsorted_indexed_access(const TestAccessParameters& params, const tatami::Matrix<Value_, Index_>* ptr, int start, int step) {
    if (params.use_oracle) {
        internal::test_unsorted_indexed_access<true>(params, ptr, start, step);
    } else {
        internal::test_unsorted_indexed_access<false>(params, ptr, start, step);
    }
}

}

#endif
