#ifndef TATAMI_TEST_ACCESS_HPP
#define TATAMI_TEST_ACCESS_HPP

#include <gtest/gtest.h>

#include "fetch.hpp"
#include "../tatami/utils/new_extractor.hpp"
#include "../tatami/utils/ConsecutiveOracle.hpp"
#include "../tatami/utils/FixedOracle.hpp"

#include <vector>
#include <limits>
#include <random>
#include <cmath>
#include <memory>

namespace tatami_test {

enum TestAccessOrder { FORWARD, REVERSE, RANDOM };

struct TestAccessParameters {
    // Whether to use an oracle.
    bool use_oracle = false;

    // Whether to test access across rows.
    bool use_row = true;

    // Ordering of rows/columns to test.
    TestAccessOrder order = FORWARD;

    // Minimum jump between rows/columns.
    int jump = 1;

    // Whether to expect NaNs (and sanitize them in the comparison).
    bool has_nan = false;

    // Whether to check that "sparse" matrices are actually sparse.
    bool check_sparse = true;

};

typedef std::tuple<bool, bool, TestAccessOrder, int> StandardTestAccessParameters;

inline TestAccessParameters convert_access_parameters(const StandardTestAccessParameters& tup) {
    TestAccessParameters output;
    output.use_row = std::get<0>(tup);
    output.use_oracle = std::get<1>(tup);
    output.order = std::get<2>(tup);
    output.jump = std::get<3>(tup);
    return output;
}

inline auto standard_test_access_parameter_combinations() {
    return ::testing::Combine(
        ::testing::Values(true, false), /* whether to access the rows. */
        ::testing::Values(true, false), /* whether to use an oracle. */
        ::testing::Values(tatami_test::FORWARD, tatami_test::REVERSE, tatami_test::RANDOM), /* access order. */
        ::testing::Values(1, 3) /* jump between rows/columns. */
    );
}

/********************************************************
 ********************************************************/

namespace internal {

template<typename Value_>
void sanitize_nan(std::vector<Value_>& values, bool has_nan) {
    if constexpr(std::numeric_limits<Value_>::has_quiet_NaN) {
        if (has_nan) {
            for (auto& x : values) {
                if (std::isnan(x)) {
                    // Using an appropriately silly placeholder that
                    // should never occur in our test data.
                    x = 1234567890;
                }
            }
        }
    }
}

template<typename Index_>
std::vector<Index_> simulate_sequence(const TestAccessParameters& params, Index_ NR, Index_ NC) {
    std::vector<Index_> sequence;
    auto limit = (params.use_row ? NR : NC);
    if (params.order == REVERSE) {
        for (int i = limit; i > 0; i -= params.jump) {
            sequence.push_back(i - 1);
        }
    } else {
        for (int i = 0; i < limit; i += params.jump) {
            sequence.push_back(i);
        }
        if (params.order == RANDOM) {
            uint64_t seed = static_cast<uint64_t>(NR) * NC + limit * static_cast<uint64_t>(params.order);
            std::mt19937_64 rng(seed);
            std::shuffle(sequence.begin(), sequence.end(), rng);
        }
    }
    return sequence;
}

template<bool use_oracle_, typename Index_>
auto create_oracle(const TestAccessParameters& params, const std::vector<Index_>& sequence) {
    if constexpr(use_oracle_) {
        std::shared_ptr<tatami::Oracle<Index_> > oracle;
        if (params.jump == 1 && params.order == FORWARD) {
            oracle.reset(new tatami::ConsecutiveOracle<Index_>(0, sequence.size()));
        } else {
            oracle.reset(new tatami::FixedViewOracle<Index_>(sequence.data(), sequence.size()));
        }
        return oracle;
    } else {
        return false;
    }
}

template<bool use_oracle_, typename Value_, typename Index_, class DenseExtract_, class SparseExpand_, typename ...Args_>
void test_access_base(
    const TestAccessParameters& params, 
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref, 
    DenseExtract_ expector, 
    SparseExpand_ sparse_expand, 
    Args_... args) 
{
    auto NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    auto NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());

    auto sequence = simulate_sequence(params, NR, NC);
    auto oracle = create_oracle<use_oracle_>(params, sequence);

    auto pwork = tatami::new_extractor<false, use_oracle_>(ptr, params.use_row, oracle, args...);
    auto swork = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args...);

    tatami::Options opt;
    opt.sparse_extract_index = false;
    auto swork_v = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    opt.sparse_extract_value = false;
    auto swork_n = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    opt.sparse_extract_index = true;
    auto swork_i = tatami::new_extractor<true, use_oracle_>(ptr, params.use_row, oracle, args..., opt);

    size_t sparse_counter = 0;

    // Looping over rows/columns and checking extraction against the reference.
    for (auto i : sequence) {
        auto expected = expector(i);
        sanitize_nan(expected, params.has_nan);
        auto extent = expected.size(); // using the dense expected size to determine the expected extraction length.

        // Checking dense retrieval first.
        {
            auto observed = [&]() {
                if constexpr(use_oracle_) {
                    return fetch(pwork.get(), extent);
                } else {
                    return fetch(pwork.get(), i, extent);
                }
            }();
            sanitize_nan(observed, params.has_nan);
            ASSERT_EQ(expected, observed);
        }

        // Various flavors of sparse retrieval.
        {
            auto observed = [&]() {
                if constexpr(use_oracle_) {
                    return fetch(swork.get(), extent);
                } else {
                    return fetch(swork.get(), i, extent);
                }
            }();

            sanitize_nan(observed.value, params.has_nan);
            ASSERT_EQ(expected, sparse_expand(observed));

            sparse_counter += observed.value.size();
            {
                bool is_increasing = true;
                for (size_t i = 1; i < observed.index.size(); ++i) {
                    if (observed.index[i] <= observed.index[i-1]) {
                        is_increasing = false;
                        break;
                    }
                }
                ASSERT_TRUE(is_increasing);
            }

            std::vector<Index_> indices(extent);
            auto observed_i = [&]() {
                if constexpr(use_oracle_) {
                    return swork_i->fetch(NULL, indices.data());
                } else {
                    return swork_i->fetch(i, NULL, indices.data());
                }
            }();
            ASSERT_TRUE(observed_i.value == NULL);
            tatami::copy_n(observed_i.index, observed_i.number, indices.data());
            indices.resize(observed_i.number);
            ASSERT_EQ(observed.index, indices);

            std::vector<Value_> vbuffer(extent);
            auto observed_v = [&]() {
                if constexpr(use_oracle_) {
                    return swork_v->fetch(vbuffer.data(), NULL);
                } else {
                    return swork_v->fetch(i, vbuffer.data(), NULL);
                }
            }();
            ASSERT_TRUE(observed_v.index == NULL);
            tatami::copy_n(observed_v.value, observed_v.number, vbuffer.data());
            vbuffer.resize(observed_v.number);
            sanitize_nan(vbuffer, params.has_nan);
            ASSERT_EQ(observed.value, vbuffer);

            auto observed_n = [&]() {
                if constexpr(use_oracle_) {
                    return swork_n->fetch(NULL, NULL);
                } else {
                    return swork_n->fetch(i, NULL, NULL);
                }
            }();
            ASSERT_TRUE(observed_n.value == NULL);
            ASSERT_TRUE(observed_n.index == NULL);
            ASSERT_EQ(observed.value.size(), observed_n.number);
        } 
    }

    if (params.check_sparse && ptr->is_sparse()) {
        EXPECT_TRUE(sparse_counter < static_cast<size_t>(NR) * static_cast<size_t>(NC));
    }
}

template<bool use_oracle_, typename Value_, typename Index_>
void test_full_access(
    const TestAccessParameters& params, 
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref)
{
    int nsecondary = (params.use_row ? ref->ncol() : ref->nrow());
    auto refwork = (params.use_row ? ref->dense_row() : ref->dense_column());
    test_access_base<use_oracle_>(
        params,
        ptr, 
        ref, 
        [&](int i) -> auto { 
            return fetch(refwork.get(), i, nsecondary);
        },
        [&](const auto& svec) -> auto {
            std::vector<Value_> output(nsecondary);
            for (size_t i = 0; i < svec.index.size(); ++i) {
                output[svec.index[i]] = svec.value[i];
            }
            return output;
        }
    );
}

template<bool use_oracle_, typename Value_, typename Index_>
void test_block_access(
    const TestAccessParameters& params,
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref,
    int start,
    int end) 
{
    int nsecondary = (params.use_row ? ref->ncol() : ref->nrow());
    auto refwork = (params.use_row ? ref->dense_row() : ref->dense_column());
    test_access_base<use_oracle_>(
        params,
        ptr, 
        ref, 
        [&](int i) -> auto { 
            auto raw_expected = fetch(refwork.get(), i, nsecondary);
            return std::vector<Value_>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& svec) -> auto {
            std::vector<Value_> output(end - start);
            for (size_t i = 0; i < svec.value.size(); ++i) {
                output[svec.index[i] - start] = svec.value[i];
            }
            return output;
        },
        start,
        end - start
    );
}

template<bool use_oracle_, typename Value_, typename Index_>
void test_indexed_access(
    const TestAccessParameters& params, 
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref,
    int start,
    int step)
{
    int nsecondary = (params.use_row ? ref->ncol() : ref->nrow());
    auto refwork = (params.use_row ? ref->dense_row() : ref->dense_column());

    std::vector<Index_> indices;
    {
        int counter = start;
        while (counter < nsecondary) {
            indices.push_back(counter);
            counter += step;
        }
    }

    std::vector<size_t> reposition(nsecondary, -1);
    for (size_t i = 0, end = indices.size(); i < end; ++i) {
        reposition[indices[i]] = i;
    }

    test_access_base<use_oracle_>(
        params,
        ptr, 
        ref, 
        [&](int i) -> auto { 
            auto raw_expected = fetch(refwork.get(), i, nsecondary);
            std::vector<Value_> expected;
            expected.reserve(indices.size());
            for (auto idx : indices) {
                expected.push_back(raw_expected[idx]);
            }
            return expected;
        }, 
        [&](const auto& svec) -> auto {
            std::vector<Value_> expected(indices.size());
            for (size_t i = 0, end = svec.value.size(); i < end; ++i) {
                expected[reposition[svec.index[i]]] = svec.value[i];
            }
            return expected;
        },
        indices
    );
}

}

/********************************************************
 ********************************************************/

// Finally, some user-visible functions that convert compile-time parameters to run-time
// parameters, to make it easier to define parametrized tests across oracle/row status.
template<typename Value_, typename Index_>
void test_full_access(
    const TestAccessParameters& params, 
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref) 
{
    if (params.use_oracle) {
        internal::test_full_access<true>(params, ptr, ref);
    } else {
        internal::test_full_access<false>(params, ptr, ref);
    }
}

template<typename Value_, typename Index_>
void test_block_access(
    const TestAccessParameters& params,
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref,
    int start,
    int end) 
{
    if (params.use_oracle) {
        internal::test_block_access<true>(params, ptr, ref, start, end);
    } else {
        internal::test_block_access<false>(params, ptr, ref, start, end);
    }
}

template<typename Value_, typename Index_>
void test_indexed_access(
    const TestAccessParameters& params,
    const tatami::Matrix<Value_, Index_>* ptr, 
    const tatami::Matrix<Value_, Index_>* ref,
    int start,
    int step)
{
    if (params.use_oracle) {
        internal::test_indexed_access<true>(params, ptr, ref, start, step);
    } else {
        internal::test_indexed_access<false>(params, ptr, ref, start, step);
    }
}

/********************************************************
 ********************************************************/

// A couple of very simple helpers for just iterating through the matrix.
template<class Matrix_, class Matrix2_>
void test_simple_column_access(const Matrix_* ptr, const Matrix2_* ref) {
    TestAccessParameters params;
    params.use_row = false;
    test_full_access(params, ptr, ref);
}

template<class Matrix_, class Matrix2_>
void test_simple_row_access(const Matrix_* ptr, const Matrix2_* ref) {
    TestAccessParameters params;
    params.use_row = true;
    test_full_access(params, ptr, ref);
}

}

#endif
