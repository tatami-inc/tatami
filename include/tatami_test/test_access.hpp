#ifndef TATAMI_TEST_ACCESS_BASE_HPP
#define TATAMI_TEST_ACCESS_BASE_HPP

#include "utils.hpp"
#include <type_traits>
#include "../tatami/base/utils.hpp"
#include "../tatami/utils/ConsecutiveOracle.hpp"
#include "../tatami/utils/FixedOracle.hpp"

namespace tatami_test {

enum TestAccessOrder { FORWARD, REVERSE, RANDOM };

struct TestAccessParameters {
    // Whether to use an oracle.
    bool use_oracle = false;

    // Whether to test access across rows.
    bool use_row = true;

    // Whether to expect NaNs (and sanitize them in the comparison).
    bool has_nan = false;

    // Whether the input matrix is actually sparse and should have fewer values from the sparse extractors. 
    bool less_sparse = false; 

    // Ordering of rows/columns to test.
    TestAccessOrder order = FORWARD;

    // Minimum jump between rows/columns.
    int jump = 1;
};

/********************************************************
 ********************************************************/

// All functions in this namespace will ignore 'params.use_oracle' and
// 'params.use_row', as these are handled by template arguments for brevity.

namespace base {

template<bool use_row_, bool use_oracle_, class TestMatrix_, class RefMatrix_, class DenseExtract_, class SparseExpand_, class PropCheck_, typename ...Args_>
void test_access_base(const TestAccessParameters& params, const TestMatrix_* ptr, const RefMatrix_* ref, DenseExtract_ expector, SparseExpand_ sparse_expand, PropCheck_ check_properties, Args_... args) {
    int NR = ptr->nrow();
    ASSERT_EQ(NR, ref->nrow());
    int NC = ptr->ncol();
    ASSERT_EQ(NC, ref->ncol());
    int limit = (use_row_ ? NR : NC);

    // Setting up the sequence of observations.
    typedef typename TestMatrix_::index_type Index_;
    std::vector<Index_> sequence;
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

    // Setting up the oracle.
    std::shared_ptr<tatami::Oracle<Index_> > oracle;
    if constexpr(use_oracle_) {
        if (params.jump == 1 && params.order == FORWARD) {
            oracle.reset(new tatami::ConsecutiveOracle<Index_>(0, limit));
        } else {
            oracle.reset(new tatami::FixedViewOracle<Index_>(sequence.data(), sequence.size()));
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

    for (auto i : sequence) {
        auto expected = expector(i);
        sanitize_nan(expected, params.has_nan);
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
            sanitize_nan(observed, params.has_nan);
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
            sanitize_nan(observed.value, params.has_nan);
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
            sanitize_nan(values_only, params.has_nan);
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

            auto observed_uns = [&]() {
                if constexpr(use_oracle_) {
                    Index_ j = limit;
                    auto output = swork_uns->fetch(j);
                    EXPECT_EQ(i, j);
                    return output;
                } else {
                    return swork_uns->fetch(i);
                }
            }();
            sanitize_nan(observed_uns.value, params.has_nan);
            ASSERT_EQ(expected, sparse_expand(observed_uns));
        } 
    }

    // Check that sparse extraction isn't just the same as the dense extractor. 
    if (params.less_sparse) {
        if (ptr->sparse()) {
            ASSERT_TRUE(dense_extracted > sparse_extracted);
        }
    }

    // Checking the expected lengths.
    check_properties(pwork.get());
    check_properties(swork.get());
}

template<bool use_row_, bool use_oracle_, class Matrix_, class Matrix2_>
void test_full_access(const TestAccessParameters& params, const Matrix_* ptr, const Matrix2_* ref) {
    int nsecondary = (use_row_ ? ref->ncol() : ref->nrow());
    auto refwork = [&]() { 
        if constexpr(use_row_) {
            return ref->dense_row();
        } else {
            return ref->dense_column();
        }
    }();

    test_access_base<use_row_, use_oracle_>(
        params,
        ptr, 
        ref, 
        [&](int i) -> auto { 
            auto expected = refwork->fetch(i);
            EXPECT_EQ(expected.size(), nsecondary);
            return expected;
        },
        [&](const auto& range) -> auto {
            return expand(range, nsecondary);
        },
        [&](const auto* work) {
            EXPECT_EQ(work->full_length, nsecondary);
        }
    );
}

template<bool use_row_, bool use_oracle_, class Matrix_, class Matrix2_>
void test_block_access(const TestAccessParameters& params, const Matrix_* ptr, const Matrix2_* ref, int start, int end) {
    int nsecondary = (use_row_ ? ref->ncol() : ref->nrow());
    auto refwork = [&]() { 
        if constexpr(use_row_) {
            return ref->dense_row();
        } else {
            return ref->dense_column();
        }
    }();

    test_access_base<use_row_, use_oracle_>(
        params,
        ptr, 
        ref, 
        [&](int i) -> auto { 
            auto raw_expected = refwork->fetch(i);
            return std::vector<typename Matrix_::value_type>(raw_expected.begin() + start, raw_expected.begin() + end);
        }, 
        [&](const auto& range) -> auto {
            return expand(range, start, end);
        },
        [&](const auto* work) {
            EXPECT_EQ(work->block_start, start);
            EXPECT_EQ(work->block_length, end - start);
        },
        start,
        end - start
    );
}

template<bool use_row_, bool use_oracle_, class Matrix_, class Matrix2_>
void test_indexed_access(const TestAccessParameters& params, const Matrix_* ptr, const Matrix2_* ref, int start, int step) {
    int nsecondary = (use_row_ ? ref->ncol() : ref->nrow());
    auto refwork = [&]() { 
        if constexpr(use_row_) {
            return ref->dense_row();
        } else {
            return ref->dense_column();
        }
    }();

    std::vector<typename Matrix_::index_type> indices;
    {
        int counter = start;
        while (counter < nsecondary) {
            indices.push_back(counter);
            counter += step;
        }
    }

    test_access_base<use_row_, use_oracle_>(
        params,
        ptr, 
        ref, 
        [&](int i) -> auto { 
            auto raw_expected = refwork->fetch(i);
            std::vector<typename Matrix_::value_type> expected;
            expected.reserve(indices.size());
            for (auto idx : indices) {
                expected.push_back(raw_expected[idx]);
            }
            return expected;
        }, 
        [&](const auto& range) -> auto {
            auto full = expand(range, nsecondary);
            std::vector<double> sub;
            sub.reserve(indices.size());
            for (auto idx : indices) {
                sub.push_back(full[idx]);
            }
            return sub;
        },
        [&](const auto* work) {
            EXPECT_EQ(std::vector<int>(work->index_start(), work->index_start() + work->index_length), indices);
        },
        indices
    );
}

}

/********************************************************
 ********************************************************/

// Finally, some user-visible functions that convert compile-time parameters to run-time
// parameters, to make it easier to define parametrized tests across oracle/row status.
template<class Matrix_, class Matrix2_>
void test_full_access(const TestAccessParameters& params, const Matrix_* ptr, const Matrix2_* ref) {
    if (params.use_row) {
        if (params.use_oracle) {
            base::test_full_access<true, true>(params, ptr, ref);
        } else {
            base::test_full_access<true, false>(params, ptr, ref);
        }
    } else {
        if (params.use_oracle) {
            base::test_full_access<false, true>(params, ptr, ref);
        } else {
            base::test_full_access<false, false>(params, ptr, ref);
        } 
    }
}

template<class Matrix_, class Matrix2_>
void test_block_access(const TestAccessParameters& params, const Matrix_* ptr, const Matrix2_* ref, int start, int end) {
    if (params.use_row) {
        if (params.use_oracle) {
            base::test_block_access<true, true>(params, ptr, ref, start, end);
        } else {
            base::test_block_access<true, false>(params, ptr, ref, start, end);
        }
    } else {
        if (params.use_oracle) {
            base::test_block_access<false, true>(params, ptr, ref, start, end);
        } else {
            base::test_block_access<false, false>(params, ptr, ref, start, end);
        } 
    }
}

template<class Matrix_, class Matrix2_>
void test_indexed_access(const TestAccessParameters& params, const Matrix_* ptr, const Matrix2_* ref, int start, int step) {
    if (params.use_row) {
        if (params.use_oracle) {
            base::test_indexed_access<true, true>(params, ptr, ref, start, step);
        } else {
            base::test_indexed_access<true, false>(params, ptr, ref, start, step);
        }
    } else {
        if (params.use_oracle) {
            base::test_indexed_access<false, true>(params, ptr, ref, start, step);
        } else {
            base::test_indexed_access<false, false>(params, ptr, ref, start, step);
        } 
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