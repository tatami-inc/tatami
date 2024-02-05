#ifndef TATAMI_TEST_ORACLE_ACCESS_HPP
#define TATAMI_TEST_ORACLE_ACCESS_HPP

#include <gtest/gtest.h>

#include "../tatami/base/Matrix.hpp"
#include "../tatami/utils/Oracles.hpp"
#include <random>
#include <deque>

namespace tatami_test {

template<class Matrix, bool row_, typename ... Args_>
void test_oracle_access(const Matrix* ptr, const Matrix* ref, bool randomized, Args_... args) {
    int NR = ptr->nrow();
    int NC = ptr->ncol();

    auto pwork = [&]() {
        if constexpr(row_) {
            return ref->dense_row(args...);
        } else {
            return ref->dense_column(args...);
        }
    }();
    auto swork = [&]() {
        if constexpr(row_) {
            return ref->sparse_row(args...);
        } else {
            return ref->sparse_column(args...);
        }
    }();

    typedef typename Matrix::index_type Index_;

    if (randomized) {
        std::mt19937_64 rng(NR + NC * 10); // making up an interesting seed.
        std::vector<Index_> fixed(NC * 2);
        for (auto& x : fixed) {
            x = rng() % NC;
        }

        auto oracle = std::make_shared<tatami::FixedOracle<Index_> >(fixed.data(), fixed.size());
        auto pwork_o = [&]() {
            if constexpr(row_) {
                return ptr->dense_row(oracle, args...);
            } else {
                return ptr->dense_column(oracle, args...);
            }
        }();
        auto swork_o = [&]() {
            if constexpr(row_) {
                return ptr->sparse_row(oracle, args...);
            } else {
                return ptr->sparse_column(oracle, args...);
            }
        }();

        for (auto i : fixed) {
            auto expected = pwork->fetch(i);
            auto sexpected = swork->fetch(i);

            Index_ jd;
            auto observed = pwork_o->fetch(jd);
            EXPECT_EQ(jd, i);
            EXPECT_EQ(expected, observed);

            Index_ js;
            auto sobserved = swork_o->fetch(js);
            EXPECT_EQ(js, i);
            EXPECT_EQ(sexpected.value, sobserved.value);
        }

    } else {
        auto oracle = std::make_unique<tatami::ConsecutiveOracle<Index_> >(0, NC);
        auto pwork_o = [&]() {
            if constexpr(row_) {
                return ptr->dense_row(oracle, args...);
            } else {
                return ptr->dense_column(oracle, args...);
            }
        }();
        auto swork_o = [&]() {
            if constexpr(row_) {
                return ptr->sparse_row(oracle, args...);
            } else {
                return ptr->sparse_column(oracle, args...);
            }
        }();

        for (Index_ i = 0; i < NC; ++i) {
            auto expected = pwork->fetch(i);
            auto sexpected = swork->fetch(i);

            Index_ jd;
            auto observed = pwork_o->fetch(jd);
            EXPECT_EQ(jd, i);
            EXPECT_EQ(expected, observed);

            Index_ js;
            auto sobserved = swork_o->fetch(js);
            EXPECT_EQ(js, i);
            EXPECT_EQ(sexpected.value, sobserved.value);
        }
    }
}

}

#endif
