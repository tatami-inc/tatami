#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <cstdint>
#include <random>
#include <algorithm>

#include "tatami/bind/split_oracles.hpp"

class DelayedBindSplitOraclesTest : public ::testing::TestWithParam<std::tuple<std::vector<int>, std::vector<int> > > {};

TEST_P(DelayedBindSplitOraclesTest, Basic) {
    std::vector<int> cumulative{ 0, 5, 13, 22 };

    auto params = GetParam();
    auto predictions = std::get<0>(params);
    auto types = std::get<1>(params);

    tatami::DelayedBind_internal::CompressedMapping<int> mapping(5, cumulative);
    auto ora = std::make_shared<const tatami::FixedVectorOracle<int> >(predictions);

    const std::size_t nmats = cumulative.size() - 1;
    std::vector<std::shared_ptr<const tatami::Oracle<int> > > split(nmats);
    tatami::DelayedBind_internal::split_oracles<int>(
        cumulative,
        ora,
        tatami::DelayedBind_internal::split_predictions(cumulative, mapping, *ora),
        [&](const int x, std::shared_ptr<const tatami::Oracle<int> > subora) -> void {
            split[x] = std::move(subora);
        }
    );

    std::vector<std::vector<int> > fragmented(cumulative.size());
    for (auto p : predictions) {
        auto chosen = mapping[p];
        fragmented[chosen].push_back(p - cumulative[chosen]);
    }

    for (std::size_t x = 0; x < nmats; ++x) {
        if (fragmented[x].size()) {
            ASSERT_TRUE(split[x]);
            const auto npreds = split[x]->total();
            EXPECT_EQ(fragmented[x].size(), npreds);
            for (tatami::PredictionIndex p = 0; p < npreds; ++p) {
                EXPECT_EQ(fragmented[x][p], split[x]->get(p));
            }
        } else {
            EXPECT_FALSE(split[x]);
            continue;
        }

        // Check that the correct classes are used.
        if (types[x] == 0) {
            EXPECT_TRUE(dynamic_cast<const tatami::ConsecutiveOracle<int>*>(split[x].get()));
        } else if (types[x] == 1) {
            EXPECT_TRUE(dynamic_cast<const tatami::DelayedBind_internal::OffsetOracleWrapper<int>*>(split[x].get()));
        } else if (types[x] == 2) {
            EXPECT_TRUE(dynamic_cast<const tatami::DelayedBind_internal::ConsecutiveOracleWithExtras<int>*>(split[x].get()));
        } else if (types[x] == 3) {
            EXPECT_TRUE(dynamic_cast<const tatami::DelayedBind_internal::OffsetOracleWrapperWithExtras<int>*>(split[x].get()));
        } else {
            EXPECT_TRUE(dynamic_cast<const tatami::FixedVectorOracle<int>*>(split[x].get()));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    DelayedBind,
    DelayedBindSplitOraclesTest,
    ::testing::Values(
        // Hard-coded boundaries are { 0, 5, 13, 22 }

        // simple consecutive predictions
        std::make_tuple((std::vector<int>{ 1,2,3, 8,9,10,11,12, 15,16 }), (std::vector<int>{ 0, 0, 0 })),

        // interspersed consecutive predictions
        std::make_tuple((std::vector<int>{ 8, 1, 9, 2, 15, 16, 10, 3, 17, 11, 4 }), (std::vector<int>{ 0, 0, 0 })),

        // consecutive calls, non-consecutive predictions.
        std::make_tuple((std::vector<int>{ 6,8,10,12, 4,3,2,1,0, 20,18,19,14,15,21,13 }), (std::vector<int>{ 1, 1, 1 })),

        // partial consecutive predictions, interspersed so that runs of consecutive predictions are longer than the runs of consecutive calls.
        std::make_tuple((std::vector<int>{ 5,6, 1, 7,6,5, 13,14,15, 2,4, 16,17,19,21 }), (std::vector<int>{ 2, 2, 2 })),

        // partial consecutive predictions, interspersed so that runs of consecutive predictions are equal to the runs of consecutive calls.
        std::make_tuple((std::vector<int>{ 5,6,7, 1,2, 6,5, 13,14,15,16,17, 4, 19,21 }), (std::vector<int>{ 2, 2, 2 })),

        // partial consecutive predictions, interspersed so that runs of consecutive predictions are less than the runs of consecutive calls.
        std::make_tuple((std::vector<int>{ 5,6,7,6,6, 1,2,4, 8,5, 13,14,15,16,17,18,20, 4,3, 19 }), (std::vector<int>{ 3, 3, 3 })),

        // fully interspersed non-consecutive predictions.
        std::make_tuple((std::vector<int>{ 1, 10, 20, 3, 8, 15, 4, 2, 6, 11, 20, 5, 15, 2 }), (std::vector<int>{ 4, 4, 4 })),

        // mixed (I)
        std::make_tuple((std::vector<int>{ 0, 10, 1,2, 8, 20,20,15,12, 3, 7,6, 4, }), (std::vector<int>{ 0, 4, 1 })),

        // mixed (II)
        std::make_tuple((std::vector<int>{ 20,19,17,16,15, 1,2,2,1, 7,8,9,10,11, 3,4, 13,16,20, 6,5 }), (std::vector<int>{ 3, 2, 3 })),

        // some empty
        std::make_tuple((std::vector<int>{ 7,8,9,10,12,11 }), (std::vector<int>{ 0, 1, 0 }))
    )
);

TEST(DelayedBind, UsesLoneComponent) {
    std::vector<int> cumulative{ 0, 5, 13, 22 };
    tatami::DelayedBind_internal::CompressedMapping<int> mapping(5, cumulative);

    {
        std::vector<int> predictions{ 7,8,9,10,12,11 };
        auto split = tatami::DelayedBind_internal::split_predictions(cumulative, mapping, tatami::FixedVectorOracle<int>(std::move(predictions)));
        auto comp = tatami::DelayedBind_internal::uses_lone_component(split);
        EXPECT_TRUE(comp.has_value());
        EXPECT_EQ(*comp, 1);
    }

    {
        std::vector<int> predictions{ 7,8,9,20,12,11 };
        auto split = tatami::DelayedBind_internal::split_predictions(cumulative, mapping, tatami::FixedVectorOracle<int>(std::move(predictions)));
        auto comp = tatami::DelayedBind_internal::uses_lone_component(split);
        EXPECT_FALSE(comp.has_value());
    }

    {
        auto split = tatami::DelayedBind_internal::split_predictions(cumulative, mapping, tatami::FixedVectorOracle<int>({}));
        auto comp = tatami::DelayedBind_internal::uses_lone_component(split);
        EXPECT_FALSE(comp.has_value());
    }
}
