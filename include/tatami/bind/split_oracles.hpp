#ifndef TATAMI_DELAYED_BIND_SPLIT_ORACLES_HPP
#define TATAMI_DELAYED_BIND_SPLIT_ORACLES_HPP

#include <memory>
#include <vector>
#include <cassert>
#include <optional>

#include "sanisizer/sanisizer.hpp"

#include "../base/Oracle.hpp"
#include "../utils/ConsecutiveOracle.hpp"
#include "../utils/FixedOracle.hpp"
#include "../utils/copy.hpp"

#include "CompressedMapping.hpp"

namespace tatami {

namespace DelayedBind_internal {

template<typename Index_>
class OffsetOracleWrapper : public Oracle<Index_> {
public:
    OffsetOracleWrapper(std::shared_ptr<const Oracle<Index_> > ptr, PredictionIndex start, PredictionIndex number, Index_ offset) :
        my_ptr(std::move(ptr)),
        my_start(start),
        my_number(number),
        my_offset(offset)
    {}

    PredictionIndex total() const {
        return my_number;
    }

    Index_ get(PredictionIndex i) const {
        return my_ptr->get(i + my_start) - my_offset;
    }
private:
    std::shared_ptr<const Oracle<Index_> > my_ptr;
    PredictionIndex my_start, my_number;
    Index_ my_offset;
};

template<typename Index_>
class OffsetOracleWrapperWithExtras : public Oracle<Index_> {
public:
    OffsetOracleWrapperWithExtras(std::shared_ptr<const Oracle<Index_> > ptr, PredictionIndex start, PredictionIndex number, Index_ offset, std::vector<Index_> extras) :
        my_ptr(std::move(ptr)),
        my_start(start),
        my_number(number),
        my_offset(offset),
        my_extras(std::move(extras))
    {
        assert(sanisizer::sum<PredictionIndex>(my_number, my_extras.size()));
    }

    PredictionIndex total() const {
        // This addition should be safe under the assumption that this is just a subset of predictions from the parent oracle.
        return my_number + my_extras.size();
    }

    Index_ get(PredictionIndex i) const {
        if (sanisizer::is_greater_than_or_equal(i, my_number)) {
            return my_extras[i - my_number];
        } else {
            return my_ptr->get(i + my_start) - my_offset;
        }
    }
private:
    std::shared_ptr<const Oracle<Index_> > my_ptr;
    PredictionIndex my_start, my_number;
    Index_ my_offset;
    std::vector<Index_> my_extras;
};

template<typename Index_>
class ConsecutiveOracleWithExtras : public Oracle<Index_> {
public:
    ConsecutiveOracleWithExtras(Index_ start, Index_ number, std::vector<Index_> extras) :
        my_start(start),
        my_number(number),
        my_extras(std::move(extras))
    {
        assert(sanisizer::sum<PredictionIndex>(my_number, my_extras.size()));
    }

    PredictionIndex total() const {
        // This addition should be safe under the assumption that this is just a subset of predictions from the parent oracle.
        return my_number + my_extras.size();
    }

    Index_ get(PredictionIndex i) const {
        if (sanisizer::is_greater_than_or_equal(i, my_number)) {
            return my_extras[i - my_number];
        } else {
            return my_start + i;
        }
    }
private:
    Index_ my_start, my_number;
    std::vector<Index_> my_extras;
};

template<typename Index_>
struct SplitPredictions {
    bool consecutive_predictions = true;
    Index_ start_prediction = 0;
    Index_ number_of_predictions = 0;

    bool consecutive_calls = true;
    PredictionIndex start_call = 0;
    PredictionIndex number_of_calls = 0;

    std::vector<Index_> realized_predictions;
};

template<typename Index_>
std::vector<SplitPredictions<Index_> > split_predictions(
    const std::vector<Index_>& cumulative,
    const CompressedMapping<Index_>& mapping, 
    const Oracle<Index_>& oracle
) {
    const auto nmats = cumulative.size() - 1;
    auto predictions = create_container_of_Index_size<std::vector<SplitPredictions<Index_> > >(nmats); // nmats should fit in an Index_, so this call is legal.

    const auto ntotal = oracle.total();
    for (I<decltype(ntotal)> i = 0; i < ntotal; ++i) {
        const auto raw_prediction = oracle.get(i);
        const Index_ choice = mapping[raw_prediction];
        const auto pred = raw_prediction - cumulative[choice];

        auto& current = predictions[choice];
        if (current.consecutive_predictions) {
            if (current.number_of_predictions == 0) {
                current.start_prediction = pred;
                current.number_of_predictions = 1;
            } else if (current.number_of_predictions + current.start_prediction == pred) {
                ++current.number_of_predictions;
            } else {
                current.consecutive_predictions = false;
            }
        }

        if (current.consecutive_calls) {
            if (current.number_of_calls == 0) {
                current.start_call = i;
                current.number_of_calls = 1;
            } else if (current.number_of_calls + current.start_call == i) {
                ++current.number_of_calls;
            } else {
                current.consecutive_calls = false;
            }
        }

        // We only need to start recording the realized predictions when both runs have failed,
        // as we will always use the longer run to initialize an Oracle.
        if (!current.consecutive_predictions && !current.consecutive_calls) {
            current.realized_predictions.push_back(pred);
        }
    }

    return predictions;
}

// Check whether the predictions only touch a single component matrix.
template<typename Index_>
std::optional<Index_> uses_lone_component(const std::vector<SplitPredictions<Index_> >& split_pred) {
    std::optional<Index_> out;
    const auto nmats = split_pred.size();
    for (I<decltype(nmats)> i = 0; i < nmats; ++i) {
        if (split_pred[i].number_of_predictions) {
            if (out.has_value()) {
                out.reset();
                break;
            }
            out = i;
        }
    }
    return out;
}

template<typename Index_, class Initialize_>
void split_oracles(
    const std::vector<Index_>& cumulative, 
    const std::shared_ptr<const Oracle<Index_> >& ora, 
    std::vector<SplitPredictions<Index_> > split_pred,
    Initialize_ init
) {
    const auto nmats = cumulative.size() - 1;
    assert(sanisizer::is_equal(nmats, split_pred.size()));

    for (I<decltype(nmats)> x = 0; x < nmats; ++x) {
        auto& sp = split_pred[x];

        if (sp.consecutive_predictions) {
            if (sp.number_of_predictions) {
                init(x, std::make_shared<ConsecutiveOracle<Index_> >(sp.start_prediction, sp.number_of_predictions));
            }

        } else if (sp.consecutive_calls) {
            if (sp.number_of_calls) {
                init(x, std::make_shared<OffsetOracleWrapper<Index_> >(ora, sp.start_call, sp.number_of_calls, cumulative[x]));
            }

        } else if (sanisizer::is_greater_than(sp.number_of_predictions, sp.number_of_calls)) {
            init(x, std::make_shared<ConsecutiveOracleWithExtras<Index_> >(sp.start_prediction, sp.number_of_predictions, std::move(sp.realized_predictions)));

        } else if (sanisizer::is_less_than(sp.number_of_predictions, sp.number_of_calls)) {
            init(x, std::make_shared<OffsetOracleWrapperWithExtras<Index_> >(ora, sp.start_call, sp.number_of_calls, cumulative[x], std::move(sp.realized_predictions)));

        } else if (sp.number_of_predictions > 1) {
            // i.e., if the runs are equal and are not of length 1, we prefer the ConsecutiveOracleWithExtras as that doesn't need to copy 'ora'.
            init(x, std::make_shared<ConsecutiveOracleWithExtras<Index_> >(sp.start_prediction, sp.number_of_predictions, std::move(sp.realized_predictions)));

        } else {
            // Need to add back the first element from a length-1 run.
            sp.realized_predictions.insert(sp.realized_predictions.begin(), 1, sp.start_prediction);
            init(x, std::make_shared<FixedVectorOracle<Index_> >(std::move(sp.realized_predictions)));
        }
    }
}

}

}

#endif
