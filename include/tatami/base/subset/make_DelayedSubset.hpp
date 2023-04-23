#ifndef TATAMI_MAKE_DELAYED_SUBSET_HPP
#define TATAMI_MAKE_DELAYED_SUBSET_HPP

#include "DelayedSubsetSortedUnique.hpp"
#include "DelayedSubsetSorted.hpp"
#include "DelayedSubsetUnique.hpp"
#include "DelayedSubset.hpp"
#include "DelayedSubsetBlock.hpp"

#include <algorithm>
#include <memory>

/**
 * @file make_DelayedSubset.hpp
 *
 * @brief Make a delayed subset wrapper based on row/column indices.
 */

namespace tatami {

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 * This will automatically dispatch to `DelayedSubsetSortedUnique`, `DelayedSubsetUnique`, `DelayedSubsetSorted` or `DelayedSubset`, depending on the values in `idx`.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Matrix_ A specialized `Matrix`, to be automatically deducted.
 * @tparam IndexStsorage_ Vector containing the subset indices, to be automatically deduced.
 *
 * @param p Pointer to a `Matrix`.
 * @param idx Instance of the index vector.
 *
 * @return A pointer to a `DelayedSubset` instance.
 */
template<int margin_, class Matrix_, class IndexStorage_>
std::shared_ptr<Matrix_> make_DelayedSubset(std::shared_ptr<Matrix_> p, IndexStorage_ idx) {
    typedef typename Matrix_::value_type Value_;
    typedef typename Matrix_::index_type Index_;
    typedef typename std::remove_reference<IndexStorage_>::type PureIndexStorage_;
    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<IndexStorage_>()[0])>::type>::type storage_type;

    bool is_unsorted = false;
    for (Index_ i = 0, end = idx.size(); i < end; ++i) {
        if (i) {
            if (idx[i] < idx[i-1]) {
                is_unsorted = true;
                break;
            }
        }
    }

    if (!is_unsorted) {
        bool has_duplicates = false;
        for (Index_ i = 0, end = idx.size(); i < end; ++i) {
            if (i) {
                if (idx[i] == idx[i-1]) {
                    has_duplicates = true;
                    break;
                }
            }
        }

        if (!has_duplicates) {
            return std::shared_ptr<Matrix_>(
                new DelayedSubsetSortedUnique<margin_, Value_, Index_, PureIndexStorage_>(std::move(p), std::move(idx), false)
            );
        } else {
            return std::shared_ptr<Matrix_>(
                new DelayedSubsetSorted<margin_, Value_, Index_, PureIndexStorage_>(std::move(p), std::move(idx), false)
            );
        }
    }

    std::vector<std::pair<storage_type, Index_> > collected;
    collected.reserve(idx.size());
    for (Index_ i = 0, end = idx.size(); i < end; ++i) {
        collected.emplace_back(idx[i], i);
    }
    std::sort(collected.begin(), collected.end());

    bool has_duplicates = false;
    for (Index_ i = 1, end = collected.size(); i < end; ++i) {
        if (collected[i].first == collected[i-1].first) {
            has_duplicates = true;
            break;
        }
    }

    if (!has_duplicates) {
        return std::shared_ptr<Matrix_>(
            new DelayedSubsetUnique<margin_, Value_, Index_, PureIndexStorage_>(std::move(p), collected, std::move(idx))
        );
    } else {
        return std::shared_ptr<Matrix_>(
            new DelayedSubset<margin_, Value_, Index_, PureIndexStorage_>(std::move(p), collected, std::move(idx))
        );
    }
}

}

#endif
