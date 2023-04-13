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
 * @tparam MARGIN Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam MAT A specialized `Matrix`, to be automatically deducted.
 * @tparam V Vector containing the subset indices, to be automatically deducted.
 *
 * @param p Pointer to a `Matrix`.
 * @param idx Instance of the index vector.
 *
 * @return A pointer to a `DelayedSubset` instance.
 */
template<int MARGIN, class MAT, class V>
std::shared_ptr<MAT> make_DelayedSubset(std::shared_ptr<MAT> p, V idx) {
    bool is_unsorted = false;
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        if (i) {
            if (idx[i] < idx[i-1]) {
                is_unsorted = true;
                break;
            }
        }
    }

    if (!is_unsorted) {
        bool has_duplicates = false;
        for (size_t i = 0, end = idx.size(); i < end; ++i) {
            if (i) {
                if (idx[i] == idx[i-1]) {
                    has_duplicates = true;
                    break;
                }
            }
        }

        if (!has_duplicates) {
            return std::shared_ptr<MAT>(
                new DelayedSubsetSortedUnique<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), std::move(idx), false)
            );
        } else {
            return std::shared_ptr<MAT>(
                new DelayedSubsetSorted<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), std::move(idx), false)
            );
        }
    }

    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<V>()[0])>::type>::type V_type;
    std::vector<std::pair<V_type, size_t> > collected;
    collected.reserve(idx.size());
    for (size_t i = 0, end = idx.size(); i < end; ++i) {
        collected.emplace_back(idx[i], i);
    }
    std::sort(collected.begin(), collected.end());

    bool has_duplicates = false;
    for (size_t i = 1, end = collected.size(); i < end; ++i) {
        if (collected[i].first == collected[i-1].first) {
            has_duplicates = true;
            break;
        }
    }

    if (!has_duplicates) {
        return std::shared_ptr<MAT>(
            new DelayedSubsetUnique<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), collected, std::move(idx))
        );
    } else {
        return std::shared_ptr<MAT>(
            new DelayedSubset<MARGIN, typename MAT::data_type, typename MAT::index_type, typename std::remove_reference<V>::type>(std::move(p), collected, std::move(idx))
        );
    }
}

}

#endif
