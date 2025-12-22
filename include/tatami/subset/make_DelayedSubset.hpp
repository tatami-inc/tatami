#ifndef TATAMI_MAKE_DELAYED_SUBSET_HPP
#define TATAMI_MAKE_DELAYED_SUBSET_HPP

#include "DelayedSubsetSortedUnique.hpp"
#include "DelayedSubsetSorted.hpp"
#include "DelayedSubsetUnique.hpp"
#include "DelayedSubset.hpp"
#include "DelayedSubsetBlock.hpp"
#include "../utils/ArrayView.hpp"
#include "../utils/copy.hpp"

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
 * This will automatically dispatch to `DelayedSubsetSortedUnique`, `DelayedSubsetUnique`, `DelayedSubsetSorted` or `DelayedSubset`, depending on the values in `subset`.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam SubsetStorage_ Vector containing the subset indices, to be automatically deduced.
 * Any class implementing `[`, `size()`, `begin()` and `end()` can be used here.
 *
 * @param matrix Pointer to a (possibly `const`) `Matrix`.
 * @param subset Instance of the subset index vector.
 * @param by_row Whether to apply the subset to the rows.
 * If false, the subset is applied to the columns.
 *
 * @return A pointer to a `DelayedSubset` instance.
 */
template<typename Value_, typename Index_, class SubsetStorage_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > matrix, SubsetStorage_ subset, const bool by_row) {
    const auto nsub = subset.size();
    typedef I<decltype(nsub)> Subset;

    bool is_unsorted = false;
    for (Subset i = 1; i < nsub; ++i) {
        if (subset[i] < subset[i-1]) {
            is_unsorted = true;
            break;
        }
    }

    if (!is_unsorted) {
        bool has_duplicates = false;
        for (Subset i = 1; i < nsub; ++i) {
            if (subset[i] == subset[i-1]) {
                has_duplicates = true;
                break;
            }
        }

        if (!has_duplicates) {
            bool consecutive = true;
            for (Subset i = 1; i < nsub; ++i) {
                if (subset[i] > subset[i-1] + 1) {
                    consecutive = false;
                    break;
                }
            }

            if (consecutive) {
                const auto start = (nsub ? subset[0] : 0);
                return std::shared_ptr<Matrix<Value_, Index_> >(
                    new DelayedSubsetBlock<Value_, Index_>(std::move(matrix), start, subset.size(), by_row)
                );
            } else {
                return std::shared_ptr<Matrix<Value_, Index_> >(
                    new DelayedSubsetSortedUnique<Value_, Index_, SubsetStorage_>(std::move(matrix), std::move(subset), by_row, false)
                );
            }
        } else {
            return std::shared_ptr<Matrix<Value_, Index_> >(
                new DelayedSubsetSorted<Value_, Index_, SubsetStorage_>(std::move(matrix), std::move(subset), by_row, false)
            );
        }
    }

    bool has_duplicates = false;
    auto accumulated = create_container_of_Index_size<std::vector<unsigned char> >(by_row ? matrix->nrow() : matrix->ncol());
    for (Subset i = 0; i < nsub; ++i) {
        auto& found = accumulated[subset[i]];
        if (found) {
            has_duplicates = true;
            break;
        } else {
            found = 1;
        }
    }

    if (!has_duplicates) {
        return std::shared_ptr<Matrix<Value_, Index_> >(
            new DelayedSubsetUnique<Value_, Index_, SubsetStorage_>(std::move(matrix), std::move(subset), by_row, false)
        );
    } else {
        return std::shared_ptr<Matrix<Value_, Index_> >(
            new DelayedSubset<Value_, Index_, SubsetStorage_>(std::move(matrix), std::move(subset), by_row)
        );
    }
}

/**
 * @cond
 */
template<typename Value_, typename Index_, class SubsetStorage_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubset(std::shared_ptr<Matrix<Value_, Index_> > matrix, SubsetStorage_ subset, const bool by_row) {
    return make_DelayedSubset<Value_, Index_, SubsetStorage_>(std::shared_ptr<const Matrix<Value_, Index_> >(std::move(matrix)), std::move(subset), by_row);
}

// Back-compatibility only.
template<int margin_, typename Value_, typename Index_, class SubsetStorage_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > matrix, SubsetStorage_ subset) {
    return make_DelayedSubset<Value_, Index_, SubsetStorage_>(std::move(matrix), std::move(subset), margin_ == 0);
}

template<int margin_, typename Value_, typename Index_, class SubsetStorage_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubset(std::shared_ptr<Matrix<Value_, Index_> > matrix, SubsetStorage_ subset) {
    return make_DelayedSubset<Value_, Index_, SubsetStorage_>(std::move(matrix), std::move(subset), margin_ == 0);
}
/**
 * @endcond
 */

}

#endif
