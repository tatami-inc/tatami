#ifndef TATAMI_DELAYED_ISOMETRIC_OPERATION_UTILS_HPP
#define TATAMI_DELAYED_ISOMETRIC_OPERATION_UTILS_HPP

#include <type_traits>
#include <cstddef>

#include "../utils/new_extractor.hpp"

namespace tatami {

namespace DelayedIsometricOperation_internal {

// Check if the oracle requires the row or column ID.
template<bool oracle_, class Operation_, typename Index_>
class MaybeOracleDepends {
public:
    MaybeOracleDepends(const MaybeOracle<oracle_, Index_>& oracle, const Operation_& op, bool row) {
        if constexpr(oracle_) {
            if ([&]{
                // Only storing if the oracle if we need the row/column index
                // on the target dimension to apply the operation.
                if (row) {
                    if (op.non_zero_depends_on_row()) {
                        return true;
                    } else if (!op.is_sparse() && op.zero_depends_on_row()) { // if it's sparse, zeros remain zero so they can't depend on the row.
                        return true;
                    }
                } else {
                    if (op.non_zero_depends_on_column()) {
                        return true;
                    } else if (!op.is_sparse() && op.zero_depends_on_column()) { // ditto for the columns.
                        return true;
                    }
                }
                return false;
            }()) 
            {
                my_oracle = oracle;
            }
        }
    }

    Index_ get(Index_ i) {
        if constexpr(oracle_) {
            if (my_oracle) {
                return my_oracle->get(my_used++);
            }
        }
        return i;
    }

private:
    MaybeOracle<oracle_, Index_> my_oracle;
    typename std::conditional<oracle_, PredictionIndex, bool>::type my_used = 0;
};

template<class Operation_>
bool can_dense_expand(const Operation_& op, bool row) {
    if (!op.is_sparse()) { // as I said before: if it's sparse, zeros remain zero so they can't depend on either the row or column.
        if (row) {
            // If zero processing doesn't depend on the column identity, then
            // during row extraction, we can use a constant value to fill in the
            // zero values for the entire row. 
            return !op.zero_depends_on_column();
        } else {
            // Similarly, if zero processing doesn't depend on
            // columns, then during row extraction, we can fill in the
            // zero values across columns with a constant value.
            return !op.zero_depends_on_row();
        }
    }

    return true;
}

template<class Operation_>
bool needs_sparse_indices(const Operation_& op, bool row) {
    if (row) {
        // If we depend on columns, then we need column indices
        // when the rows are the target dimension.
        return op.non_zero_depends_on_column();
    } else {
        // Similarly, if we depend on the rows, then we need row
        // indices when the columns are the target dimension.
        return op.non_zero_depends_on_row();
    }
}

}

}

#endif
