#ifndef TATAMI_DELAYED_ISOMETRIC_OPERATION_UTILS_HPP
#define TATAMI_DELAYED_ISOMETRIC_OPERATION_UTILS_HPP

#include <type_traits>
#include "../utils/new_extractor.hpp"

namespace tatami {

namespace DelayedIsometricOperation_internal {

// Some SFINAE shenanigans for checking if zero handling depends on the row ID of the zero value.
template<class Operation_, typename = int>
struct has_zero_depends_on_row {
    static constexpr bool value = false;
};

template<class Operation_>
struct has_zero_depends_on_row<Operation_, decltype((void) std::declval<Operation_>().zero_depends_on_row(), 0)> {
    static constexpr bool value = true;
};

template<class Operation_>
bool zero_depends_on_row(const Operation_& op) {
    if constexpr(!has_zero_depends_on_row<Operation_>::value) {
        return false;
    } else {
        return op.zero_depends_on_row();
    }
}

// ... Checking if zero handling depends on the column ID of the zero value.
template<class Operation_, typename = int>
struct has_zero_depends_on_column {
    static constexpr bool value = false;
};

template<class Operation_>
struct has_zero_depends_on_column<Operation_, decltype((void) std::declval<Operation_>().zero_depends_on_column(), 0)> {
    static constexpr bool value = true;
};

template<class Operation_>
bool zero_depends_on_column(const Operation_& op) {
    if constexpr(!has_zero_depends_on_column<Operation_>::value) {
        return false;
    } else {
        return op.zero_depends_on_column();
    }
}

// ... Checking if non-zero processing depends on the row ID.
template<class Operation_, typename = int>
struct has_non_zero_depends_on_row {
    static constexpr bool value = false;
};

template<class Operation_>
struct has_non_zero_depends_on_row<Operation_, decltype((void) std::declval<Operation_>().non_zero_depends_on_row(), 0)> {
    static constexpr bool value = true;
};

template<class Operation_>
bool non_zero_depends_on_row(const Operation_& op) {
    if constexpr(!has_non_zero_depends_on_row<Operation_>::value) {
        return false;
    } else {
        return op.non_zero_depends_on_row();
    }
}

// ... Checking if non-zero processing depends on the column ID.
template<class Operation_, typename = int>
struct has_non_zero_depends_on_column {
    static constexpr bool value = false;
};

template<class Operation_>
struct has_non_zero_depends_on_column<Operation_, decltype((void) std::declval<Operation_>().non_zero_depends_on_column(), 0)> {
    static constexpr bool value = true;
};

template<class Operation_>
bool non_zero_depends_on_column(const Operation_& op) {
    if constexpr(!has_non_zero_depends_on_column<Operation_>::value) {
        return false;
    } else {
        return op.non_zero_depends_on_column();
    }
}

// Check if the oracle requires the row or column ID.
template<bool oracle_, class Operation_, typename Index_>
class MaybeOracleDepends {
public:
    MaybeOracleDepends(const MaybeOracle<oracle_, Index_>& oracle, const Operation_& op, bool row) {
        if constexpr(oracle_) {
            if constexpr(Operation_::is_basic) {
                my_oracle = oracle;
            } else if ([&]() {
                // Only storing if the oracle if we need the row/column index
                // on the target dimension to apply the operation.
                if (row) {
                    if (non_zero_depends_on_row(op)) {
                        return true;
                    } else if (!op.is_sparse() && zero_depends_on_row(op)) { // if it's sparse, zeros remain zero so they can't depend on the row.
                        return true;
                    }
                } else {
                    if (non_zero_depends_on_column(op)) {
                        return true;
                    } else if (!op.is_sparse() && zero_depends_on_column(op)) { // ditto for the columns.
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
            if constexpr(Operation_::is_basic) {
                return my_oracle->get(my_used++);
            } else if constexpr(
                // If none of these are present, we can assume that they're all
                // false, allowing us to elide the my_oracle=NULL check.
                has_zero_depends_on_row<Operation_>::value || 
                has_zero_depends_on_column<Operation_>::value ||
                has_non_zero_depends_on_row<Operation_>::value || 
                has_non_zero_depends_on_column<Operation_>::value) 
            {
                if (my_oracle) {
                    return my_oracle->get(my_used++);
                }
            }
        }
        return i;
    }

private:
    MaybeOracle<oracle_, Index_> my_oracle;
    typename std::conditional<oracle_, size_t, bool>::type my_used = 0;
};

template<class Operation_>
bool can_dense_expand(const Operation_& op, bool row) {
    if (!op.is_sparse()) { // as I said before: if it's sparse, zeros remain zero so they can't depend on either the row or column.
        if (row) {
            // If zero processing doesn't depend on the column identity, then
            // during row extraction, we can use a constant value to fill in the
            // zero values for the entire row. 
            return !zero_depends_on_column(op);
        } else {
            // Similarly, if zero processing doesn't depend on
            // columns, then during row extraction, we can fill in the
            // zero values across columns with a constant value.
            return !zero_depends_on_row(op);
        }
    }

    return true;
}

template<class Operation_>
bool needs_sparse_indices(const Operation_& op, bool row) {
    if (row) {
        // If we depend on columns, then we need column indices
        // when the rows are the target dimension.
        return non_zero_depends_on_column(op);
    } else {
        // Similarly, if we depend on the rows, then we need row
        // indices when the columns are the target dimension.
        return non_zero_depends_on_row(op);
    }
}

}

}

#endif
