#ifndef TATAMI_CONSEUCTIVE_EXTRACTOR_HPP
#define TATAMI_CONSEUCTIVE_EXTRACTOR_HPP

#include <memory>
#include "../base/Matrix.hpp"
#include "new_extractor.hpp"
#include "ConsecutiveOracle.hpp"

/**
 * @file consecutive_extractor.hpp
 * @brief Templated construction of a new consecutive extractor.
 */

namespace tatami {

/**
 * This function creates an extractor object with a `ConsecutiveOracle` instance spanning a range of rows or columns.
 * `Matrix` implementations that are oracle-aware can then perform pre-fetching of future accesses for greater performance.
 *
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column index.
 * @tparam Args_ Types of further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @param matrix A `tatami::Matrix` to iterate over.
 * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension. 
 * @param iter_start Index of the first row (if `row = true`) or column (otherwise) of the iteration range.
 * @param iter_length Number of rows (if `row = true`) or columns (otherwise) in the iteration range.
 * @param args Further arguments to pass to `new_extractor()`.
 *
 * @return An extractor for iteration over consecutive rows/columns in `[iter_start, iter_start + iter_length)`.
 * This may be either an `OracularDenseExtractor` or `OracularSparseExtractor` depending on `sparse_`.
 */
template<bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>& matrix, const bool row, const Index_ iter_start, const Index_ iter_length, Args_&&... args) {
    return new_extractor<sparse_, true>(
        matrix,
        row,
        std::make_shared<ConsecutiveOracle<Index_> >(iter_start, iter_length),
        std::forward<Args_>(args)...
    );
}

/**
 * @cond
 */
// Provided for back-compatibility only.
template<bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>* matrix, bool row, Index_ iter_start, Index_ iter_length, Args_&&... args) {
    return consecutive_extractor<sparse_>(*matrix, row, iter_start, iter_length, std::forward<Args_>(args)...);
}

template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>* ptr, Args_&&... args) {
    return consecutive_extractor<sparse_, Value_, Index_>(ptr, row_, std::forward<Args_>(args)...);
}
/**
 * @endcond
 */

}

#endif
