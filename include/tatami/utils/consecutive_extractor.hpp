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
 * @tparam row_ Whether to perform extraction on rows.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column index.
 * @tparam Args_ Types of further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @param mat Matrix to iterate over.
 * @param iter_start Index of the first row/column of the iteration range.
 * @param iter_length Number of rows/columns in the iteration range.
 * @param args Further arguments to pass to `Matrix::dense_row` or `Matrix::dense_column`.
 *
 * @return An extractor for iteration over consecutive rows/columns in `[iter_start, iter_start + iter_length)`.
 * This may be either an `OracularDenseExtractor` or `OracularSparseExtractor` depending on `sparse_`.
 *
 * This function is equivalent to `new_extractor()` but creates an a `ConsecutiveOracle` instance.
 * `Matrix` implementations that are oracle-aware can then perform pre-fetching of future accesses for greater performance.
 * Of course, this assumes that the iteration over the target dimension does actually involve consecutive elements from `iter_start` to `iter_start + iter_length`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto consecutive_extractor(const Matrix<Value_, Index_>* mat, Index_ iter_start, Index_ iter_length, Args_&&... args) {
    return new_extractor<row_, sparse_>(mat, 
        std::make_shared<ConsecutiveOracle<Index_> >(iter_start, iter_length),
        std::forward<Args_>(args)...
    );
}

}

#endif
