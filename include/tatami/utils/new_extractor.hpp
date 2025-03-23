#ifndef TATAMI_NEW_EXTRACTOR_HPP
#define TATAMI_NEW_EXTRACTOR_HPP

#include "../base/Matrix.hpp"

/**
 * @file new_extractor.hpp
 * @brief Templated construction of a new extractor.
 */

namespace tatami {

/**
 * @tparam oracle_ Whether an oracle should be supplied.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Maybe an `Oracle`, maybe a placeholder boolean, depending on `oracle_`.
 */
template<bool oracle_, typename Index_>
using MaybeOracle = typename std::conditional<oracle_, std::shared_ptr<const Oracle<Index_> >, bool>::type;

/**
 * This utility makes it easier for developers to write a single templated function that works with and without oracles.
 * A boolean placeholder should be provided as the "oracle" in the myopic extractor case.
 * 
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam oracle_ Whether an oracle should be supplied.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @tparam Args_ Further arguments.
 *
 * @param[in] matrix A `tatami::Matrix` object to iterate over.
 * @param row Whether to create a row-wise extractor, i.e., the rows are the target dimension.
 * @param oracle Pointer to an oracle if `oracle_ = true`, otherwise a placeholder boolean that is ignored.
 * @param args Zero or more additional arguments to pass to methods like `Matrix::dense_row()`.
 *
 * @return An extractor to access elements of the requested dimension of `ptr`.
 * This may be any of `MyopicDenseExtractor`, `MyopicSparseExtractor`, `OracularDenseExtractor` or `OracularSparseExtractor`,
 * depending on `sparse_` and `oracle_`.
 */
template<bool sparse_, bool oracle_, typename Value_, typename Index_, typename ... Args_>
auto new_extractor(const Matrix<Value_, Index_>& matrix, bool row, MaybeOracle<oracle_, Index_> oracle, Args_&&... args) {
    // We could use 'sparse()' and 'dense()' directly here, but that assumes
    // that 'opt' is always supplied at the end of 'args', which might not be
    // the case... so we just spell it out and save ourselves the trouble
    // of implementing default option methods for 'sparse()' and 'dense()'.
    if constexpr(sparse_) {
        if constexpr(oracle_) {
            if (row) {
                return matrix.sparse_row(std::move(oracle), std::forward<Args_>(args)...);
            } else {
                return matrix.sparse_column(std::move(oracle), std::forward<Args_>(args)...);
            }
        } else {
            if (row) {
                return matrix.sparse_row(std::forward<Args_>(args)...);
            } else {
                return matrix.sparse_column(std::forward<Args_>(args)...);
            }
        }
    } else {
        if constexpr(oracle_) {
            if (row) {
                return matrix.dense_row(std::move(oracle), std::forward<Args_>(args)...);
            } else {
                return matrix.dense_column(std::move(oracle), std::forward<Args_>(args)...);
            }
        } else {
            if (row) {
                return matrix.dense_row(std::forward<Args_>(args)...);
            } else {
                return matrix.dense_column(std::forward<Args_>(args)...);
            }
        }
    }
}

/**
 * @cond
 */
// Provided for back-compatibility only.
template<bool sparse_, bool oracle_, typename Value_, typename Index_, typename ... Args_>
auto new_extractor(const Matrix<Value_, Index_>* ptr, bool row, tatami::MaybeOracle<oracle_, Index_> oracle, Args_&&... args) {
    return new_extractor<sparse_, oracle_, Value_, Index_>(*ptr, row, std::move(oracle), std::forward<Args_>(args)...);
}

template<bool row_, bool sparse_, typename Value_, typename Index_>
auto new_extractor(const Matrix<Value_, Index_>* ptr) {
    return new_extractor<sparse_, false, Value_, Index_>(ptr, row_, false);
}

template<bool row_, bool sparse_, typename Value_, typename Index_>
auto new_extractor(const Matrix<Value_, Index_>* ptr, const Options& opt) {
    return new_extractor<sparse_, false, Value_, Index_>(ptr, row_, false, opt);
}

template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto new_extractor(const Matrix<Value_, Index_>* ptr, Index_ start, Index_ len, Args_&&... args) {
    return new_extractor<sparse_, false, Value_, Index_>(ptr, row_, false, start, len, std::forward<Args_>(args)...);
}

template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
auto new_extractor(const Matrix<Value_, Index_>* ptr, std::vector<Index_> i, Args_&&... args) {
    return new_extractor<sparse_, false, Value_, Index_>(ptr, row_, false, std::move(i), std::forward<Args_>(args)...);
}
/**
 * @endcond
 */

}

#endif
