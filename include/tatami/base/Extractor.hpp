#ifndef TATAMI_EXTRACTOR_HPP
#define TATAMI_EXTRACTOR_HPP

#include <vector>
#include <type_traits>
#include "SparseRange.hpp"
#include "Options.hpp"
#include "Oracle.hpp"

/**
 * @file Extractor.hpp
 *
 * @brief Virtual classes for extracting matrix data.
 */

namespace tatami {

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * @brief Extract a dimension element in dense form without an oracle.
 */
template<typename Value_, typename Index_>
struct MyopicDenseExtractor {
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param i Index of the desired dimension element,
     * i.e., the row or column index for instances created with `Matrix::dense_row()` or `Matrix::dense_column()`, respectively.
     * @param[out] buffer Pointer to an array of length no less than `N`, where `N` is defined as:
     * - the number of indices, for instances created by any `Matrix::dense_row()` or `Matrix::dense_column()` method that accepts a vector of indices.
     * - the block length, for instances created by any `Matrix::dense_row()` or `Matrix::dense_column()` method that accepts a block start and length.
     * - `Matrix::ncol()` or `Matrix::nrow()`, for all other instances created with `Matrix::dense_row()` and `Matrix::dense_column()` respectively.
     *
     * @return Pointer to an array containing the values from the `i`-th dimension element.
     * This is guaranteed to hold `N` values.
     */
    virtual const Value_* fetch(Index_ i, Value_* buffer) = 0;

    /**
     * @cond
     */
    ~MyopicDenseExtractor() = default;
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract in dense form with an oracle.
 */
template<typename Value_, typename Index_>
struct OracularDenseExtractor {
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param[out] buffer Pointer to an array of length no less than `N`,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     *
     * @return Pointer to an array containing the contents of the next dimension element,
     * as predicted by the `Oracle` used to construct this instance.
     * This is guaranteed to have `N` values.
     */
    virtual const Value_* fetch(Value_* buffer) = 0;

    /**
     * @cond
     */
    ~OracularDenseExtractor() = default;
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract in sparse form without an oracle.
 */
template<typename Value_, typename Index_>
struct MyopicSparseExtractor {
    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled.
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * If `Options::sparse_extract_value` was set to `false` during construction of this instance,
     * `vbuffer` is ignored and `SparseRange::value` is set to `NULL` in the output.
     * Similarly, if `Options::sparse_extract_index` was set to `false` during construction of this instance,
     * `ibuffer` is ignored and `SparseRange::index` is set to `NULL` in the output.
     *
     * @param i Index of the desired dimension element,
     * i.e., the row or column index for instances created with `Matrix::sparse_row()` or `Matrix::sparse_column()`, respectively.
     * @param[out] vbuffer Pointer to an array with enough space for at least `N` values, where `N` is defined as:
     * - the number of indices, for instances created by any `Matrix::sparse_row()` or `Matrix::sparse_column()` method that accepts a vector of indices.
     * - the block length, for instances created by any `Matrix::sparse_row()` or `Matrix::sparse_column()` method that accepts a block start and length.
     * - `Matrix::ncol()` or `Matrix::nrow()`, for all other instances created with `Matrix::sparse_row()` and `Matrix::sparse_column()` respectively.
     * @param[out] ibuffer Pointer to an array with enough space for at least `N` indices,
     * where `N` is defined as described for `vbuffer`
     *
     * @return A `SparseRange` object describing the contents of the `i`-th dimension element.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) = 0;

    /**
     * @cond
     */
    ~MyopicSparseExtractor() = default;
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract in sparse form with an oracle.
 */
template<typename Value_, typename Index_>
struct OracularSparseExtractor {
    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled.
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * If `Options::sparse_extract_value` was set to `false` during construction of this instance,
     * `vbuffer` is ignored and `SparseRange::value` is set to `NULL` in the output.
     * Similarly, if `Options::sparse_extract_index` was set to `false` during construction of this instance,
     * `ibuffer` is ignored and `SparseRange::index` is set to `NULL` in the output.
     *
     * @param[out] i Index of the desired element on the iteration dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `N` values,
     * where `N` is defined as described for `MyopicSparseExtractor::fetch()`.
     * @param[out] ibuffer Pointer to an array with enough space for at least `N` indices,
     * where `N` is defined as described for `MyopicSparseExtractor::fetch()`.
     *
     * @return A `SparseRange` object describing the contents of the next dimension element,
     * as predicted by the `Oracle` used to construct this instance.
     */
    virtual SparseRange<Value_, Index_> fetch(Value_* vbuffer, Index_* ibuffer) = 0;

    /**
     * @cond
     */
    ~OracularSparseExtractor() = default;
    /**
     * @endcond
     */
};

}

#endif
