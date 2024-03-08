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
 *
 * We denote the "iteration" dimension of the `Matrix` as the one that is being iteratively accessed across its elements.
 * The other dimension is subsequently denoted as the "extraction" dimension.
 * For example, when iterating across the rows of a matrix, the rows are the iteration dimension, the columns are the extraction dimension,
 * and the current element of the iteration dimension is the specific row that is being accessed at any loop iteration.
 */

namespace tatami {

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract in dense form without an oracle.
 */
template<typename Value_, typename Index_>
struct MyopicDenseExtractor {
    /**
     * @return Number of elements to extract at each iteration. 
     */
    virtual Index_ number() const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param i Index of the desired element on the iteration dimension.
     * @param[out] buffer Pointer to an array of length no less than `number()`.
     *
     * @return Pointer to an array containing the values from the `i`-th element of the iteration dimension.
     * This is guaranteed to hold `number()` values.
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
     * @return Number of elements to extract at each iteration. 
     */
    virtual Index_ number() const = 0;

    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param[out] i Index of the predicted element on the iteration dimension.
     * @param[out] buffer Pointer to an array of length no less than `number()`.
     * This does not have to be filled if a pointer can be directly returned.
     *
     * @return Pointer to an array containing the `i`-th element of the iteration dimension.
     * This is guaranteed to have `number()` values.
     * `i` is filled with the index of the current prediction.
     */
    virtual const Value_* fetch(Index_& i, Value_* buffer) = 0;

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
     * @return Number of elements to extract at each iteration. 
     */
    virtual Index_ number() const = 0;

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * @param i Index of the desired element on the iteration dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `number()` values.
     * Ignored if `Options::sparse_extract_value` was set to `false` during construction of this instance.
     * @param[out] ibuffer Pointer to an array with enough space for at least `number()` indices.
     * Ignored if `Options::sparse_extract_index` was set to `false` during construction of this instance.
     *
     * @return A `SparseRange` object describing the contents of the desired dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     * Both arrays are guaranteed to have no more than `number()` values.
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
     * @return Number of elements to extract at each iteration. 
     */
    virtual Index_ number() const = 0;

    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * @param[out] i Index of the desired element on the iteration dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `number()` values.
     * Ignored if `Options::sparse_extract_value` was set to `false` during construction of this instance.
     * @param[out] ibuffer Pointer to an array with enough space for at least `number()` indices.
     * Ignored if `Options::sparse_extract_index` was set to `false` during construction of this instance.
     *
     * @return A `SparseRange` object describing the contents of the desired dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     * Both arrays are guaranteed to have no more than `number()` values.
     * `i` is filled with the index of the current prediction.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) = 0;

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
