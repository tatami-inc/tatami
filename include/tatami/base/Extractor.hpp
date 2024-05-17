#ifndef TATAMI_EXTRACTOR_HPP
#define TATAMI_EXTRACTOR_HPP

#include <vector>
#include <type_traits>
#include "SparseRange.hpp"
#include "Options.hpp"

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
 * @brief Extract an element of the target dimension in dense form without an oracle.
 */
template<typename Value_, typename Index_>
class MyopicDenseExtractor {
public:
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param i Index of the target dimension element, i.e., the row or column index.
     * @param[out] buffer Pointer to an array of length no less than `N`, where `N` is defined as:
     * - the number of columns, when extracting each row.
     * - the number of rows, when extracting each column.
     * - the block length, when extracting a contiguous block from each row/column.
     * - the number of indices, when extracting an indexed subset of each row/column.
     *
     * @return Pointer to an array containing the values from the `i`-th dimension element.
     * This is guaranteed to hold `N` values.
     */
    virtual const Value_* fetch(Index_ i, Value_* buffer) = 0;

    /**
     * @cond
     */
    MyopicDenseExtractor() = default;
    MyopicDenseExtractor(const MyopicDenseExtractor&) = default;
    MyopicDenseExtractor& operator=(const MyopicDenseExtractor&) = default;
    MyopicDenseExtractor(MyopicDenseExtractor&&) = default;
    MyopicDenseExtractor& operator=(MyopicDenseExtractor&&) = default;
    virtual ~MyopicDenseExtractor() = default;
    /**
     * @endcond
     */

    /**
     * @cond
     */
    // No-op for back-compatibility only.
    template<class Oracle_>
    void set_oracle(Oracle_) {}
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract an element of the target dimension in dense form with an oracle.
 */
template<typename Value_, typename Index_>
class OracularDenseExtractor {
public:
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param[out] buffer Pointer to an array of length no less than `N`,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     *
     * @return Pointer to an array containing the contents of the next element of the target dimension,
     * as predicted by the `Oracle` used to construct this instance.
     * This is guaranteed to have `N` values.
     */
    const Value_* fetch(Value_* buffer) {
        return fetch(0, buffer);
    }

    /**
     * This overload is intended for developers only.
     * It introduces the `i` argument so that the signature is the same as that of `MyopicDenseExtractor::fetch()`.
     * This makes it easier to define `MyopicDenseExtractor` and `OracularDenseExtractor` subclasses from a single template,
     * avoiding code duplication that would otherwise occur when defining methods with and without `i`.
     * Of course, implementations are expected to ignore `i` in oracle-aware extraction.
     *
     * Other than the extra `i` argument, all other behaviors of the two overloads are the same.
     * To avoid confusion, most users should just use the `fetch()` overload that does not accept `i`,
     * given that the value of `i` is never actually used.
     *
     * @param i Ignored, only provided for consistency with `MyopicDenseExtractor::fetch()`,
     * @param[out] buffer Pointer to an array of length no less than `N`,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     *
     * @return Pointer to an array containing the contents of the next dimension element,
     * as predicted by the `Oracle` used to construct this instance.
     * This is guaranteed to have `N` values.
     */
    virtual const Value_* fetch(Index_ i, Value_* buffer) = 0;

    /**
     * @cond
     */
    OracularDenseExtractor() = default;
    OracularDenseExtractor(const OracularDenseExtractor&) = default;
    OracularDenseExtractor& operator=(const OracularDenseExtractor&) = default;
    OracularDenseExtractor(OracularDenseExtractor&&) = default;
    OracularDenseExtractor& operator=(OracularDenseExtractor&&) = default;
    virtual ~OracularDenseExtractor() = default;
    /**
     * @endcond
     */

    /**
     * @cond
     */
    // No-op for back-compatibility only.
    template<class Oracle_>
    void set_oracle(Oracle_) {}
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract an element of the target dimension in sparse form without an oracle.
 */
template<typename Value_, typename Index_>
class MyopicSparseExtractor {
public:
    /**
     * `value_buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `value_buffer`;
     * if they are the same, `value_buffer` has been filled with `SparseRange::number` values.
     * The same applies for `index_buffer` and the returned `SparseRange::index` pointer.
     *
     * If `Options::sparse_extract_value` was set to `false` during construction of this instance,
     * `value_buffer` is ignored and `SparseRange::value` is set to `NULL` in the output.
     * Similarly, if `Options::sparse_extract_index` was set to `false` during construction of this instance,
     * `index_buffer` is ignored and `SparseRange::index` is set to `NULL` in the output.
     *
     * @param i Index of the target dimension element, i.e., the row or column index.
     * @param[out] value_buffer Pointer to an array with enough space for at least `N` values,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     * @param[out] index_buffer Pointer to an array with enough space for at least `N` indices,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     *
     * @return A `SparseRange` object describing the contents of the `i`-th dimension element.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) = 0;

    /**
     * @cond
     */
    MyopicSparseExtractor() = default;
    MyopicSparseExtractor(const MyopicSparseExtractor&) = default;
    MyopicSparseExtractor& operator=(const MyopicSparseExtractor&) = default;
    MyopicSparseExtractor(MyopicSparseExtractor&&) = default;
    MyopicSparseExtractor& operator=(MyopicSparseExtractor&&) = default;
    virtual ~MyopicSparseExtractor() = default;
    /**
     * @endcond
     */

    /**
     * @cond
     */
    // No-op for back-compatibility only.
    template<class Oracle_>
    void set_oracle(Oracle_) {}
    /**
     * @endcond
     */
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Extract an element of the target dimension in sparse form with an oracle.
 */
template<typename Value_, typename Index_>
class OracularSparseExtractor {
public:
    /**
     * `value_buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `value_buffer`; 
     * if they are the same, `value_buffer` has been filled with `SparseRange::number` values.
     * The same applies for `index_buffer` and the returned `SparseRange::index` pointer.
     *
     * If `Options::sparse_extract_value` was set to `false` during construction of this instance,
     * `value_buffer` is ignored and `SparseRange::value` is set to `NULL` in the output.
     * Similarly, if `Options::sparse_extract_index` was set to `false` during construction of this instance,
     * `index_buffer` is ignored and `SparseRange::index` is set to `NULL` in the output.
     *
     * @param[out] value_buffer Pointer to an array with enough space for at least `N` values,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     * @param[out] index_buffer Pointer to an array with enough space for at least `N` indices,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     *
     * @return A `SparseRange` object describing the contents of the next element of the target dimension, 
     * as predicted by the `Oracle` used to construct this instance.
     */
    SparseRange<Value_, Index_> fetch(Value_* value_buffer, Index_* index_buffer) {
        return fetch(0, value_buffer, index_buffer);
    }

    /**
     * This overload is intended for developers only.
     * It introduces the `i` argument so that the signature is the same as that of `MyopicSparseExtractor::fetch()`.
     * This makes it easier to define `MyopicSparseExtractor` and `OracularSparseExtractor` subclasses from a single template,
     * avoiding code duplication that would otherwise occur when defining methods with and without `i`.
     * Of course, implementations are expected to ignore `i` in oracle-aware extraction.
     *
     * Other than the extra `i` argument, all other behaviors of the two overloads are the same.
     * To avoid confusion, most users should just use the `fetch()` overload that does not accept `i`,
     * given that the value of `i` is never actually used.
     *
     * @param i Ignored, only provided for consistency with `MyopicSparseExtractor::fetch()`.
     * @param[out] value_buffer Pointer to an array with enough space for at least `N` values,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     * @param[out] index_buffer Pointer to an array with enough space for at least `N` indices,
     * where `N` is defined as described for `MyopicDenseExtractor::fetch()`.
     *
     * @return A `SparseRange` object describing the contents of the next dimension element,
     * as predicted by the `Oracle` used to construct this instance.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) = 0;

    /**
     * @cond
     */
    OracularSparseExtractor() = default;
    OracularSparseExtractor(const OracularSparseExtractor&) = default;
    OracularSparseExtractor& operator=(const OracularSparseExtractor&) = default;
    OracularSparseExtractor(OracularSparseExtractor&&) = default;
    OracularSparseExtractor& operator=(OracularSparseExtractor&&) = default;
    virtual ~OracularSparseExtractor() = default;
    /**
     * @endcond
     */

    /**
     * @cond
     */
    // No-op for back-compatibility only.
    template<class Oracle_>
    void set_oracle(Oracle_) {}
    /**
     * @endcond
     */
};

/**
 * @tparam oracle_ Whether to use an oracle-aware interface.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Type alias that switches between `OracularDenseExtractor` and `MyopicDenseExtractor` depending on `oracle_`.
 * Intended for templated class definitions, where setting `oracle_` can define subclasses for both interfaces.
 */
template<bool oracle_, typename Value_, typename Index_>
using DenseExtractor = typename std::conditional<oracle_, OracularDenseExtractor<Value_, Index_>, MyopicDenseExtractor<Value_, Index_> >::type;

/**
 * @tparam oracle_ Whether to use an oracle-aware interface.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Type alias that switches between `OracularSparseExtractor` and `MyopicSparseExtractor` depending on `oracle_`.
 * Intended for templated class definitions, where setting `oracle_` can define subclasses for both interfaces.
 */
template<bool oracle_, typename Value_, typename Index_>
using SparseExtractor = typename std::conditional<oracle_, OracularSparseExtractor<Value_, Index_>, MyopicSparseExtractor<Value_, Index_> >::type;

}

#endif
