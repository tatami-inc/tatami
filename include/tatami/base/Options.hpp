#ifndef TATAMI_OPTIONS_HPP
#define TATAMI_OPTIONS_HPP

#include <vector>
#include <memory>

/**
 * @file Options.hpp
 *
 * @brief Options for data access.
 */

namespace tatami {

/**
 * Selection of dimension elements: 
 *
 * - `FULL`: selects the full extent of the dimension, i.e., all elements. 
 * - `BLOCK`: selects a contiguous block of elements in the dimension.
 * - `INDEX`: selects a sorted and unique array of indices of dimension elements.
 *
 * This corresponds to the various overloads of `Matrix` methods (e.g., `Matrix::dense_row()`),
 * each of which extracts a different selection of elements from the non-target dimension.
 */
enum class DimensionSelectionType : char { FULL, BLOCK, INDEX };

/**
 * @brief Options for accessing data from a `Matrix` instance. 
 */
struct Options {
    /** 
     * Whether to extract the indices of structural non-zero elements in `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * If set to `false`, this can be used to avoid unnecessary computation and copying. 
     * 
     * Note that the number of structural non-zeros reported by `fetch()` should be independent of this setting. 
     * This means that `Matrix` implementations should not try to do something overly clever that might cause the results to change depending on whether the indices are extracted -
     * for example, only reporting the structural non-zeros with odd indices -
     * though admittedly this warning is more relevant for `sparse_extract_value`.
     */
    bool sparse_extract_index = true;

    /** 
     * Whether to extract the values of structural non-zero elements in `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor::fetch()`.
     * If set to `false`, this can be used to avoid unnecessary computation and copying. 
     * 
     * Note that the number of structural non-zeros reported by `fetch()` should be independent of this setting. 
     * This means that `Matrix` implementations should not try to do something overly clever when reporting results,
     * e.g., like filtering out structural non-zeros with values of zero when `sparse_extract_value = true`;
     * doing so implies that the same filtering should be performed when `sparse_extract_value = false`,
     * which requires extraction of the value (and thus defeats the purpose of this flag).
     */
    bool sparse_extract_value = true;

    /**
     * Whether the structural non-zeros returned by `MyopicSparseExtractor::fetch()` or `OracularSparseExtractor:fetch()` should be ordered by increasing index.
     * Setting this to `false` may reduce computational work in situations where the order of non-zero elements does not matter.
     * 
     * The number and identity of structural non-zeros reported by `fetch()` should be independent of this setting, only the order is allowed to vary.
     * Note that the order may even vary across separate calls to `fetch()` for the same extractor, e.g., when pulling data parts asynchronously from a remote source.
     * However, the number and identity of structural non-zeros in each call should not change.
     */
    bool sparse_ordered_index = true;
};

}

#endif
