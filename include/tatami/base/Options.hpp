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
 * Type of selection along a dimension, typically the extraction dimension:
 *
 * - `FULL`: selects the full extent of the dimension, i.e., all elements in the dimension.
 * - `BLOCK`: selects a contiguous block of elements in the dimension.
 * - `INDEX`: selects a sorted and unique array of indices of dimension elements.
 */
enum class DimensionSelectionType : char { FULL, BLOCK, INDEX };

/**
 * @brief Options for iteration and extraction.
 */
struct Options {
    /** 
     * Whether to extract the indices in `SparseExtractor::fetch`.
     * If set to `false`, this can be used to avoid unnecessary computation and copying in `Matrix` methods.
     */
    bool sparse_extract_index = true;

    /** 
     * Whether to extract the values in `SparseExtractor::fetch`.
     * If set to `false`, this can be used to avoid unnecessary computation and copying in `Matrix` methods.
     */
    bool sparse_extract_value = true;

    /**
     * Whether the structural non-zeros returned by `SparseExtractor::fetch` should be ordered by increasing index.
     * Setting this to `false` may reduce computational work in situations where the order of non-zero elements does not matter.
     */
    bool sparse_ordered_index = true;

    /** 
     * Whether to ask implementations to cache information from every 
     * Specifically, this refers to intermediate element-specific values that can be re-used if the same dimension element is requested in a subsequent call.
     * This may enable faster iteration if the same `Extractor` object is re-used for multiple passes over the same matrix.
     */
    bool cache_for_reuse = false;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future access requests.
 *
 * This allows `Matrix` implementations to pre-fetch data for future requests to `DenseExtractor::fetch()` or `SparseExtractor::fetch()`.
 */
template<typename Index_>
struct Oracle {
    /**
     * @cond
     */
    virtual ~Oracle() = default;
    /**
     * @endcond
     */

    /**
     * Predict the indices to be accessed in future `fetch()` calls.
     *
     * @param[out] predicted Pointer to an array in which to store the predicted indices of future elements to be accessed by `fetch()`.
     * @param number Maximum number of indices to predict.
     *
     * @return Number of indices that were predicted.
     * This is guaranteed to be no greater than `number`.
     * If no more predictions are available, this method should return zero.
     */
    virtual size_t predict(Index_* predicted, size_t number) = 0;
};

}

#endif
