#ifndef TATAMI_OPTIONS_HPP
#define TATAMI_OPTIONS_HPP

#include <vector>
#include <memory>
#include "SequenceOracle.hpp"

/**
 * @file Options.hpp
 *
 * @brief Options for data access.
 */

namespace tatami {

/**
 * Type of selection along a dimension:
 *
 * - `FULL`: selects the full extent of the dimension, i.e., all elements in the dimension.
 * - `BLOCK`: selects a contiguous block of elements in the dimension.
 * - `INDEX`: selects a sorted and unique array of indices of dimension elements.
 */
enum class DimensionSelectionType : char { FULL, BLOCK, INDEX };

/**
 * @brief Options for sparse extraction.
 */
struct SparseExtractionOptions {
    /** 
     * Whether to extract the indices in `SparseExtractor::fetch`.
     * If set to `false`, this can be used to avoid unnecessary computation and copying in `Matrix` methods.
     */
    bool extract_index = true;

    /** 
     * Whether to extract the values in `SparseExtractor::fetch`.
     * If set to `false`, this can be used to avoid unnecessary computation and copying in `Matrix` methods.
     */
    bool extract_value = true;

    /**
     * Whether the structural non-zeros returned by `SparseExtractor::fetch` should be ordered by increasing index.
     * Setting this to `false` may reduce computational work in situations where the order of non-zero elements does not matter.
     */
    bool ordered_index = true;
};

/**
 * @brief Options for managing iteration.
 *
 * This refers to access to elements along the iteration dimension,
 * across multiple calls to `DenseExtractor::fetch()` or `SparseExtractor::fetch()`.
 */
struct IterationOptions {
    /** 
     * Whether to ask implementations to cache information from every 
     * Specifically, this refers to intermediate element-specific values that can be re-used if the same dimension element is requested in a subsequent call.
     * This may enable faster iteration if the same `Extractor` object is re-used for multiple passes over the same matrix.
     */
    bool cache_for_reuse = false;
};

/**
 * @brief Options for accessing elements. 
 */
struct Options {
    /**
     * Options for sparse extraction.
     */
    SparseExtractionOptions sparse;

    /**
     * Options to define the access pattern on the iteration dimension.
     */
    IterationOptions<Index_> access;
};

}

#endif
