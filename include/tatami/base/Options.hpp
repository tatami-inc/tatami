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
 * Type of selection along a dimension:
 *
 * - `FULL`: selects the full extent of the dimension, i.e., all elements in the dimension.
 * - `BLOCK`: selects a contiguous block of elements in the dimension.
 * - `INDEX`: selects a sorted and unique array of indices of dimension elements.
 */
enum class DimensionSelectionType : char { FULL, BLOCK, INDEX };

/**
 * @brief Options to define dimension elements of interest.
 *
 * Select the elements of interest in dimension of a `Matrix`.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct DimensionSelectionOptions {
    /**
     * Selection type for this `DimensionSelect` instance.
     */
    DimensionSelectionType type = DimensionSelectionType::FULL;

    /**
     * Index of the start of the contiguous block of elements.
     * Only relevant if `type = DimensionSelectionType::BLOCK`.
     */
    Index_ block_start = 0;

    /**
     * Length of the contiguous block of elements.
     * Only relevant if `type = DimensionSelectionType::BLOCK`.
     */
    Index_ block_length = 0;

    /**
     * Vector containing sorted and unique indices for dimension elements.
     * Only used if `type = DimensionSelectionType::INDEX`. 
     */
    std::vector<Index_> indices;
};

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
 * Access pattern for the elements of the iteration dimension.
 * This can be used by implementations to pre-fetch data for subsequent requests to `DenseExtractor::fetch()` or `SparseExtractor::fetch()`.
 *
 * `CONSECUTIVE` hints to the implementation that elements on the iteration dimension are accessed in consecutive order.
 * Implementations may assume that a `fetch()` request for one element will (usually) be followed by requests for the subsequent element.
 * However, this is not a strict requirement, and calls to `fetch()` may still be non-consecutive, so this should be handled appropriately.
 * 
 * For `SEQUENCE`, elements are accessed in a deterministic sequence.
 * The enables pre-fetching of future elements by `fetch()` implementations, improving the efficiency of iteration through the matrix.
 * The sequence is defined by the `AccessPatternOptions::sequencer`, and it is expected that calls to `fetch()` are in exactly the same order as those returned by `SequenceOracle::predict`.
 *
 * For `RANDOM`: elements are accessed in a random sequence.
 * This allows `Matrix` implementations to dispense with any pre-fetching or caching.
  */
enum class AccessPatternType : char { CONSECUTIVE, SEQUENCE, RANDOM };

/**
 * @brief Options for managing the accessing pattern.
 *
 * This refers to accesses along the iteration dimension.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct AccessPatternOptions {
    /** 
     * Whether to ask implementations to cache information from every call to `DenseExtractor::fetch()` or `SparseExtractor::fetch()`.
     * Specifically, this refers to intermediate element-specific values that can be re-used if the same dimension element is requested in a subsequent call.
     * This may enable faster iteration if the same `Extractor` object is re-used for multiple passes over the same matrix.
     */
    bool cache_for_reuse = false;

   /**
     * Access pattern of elements on the iteration dimension.
     *
     * If `CONSECUTIVE`, implementations should assume the user will be calling `fetch()` on the elements defined in `Options::selection`.
     * For example, if an indexed subset of the extraction dimension is of interest, implementations should assume that each call will access the indices in order.
     *
     * If `SEQUENCE`, the elements of the sequence should be a subset of those defined by `Options::selection`.
     * The `sequencer` should also be non-`NULL`.
     *
     * If `RANDOM`, the requested elements should be a subset of those defined by `Options::selection`.
     */
    AccessPatternType pattern = AccessPatternType::CONSECUTIVE;

    std::shared_ptr<SequenceOracle<Index_> > sequencer;
};

/**
 * @brief Options for accessing elements along the extraction dimension.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct Options {
    /**
     *
     */
    DimensionSelectionOptions<Index_> selection;

    SparseExtractionOptions<Index_> sparse;

    AccessPatternOptions<Index_> access;
};

}

#endif
