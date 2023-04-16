#ifndef TATAMI_OPTIONS_HPP
#define TATAMI_OPTIONS_HPP

#include <vector>

/**
 * @file Options.hpp
 *
 * @brief Options for data access.
 */

namespace tatami {

/**
 * Type of selection, see `DimensionSelection::type` for details.
 */
enum class DimensionSelectionType : char { FULL, BLOCK, INDEX };

/**
 * @brief Select dimension elements of interest.
 *
 * Select the elements of interest in a dimension of a `Matrix`.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct DimensionSelection {
    /**
     * Selection type for this `DimensionSelect` instance.
     *
     * - `ALL`: selects the full extent of the dimension, i.e., all elements in the dimension.
     * - `BLOCK`: selects a contiguous block of elements in the dimension.
     * - `INDEX`: selects a sorted and unique array of indices of dimension elements.
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
     * Pointer to an array containing sorted and unique indices for dimension elements.
     * Only relevant if `type = DimensionSelectionType::INDEX`.
     * This is intended for use by `Matrix` developers where the array is guaranteed to outlive the `DimensionSelection` object.
     */
    const Index_* index_start = NULL;

    /**
     * Length of the array pointed to by `index_start`.
     * Only relevant if `type = DimensionSelectionType::INDEX` and `index_start != NULL`.
     */
    Index_ index_length = 0;

    /**
     * Vector containing sorted and unique indices for dimension elements.
     * Only used if `type = DimensionSelectionType::INDEX` and `index_start = NULL`.
     * This is provided to allow callers to transfer ownership of the array, if `index_start` would otherwise be pointing to a temporary array.
     */
    std::vector<Index_> indices;
};

/**
 * Access pattern for the elements of the iteration dimension,
 * see `IterationOptions::access_pattern` for more details.
 */
enum class AccessPattern : char { CONSECUTIVE, SEQUENCE, RANDOM };

/**
 * @brief Options for accessing elements along the iteration dimension.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct IterationOptions {
    /**
     * Selection of elements to be accessed in the iteration dimension.
     */
    DimensionSelection<Index_> selection;

    /** 
     * Whether to cache information from every call to `DenseFormat::fetch()` or `SparseFormat::fetch()`.
     * Specifically, this refers to intermediate element-specific values that can be re-used if the same dimension element is requested in a subsequent call.
     * This may enable faster iteration if the same `ExtractFormat` object is re-used for multiple passes over the same matrix.
     */
    bool cache_for_reuse = false;

   /**
     * Access pattern of elements on the iteration dimension.
     *
     * For `CONSECUTIVE`, elements are accessed in consecutive order.
     * This enables pre-fetching of the subsequent elements by `Matrix` implementations.
     * It is expected that calls to `DenseExtractor::fetch()` or `SparseExtractor::fetch()` involve consecutive `i` across the selection defined in `selection`.
     * Repeated calls to the same `i` are allowed.
     * Methods should wrap around to the start of the selection once the end of the selection is reached.
     * 
     * For `SEQUENCE`, elements are accessed in an _a priori_ known sequence.
     * The enables pre-fetching of the (not necessarily consecutive) next elements by `Matrix` implementations.
     * It is expected that `DenseExtractor::fetch()` or `SparseExtractor::fetch()` are called with `i` in the same order as `IterationOptions::sequence` or `IterationOptions::sequence_start`.
     * Repeated calls to the same `i` are allowed.
     * Any calls to the first `i` specified in the sequence should reset the `Extractor` to its initial position.
     * Methods should wrap around to the start of the selection once the end of the sequence is reached.
     *
     * For `RANDOM`: elements are accessed in a random sequence.
     * This allows `Matrix` implementations to dispense with any pre-fetching or caching.
     * Calls to `DenseExtractor::fetch()` or `SparseExtractor::fetch()` are expected to use `i` that lie within the selection defined in `selection`.
     */
    AccessPattern access_pattern = AccessPattern::CONSECUTIVE;

    /**
     * Pointer to an array containing the indices of elements to access. 
     * Contents should be a subset of elements in `selection`.
     * Only relevant if `access_order = AccessPattern::SEQUENCE`.
     * This is intended for use by `Matrix` developers where the array is guaranteed to outlive the `IterationOptions` object.
     */
    const Index_* sequence_start = NULL;

    /**
     * Length of the array pointed to by `index_start`.
     * Only relevant if `access_order = AccessPattern::SEQUENCE` and `sequence_start != NULL`.
     */
    Index_ sequence_length = 0;

    /**
     * Vector containing the sequence of indices for elements on the iteration dimension.
     * Contents should be a subset of elements in `selection`.
     * Only relevant if `access_order = AccessPattern::SEQUENCE` and `sequence_start != NULL`.
     * This is provided to allow callers to transfer ownership of the array, if `sequence_start` would otherwise be pointing to a temporary array.
     */
    std::vector<Index_> sequence;
};

/**
 * @brief Options for accessing elements along the extraction dimension.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct ExtractionOptions {
    /**
     * Selection of elements to be extracted in the extraction dimension.
     */
    DimensionSelection<Index_> selection;

    /** 
     * Whether to extract the sparse indices.
     * If set to `false`, this can be used to avoid unnecessary computation and copying in `Matrix` methods.
     * Only used in the sparse `DimensionAccess` classes.
     */
    bool sparse_extract_index = true;

    /** 
     * Whether to extract the sparse values.
     * If set to `false`, this can be used to avoid unnecessary computation and copying in `Matrix` methods.
     * Only used in the sparse `DimensionAccess` classes.
     */
    bool sparse_extract_value = true;

    /**
     * Whether the extracted sparse output should be ordered by increasing index.
     * Setting this to `false` may reduce computational work in situations where the order of non-zero elements does not matter.
     * Only used in the sparse `DimensionAccess` classes.
     */
    bool sparse_ordered_index = true;
};

}

#endif
