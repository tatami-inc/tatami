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
 * Type of limit on the dimension - none, a contiguous block, or a sorted and unique array of indices.
 */
enum class DimensionLimitType : char { NONE, BLOCK, INDEX };

/**
 * @brief Limits on the access for a dimension.
 *
 * Specify the dimension elements to be accessed in a given `DimensionAccess` instance.
 *
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct DimensionLimit {
    /**
     * Type of limit.
     */
    DimensionLimitType type = DimensionLimitType::NONE;

    /**
     * Index of the start of the contiguous block of elements.
     * Only relevant if `type = DimensionLimitType::BLOCK`.
     */
    Index_ block_start = 0;

    /**
     * Length of the contiguous block of elements.
     * Only relevant if `type = DimensionLimitType::BLOCK`.
     */
    Index_ block_length = 0;

    /**
     * Pointer to an array containing sorted and unique indices for dimension elements.
     * Only relevant if `type = DimensionLimitType::INDEX`.
     * This is intended for use by `Matrix` developers where the array is guaranteed to outlive the `DimensionLimit` object.
     */
    const Index_* index_start = NULL;

    /**
     * Length of the array pointed to by `index_start`.
     * Only relevant if `type = DimensionLimitType::INDEX` and `index_start != NULL`.
     */
    Index_ index_length = 0;

    /**
     * Vector containing sorted and unique indices for dimension elements.
     * Only used if `type = DimensionLimitType::INDEX` and `index_start = NULL`.
     * This is provided to allow callers to transfer ownership of the array, if `index_start` would otherwise be pointing to a temporary array.
     */
    std::vector<Index_> indices;
};

/**
 * Access pattern for the elements of the iteration dimension.
 * 
 * - `CONSECUTIVE`: consecutive elements are accessed.
 *   This is the most typical iteration strategy through a matrix and likely to be the most efficient for `Matrix` implementations.
 * - `SEQUENCE`: elements are accessed in an _a priori_ known sequence.
 *   The sequence itself should be available in `IterationOptions` to enable pre-fetching of the (not necessarily consecutive) next elements by `Matrix` implementations.
 * - `RANDOM`: elements are accessed in a random sequence.
 *   This allows `Matrix` implementations to dispense with any pre-fetching or caching.
 */
enum class AccessOrder : char { CONSECUTIVE, SEQUENCE, RANDOM };

template<typename Index_>
struct IterationOptions {
    /**
     * Limits on the elements to be accessed in the iteration dimension.
     */
    DimensionLimit<Index_> limit;

    /** 
     * Whether to cache information from every call to `DenseFormat::fetch()` or `SparseFormat::fetch()`.
     * Specifically, this refers to intermediate element-specific values that can be re-used if the same dimension element is requested in a subsequent call.
     * This may enable faster iteration if the same `ExtractFormat` object is re-used for multiple passes over the same matrix.
     */
    bool cache_for_reuse = false;

   /**
     * Expected access order of elements on the iteration dimension.
     * This may be used by implementations to optimize their extraction.
     */
    AccessOrder access_order = AccessOrder::CONSECUTIVE;

    /**
     * Pointer to an array containing the indices of elements to access. 
     * Only relevant if `access_order = AccessOrder::SEQUENCE`.
     * This is intended for use by `Matrix` developers where the array is guaranteed to outlive the `IterationOptions` object.
     */
    const Index_* sequence_start = NULL;

    /**
     * Length of the array pointed to by `index_start`.
     * Only relevant if `access_order = AccessOrder::SEQUENCE` and `sequence_start != NULL`.
     */
    Index_ sequence_length = 0;

    /**
     * Vector containing the sequence of indices for elements on the iteration dimension.
     * Only relevant if `access_order = AccessOrder::SEQUENCE` and `sequence_start != NULL`.
     * This is provided to allow callers to transfer ownership of the array, if `sequence_start` would otherwise be pointing to a temporary array.
     */
    std::vector<Index_> sequence;
};

template<typename Index_>
struct ExtractionOptions {
    /**
     * Limits on the elements to be accessed in the extraction dimension.
     */
    DimensionLimit<Index_> limit;

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
