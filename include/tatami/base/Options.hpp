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
     * This is intended for use by `Matrix` developers where the array is guaranteed to outlive the `DimensionLimits` object.
     */
    const Index_* index_start = NULL;

    /**
     * Length of the array pointed to by `index_start`.
     * Only relevant if `type = DimensionLimitType::INDEX`.
     */
    Index_ index_length = 0;

    /**
     * Vector containing sorted and unique indices for dimension elements.
     * Only used if `type = DimensionLimitType::INDEX` and `index_start = NULL`.
     *
     * This is provided to allow callers to transfer ownership of the array.
     * Users should use `indices` rather than set `index_start` to avoid problems with invalidation of the latter.
     */
    std::vector<Index_> indices;
};

/**
 * Access pattern for the elements of a matrix dimension.
 * This can be sorted (i.e., non-decreasing) or random.
 */
enum class AccessOrder : char { SORTED, RANDOM };

/**
 * @brief Options for data extraction.
 * 
 * Specify how the data is to be extracted.
 */
struct ExtractOptions {
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

    /** 
     * Whether to cache information from every call to `DimensionAccess::fetch()`.
     * Specifically, this refers to intermediate element-specific values that can be re-used if the same dimension element is requested in a subsequent call.
     * This may enable faster iteration if the same `DimensionAccess` object is re-used for multiple passes over the same matrix.
     */
    bool cache_for_reuse = false;

    /**
     * Expected access order of dimension elements.
     * This may be used by implementations to optimize their extraction.
     *
     * If set to `SORTED`, successive calls to `SparseFormat::fetch()` or `DenseFormat::fetch()` should have non-decreasing `i`.
     * Any call with a lower `i` can be assumed to restart iteration from an earlier element in the dimension.
     *
     * If set to `RANDOM`, successive calls to `SparseFormat::fetch()` or `DenseFormat::fetch()` may have any `i`.
     */
    AccessOrder access_order = AccessOrder::SORTED;
};

}

#endif
