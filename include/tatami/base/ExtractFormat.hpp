#ifndef TATAMI_EXTRACT_FORMAT_HPP
#define TATAMI_EXTRACT_FORMAT_HPP

#include <type_traits>
#include "SparseRange.hpp"
#include "Options.hpp"

/**
 * @file DimensionAccess.hpp
 *
 * @brief Virtual class for acessing matrix dimensions.
 *
 * We denote the "target" dimension of the `Matrix` as the one that is being iteratively accessed across its elements.
 * The other dimension is subsequently denoted as "non-target" dimension.
 * For example, when iterating across the rows of a matrix, the rows are the target dimension, and the columns are the non-target dimension.
 * The current dimension element is the specific row that is being accessed at any loop iteration.
 */

namespace tatami {

/**
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Virtual base class for extraction formatting.
 */
template<typename Index_>
struct FormatBase {
protected:
    /**
     * @cond
     */
    FormatBase() = default;
    virtual ~FormatBase() = default;
    /**
     * @endcond
     */

public:
    /**
     * Type of limit on the extraction for the elements on the non-target dimension, when accessing a single element of the target dimension.
     * For example, when extracting a particular row, a setting of `DimensionLimitType::NONE` indicates that entries were extracted across all columns for that row.
     */
    DimensionLimitType extracted_limit = DimensionLimitType::NONE;

    /**
     * Number of extracted entries of a dimension element.
     * 
     * - For `extracted_limit = DimensionLimitType::NONE`, this should be the total extent of the non-target dimension.
     * - For `extracted_limit = DimensionLimitType::BLOCK`, this should be the size of the contiguous block along the non-target dimension.
     * - For `extracted_limit = DimensionLimitType::INDEX`, this should be the number of indices for the non-target dimension.
     */
    Index_ extracted_length = 0;

    /**
     * Index of the start of the contiguous block of entries along the non-target dimension.
     * Only relevant if `extracted_limit = DimensionLimitType::BLOCK`.
     */
    Index_ extracted_block = 0;

    /**
     * @return Pointer to the start of an array of length `extracted_length`,
     * containing the sorted and unique entry indices along the non-target dimension.
     * Only relevant if `extracted_limit = DimensionLimitType::INDEX`.
     */
    virtual const Index_* extracted_index() const = 0;
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Virtual base class for dense extraction.
 */
template<typename Value_, typename Index_>
class DenseFormat : public FormatBase<Index_> {
protected:
    /**
     * @cond
     */
    DenseFormat() = default;
    /**
     * @endcond
     */

public:
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param i Index of the current element on the target dimension.
     * @param buffer Pointer to an array of length specified by `FormatBase::extracted_length`.
     *
     * @return Pointer to the contents of the current dimension element, containing `FormatBase::extracted_length` valid entries.
     */
    virtual const Value_* fetch(Index_ i, Value_* buffer) = 0;

    /**
     * @param i Index of the current element on the target dimension.
     * @param[out] buffer Pointer to an array of length specified by `FormatBase::extracted_length`.
     * This is filled with the contents of the current dimension element.
     *
     * @return `buffer`, for consistency with `fetch()`.
     */
    const Value_* fetch_copy(Index_ i, Value_* buffer) {
        auto out = fetch(i, buffer);
        if (out != buffer) {
            std::copy(out, out + this->extracted_length, buffer);
        }
        return buffer;
    }

    /**
     * @param i Index of the current element on the target dimension.
     * @return Vector of length `FormatBase::extracted_length`, containing the contents of the current dimension element.
     */
    std::vector<Value_> fetch(Index_ i) {
        std::vector<Value_> buffer(this->extracted_length);
        fetch_copy(i, buffer.data());
        return buffer;
    }

    /**
     * Whether this class enables sparse access.
     */
    static constexpr bool sparse_ = false;
};

/**
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Virtual base class for dense extraction.
 */
template<typename Value_, typename Index_>
class SparseFormat : public FormatBase<Index_> {
protected:
    /**
     * @cond
     */
    SparseFormat() = default;
    /**
     * @endcond
     */

public:
    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * @param i Index of the current element on the target dimension.
     * @param vbuffer Pointer to an array with enough space for at least `FormatBase::extracted_length` values.
     * Ignored if `WorkspaceOptions::sparse_extract_value` was set to `false` during construction of this instance.
     * @param ibuffer Pointer to an array with enough space for at least `FormatBase::extracted_length` indices.
     * Ignored if `WorkspaceOptions::sparse_extract_index` was set to `false` during construction of this instance.
     *
     * @return A `SparseRange` object describing the contents of the current dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `WorkspaceOptions::sparse_extract_mode` used to construct this object.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) = 0;

    /**
     * @param i Index of the current element on the target dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `FormatBase::extracted_length` values.
     * On output, this is filled with the values of the structural non-zeros. 
     *
     * Ignored if `WorkspaceOptions::sparse_extract_value` was set to `false` during construction of this instance.
     * Also ignored if set to `NULL`, in which case the values are extracted but not copied to `vbuffer`.
     * @param[out] ibuffer Pointer to an array with enough space for at least `FormatBase::extracted_length` indices.
     * On output, this is filled with the indices of the structural non-zeros. 
     *
     * Ignored if `WorkspaceOptions::sparse_extract_index` was set to `false` during construction of this instance.
     * Also ignored if set to `NULL`, in which case the indices are extracted but not copied to `ibuffer`.
     *
     * @return A `SparseRange` object describing the contents of the current dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `WorkspaceOptions::sparse_extract_mode` used to construct this object.
     */
    virtual SparseRange<Value_, Index_> fetch_copy(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto output = fetch(i, vbuffer, ibuffer);

        if (vbuffer != NULL) {
            if (output.value != NULL) {
                if (output.value != vbuffer) {
                    std::copy(output.value, output.value + output.number, vbuffer);
                    output.value = vbuffer;
                }
            }
        }

        if (ibuffer != NULL) {
            if (output.index != NULL) {
                if (output.index != ibuffer) {
                    std::copy(output.index, output.index + output.number, ibuffer);
                    output.index = ibuffer;
                }
            }
        }

        return output;
    }

    /**
     * @param i Index of the current element on the target dimension.
     * @return A `SparseRangeCopy` object containing the contents of the current dimension element.
     * Either or both of `value` or `index` is empty if extraction of that field is skipped, 
     * based on the setting of `WorkspaceOptions::sparse_extract_mode` used to construct this object.
     */
    SparseRangeCopy<Value_, Index_> fetch(Index_ i) {
        SparseRangeCopy<Value_, Index_> output(this->extracted_length);
        auto range = fetch_copy(i, output.value.data(), output.index.data());
        output.number = range.number;
        output.value.resize(range.value != NULL ? range.number : 0);
        output.index.resize(range.index != NULL ? range.number : 0);
        return output;
    }

    /**
     * Whether this class enables sparse access.
     */
    static constexpr bool sparse_ = true;
};

/**
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Conditional class definition for dense/sparse formats.
 * If `sparse = true`, this will have the `SparseFormat::fetch()` method, otherwise it will have the `DenseFormat::fetch()` method.
 */
template<bool sparse_, typename Value_, typename Index_>
using ExtractFormat = typename std::conditional<sparse_, SparseFormat<Value_, Index_>, DenseFormat<Value_, Index_> >::type;

}

#endif
