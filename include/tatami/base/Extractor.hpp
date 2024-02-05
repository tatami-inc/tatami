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
 * @brief Virtual class for extracting matrix data.
 *
 * We denote the "iteration" dimension of the `Matrix` as the one that is being iteratively accessed across its elements.
 * The other dimension is subsequently denoted as the "extraction" dimension.
 * For example, when iterating across the rows of a matrix, the rows are the iteration dimension, the columns are the extraction dimension,
 * and the current element of the iteration dimension is the specific row that is being accessed at any loop iteration.
 */

namespace tatami {

/**
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Virtual base class for full access.
 * 
 * This is an interface class that provides access to the full extent of the extraction dimension.
 */
template<typename Index_>
struct FullExtractor {
    /**
     * Full extent of the extraction dimension.
     */
    Index_ full_length = 0;

    /**
     * @cond
     */
    virtual ~FullExtractor() = default;
    /**
     * @endcond
     */
};

/**
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Virtual base class for block access.
 *
 * This is an interface class that provides access to a contiguous block of elements along the extraction dimension.
 */
template<typename Index_>
struct BlockExtractor {
    /**
     * Index of the start of the contiguous block of entries along the extraction dimension.
     */
    Index_ block_start = 0;

    /**
     * Size of the contiguous block along the extraction dimension.
     */
    Index_ block_length = 0;

    /**
     * @cond
     */
    virtual ~BlockExtractor() = default;
    /**
     * @endcond
     */
};

/**
 * @tparam Index_ Row/column index type, should be integer.
 * @brief Virtual base class for indexed access. 
 *
 * This is an interface class that provides access to a subset of elements on the extraction dimension,
 * where the subset is defined by a vector of sorted and unique indices.
 */
template<typename Index_> 
struct IndexExtractor {
    /**
     * Unlike `index_length`, this is implemented as a virtual method to avoid invalidation of pointers when `IndexExtractor` instances are copied or moved.
     *
     * @return Pointer to an array containing the sorted and unique entry indices along the extraction dimension.
     */
    virtual const Index_* index_start() const = 0;

    /**
     * @return Length of the array pointed to by `indices()`.
     */
    Index_ index_length = 0;

    /**
     * @cond
     */
    virtual ~IndexExtractor() = default;
    /**
     * @endcond
     */
};


/**
 * @tparam selection_ Type of selection along the extraction dimension.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Conditional extractor interface that depends on the selection type.
 */
template<DimensionSelectionType selection_, typename Index_> 
using ConditionalSelectionExtractor = typename std::conditional<
        selection_ == DimensionSelectionType::FULL,
        FullExtractor<Index_>,
        typename std::conditional<
                selection_ == DimensionSelectionType::BLOCK,
                BlockExtractor<Index_>,
                IndexExtractor<Index_>
            >::type
    >::type;

/**
 * @tparam selection_ Type of selection along the extraction dimension.
 * @tparam Index_ Row/column index type, should be integer.
 * 
 * @param ex A `ConditionalSelectionExtractor` object.
 * @return Number of elements extracted from the extraction dimension, conditional on the selection type in `selection_`.
 */
template<DimensionSelectionType selection_, typename Index_> 
Index_ extracted_length(const ConditionalSelectionExtractor<selection_, Index_>& ex) {
    if constexpr(selection_ == DimensionSelectionType::FULL) {
        return ex.full_length;
    } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
        return ex.block_length;
    } else {
        return ex.index_length;
    }
}

/*********************************************************
 *********************************************************/

/**
 * @tparam selection_ Type of selection along the extraction dimension.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * @brief Virtual base class for dense extraction.
 */
template<DimensionSelectionType selection_, typename Value_, typename Index_>
class DenseExtractor : public ConditionalSelectionExtractor<selection_, Index_> {
protected:
    /**
     * @cond
     */
    DenseExtractor() = default;
    /**
     * @endcond
     */

public:
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * @param i Index of the desired element on the iteration dimension.
     * @param buffer Pointer to an array of length no less than `extracted_length()`.
     *
     * @return Pointer to the contents of the `i`-th element of the iteration dimension, containing `extracted_length()` entries.
     */
    virtual const Value_* fetch(Index_ i, Value_* buffer) = 0;

    /**
     * @param i Index of the desired element on the iteration dimension.
     * @param[out] buffer Pointer to an array of length specified by `extracted_length()`.
     * This is filled with the contents of the desired dimension element.
     *
     * @return `buffer`, for consistency with `fetch()`.
     */
    const Value_* fetch_copy(Index_ i, Value_* buffer) {
        auto out = fetch(i, buffer);
        if (out != buffer) {
            std::copy(out, out + extracted_length<selection_, Index_>(*this), buffer);
        }
        return buffer;
    }

    /**
     * @param i Index of the desired element on the iteration dimension.
     * @return Vector of length `extracted_length()`, containing the contents of the desired dimension element.
     */
    std::vector<Value_> fetch(Index_ i) {
        std::vector<Value_> buffer(extracted_length<selection_, Index_>(*this));
        fetch_copy(i, buffer.data());
        return buffer;
    }

    /**
     * Whether this class enables sparse access.
     */
    static constexpr bool sparse = false;

    /**
     * Type of selection on the extraction dimension.
     */
    static constexpr DimensionSelectionType selection = selection_;
};

/**
 * @tparam selection_ Type of selection along the extraction dimension.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * @brief Virtual base class for sparse extraction.
 */
template<DimensionSelectionType selection_, typename Value_, typename Index_>
class SparseExtractor : public ConditionalSelectionExtractor<selection_, Index_> {
protected:
    /**
     * @cond
     */
    SparseExtractor() = default;
    /**
     * @endcond
     */

public:
    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * @param i Index of the desired element on the iteration dimension.
     * @param vbuffer Pointer to an array with enough space for at least `extracted_length()` values.
     * Ignored if `Options::sparse_extract_value` was set to `false` during construction of this instance.
     * @param ibuffer Pointer to an array with enough space for at least `extracted_length()` indices.
     * Ignored if `Options::sparse_extract_index` was set to `false` during construction of this instance.
     *
     * @return A `SparseRange` object describing the contents of the desired dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) = 0;

    /**
     * @param i Index of the desired element on the iteration dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `extracted_length()` values.
     * On output, this is filled with the values of the structural non-zeros. 
     *
     * Ignored if `Options::sparse_extract_value` was set to `false` during construction of this instance.
     * Also ignored if set to `NULL`, in which case the values are extracted but not copied to `vbuffer`.
     * @param[out] ibuffer Pointer to an array with enough space for at least `extracted_length()` indices.
     * On output, this is filled with the indices of the structural non-zeros. 
     *
     * Ignored if `Options::sparse_extract_index` was set to `false` during construction of this instance.
     * Also ignored if set to `NULL`, in which case the indices are extracted but not copied to `ibuffer`.
     *
     * @return A `SparseRange` object describing the contents of the desired dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
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
     * @param i Index of the desired element on the iteration dimension.
     * @return A `SparseRangeCopy` object containing the contents of the desired dimension element.
     * Either or both of `value` or `index` is empty if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     */
    SparseRangeCopy<Value_, Index_> fetch(Index_ i) {
        SparseRangeCopy<Value_, Index_> output(extracted_length<selection_, Index_>(*this));
        auto range = fetch_copy(i, output.value.data(), output.index.data());
        output.number = range.number;
        output.value.resize(range.value != NULL ? range.number : 0);
        output.index.resize(range.index != NULL ? range.number : 0);
        return output;
    }

    /**
     * Whether this class enables sparse access.
     */
    static constexpr bool sparse = true;

    /**
     * Type of selection on the extraction dimension.
     */
    static constexpr DimensionSelectionType selection = selection_;
};

/**
 * @tparam selection_ Type of selection on the extraction dimension.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Conditional extractor interface that depends on the format type.
 */
template<DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_> 
using Extractor = typename std::conditional<sparse_, SparseExtractor<selection_, Value_, Index_>, DenseExtractor<selection_, Value_, Index_> >::type;

/**
 * Extractor for dense extraction of full rows.
 */
template<typename Value_, typename Index_>
using FullDenseExtractor = Extractor<DimensionSelectionType::FULL, false, Value_, Index_>;

/**
 * Extractor for dense extraction of a block of each row.
 */
template<typename Value_, typename Index_>
using BlockDenseExtractor = Extractor<DimensionSelectionType::BLOCK, false, Value_, Index_>;

/**
 * Extractor for dense extraction of an indexed subset of each row.
 */
template<typename Value_, typename Index_>
using IndexDenseExtractor = Extractor<DimensionSelectionType::INDEX, false, Value_, Index_>;

/**
 * Extractor for sparse extraction of full rows.
 */
template<typename Value_, typename Index_>
using FullSparseExtractor = Extractor<DimensionSelectionType::FULL, true, Value_, Index_>;

/**
 * Extractor for sparse extraction of a block of each column.
 */
template<typename Value_, typename Index_>
using BlockSparseExtractor = Extractor<DimensionSelectionType::BLOCK, true, Value_, Index_>;

/**
 * Extractor for sparse extraction of an indexed subset of each column.
 */
template<typename Value_, typename Index_>
using IndexSparseExtractor = Extractor<DimensionSelectionType::INDEX, true, Value_, Index_>;

/*********************************************************
 *********************************************************/

/**
 * @tparam selection_ Type of selection along the extraction dimension.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * @brief Virtual base class for oracle-aware dense extraction.
 */
template<DimensionSelectionType selection_, typename Value_, typename Index_>
class DenseOracleAwareExtractor : public ConditionalSelectionExtractor<selection_, Index_> {
protected:
    /**
     * @cond
     */
    DenseOracleAwareExtractor() = default;
    /**
     * @endcond
     */

public:
    /**
     * Number of predictions from the oracle that have already been fulfilled.
     * This can be compared to `total_predictions` to determine whether all predictions have been performed.
     * This will be incremented on any call to `fetch()` or `fetch_copy()`.
     */
    size_t used_predictions = 0;

    /**
     * Total number of predictions from the oracle.
     */
    size_t total_predictions = 0;

public:
    /**
     * `buffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This can be checked by comparing the returned pointer to `buffer`; if they are the same, `buffer` has been filled.
     *
     * Calling this function will automatically advance the `DenseOracleAwareExtractor` instance along the oracle's prediction stream,
     * i.e., `used_predictions` is incremented and the next call to `fetch()` will extract data for the next prediction.
     *
     * @param[out] i Predicted index of the element on the iteration dimension.
     * @param[out] buffer Pointer to an array of length no less than `extracted_length()`.
     *
     * @return Pointer to the contents of the predicted element of the iteration dimension, containing `extracted_length()` entries.
     * The identity of the predicted element is stored in `i`.
     */
    virtual const Value_* fetch(Index_& i, Value_* buffer) = 0;

    /**
     * Calling this function will automatically advance the `DenseOracleAwareExtractor` instance along the oracle's prediction stream,
     * i.e., `used_predictions` is incremented and the next call to `fetch()` will extract data for the next prediction.
     *
     * @param[out] i Predicted index of the element on the iteration dimension.
     * @param[out] buffer Pointer to an array of length specified by `extracted_length()`.
     * This is filled with the contents of the desired dimension element.
     *
     * @return `buffer`, for consistency with `fetch()`.
     * The identity of the predicted element is stored in `i`.
     */
    const Value_* fetch_copy(Index_& i, Value_* buffer) {
        auto out = fetch(i, buffer);
        if (out != buffer) {
            std::copy(out, out + extracted_length<selection_, Index_>(*this), buffer);
        }
        return buffer;
    }

    /**
     * Calling this function will automatically advance the `DenseOracleAwareExtractor` instance along the oracle's prediction stream,
     * i.e., `used_predictions` is incremented and the next call to `fetch()` will extract data for the next prediction.
     *
     * @param[out] i Predicted index of the element on the iteration dimension.
     * @return Vector of length `extracted_length()`, containing the contents of the desired dimension element.
     */
    std::vector<Value_> fetch(Index_& i) {
        std::vector<Value_> buffer(extracted_length<selection_, Index_>(*this));
        fetch_copy(i, buffer.data());
        return buffer;
    }

    /**
     * Whether this class enables sparse access.
     */
    static constexpr bool sparse = false;

    /**
     * Type of selection on the extraction dimension.
     */
    static constexpr DimensionSelectionType selection = selection_;
};

/**
 * @tparam selection_ Type of selection along the extraction dimension.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * @brief Virtual base class for sparse extraction.
 */
template<DimensionSelectionType selection_, typename Value_, typename Index_>
class SparseOracleAwareExtractor : public ConditionalSelectionExtractor<selection_, Index_> {
protected:
    /**
     * @cond
     */
    SparseOracleAwareExtractor() = default;
    /**
     * @endcond
     */

public:
    /**
     * Number of predictions from the oracle that have already been fulfilled.
     * This can be compared to `total_predictions` to determine whether all predictions have been performed.
     * This will be incremented on any call to `fetch()` or `fetch_copy()`.
     */
    size_t used_predictions = 0;

    /**
     * Total number of predictions from the oracle.
     */
    size_t total_predictions = 0;

public:
    /**
     * `vbuffer` may not necessarily be filled upon extraction if a pointer can be returned to the underlying data store.
     * This be checked by comparing the returned `SparseRange::value` pointer to `vbuffer`; if they are the same, `vbuffer` has been filled. 
     * The same applies for `ibuffer` and the returned `SparseRange::index` pointer.
     *
     * Calling this function will automatically advance the `SparseOracleAwareExtractor` instance along the oracle's prediction stream,
     * i.e., `used_predictions` is incremented and the next call to `fetch()` will extract data for the next prediction.
     *
     * @param[out] i Predicted index of the element on the iteration dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `extracted_length()` values.
     * Ignored if `Options::sparse_extract_value` was set to `false` during construction of this instance.
     * @param[out] ibuffer Pointer to an array with enough space for at least `extracted_length()` indices.
     * Ignored if `Options::sparse_extract_index` was set to `false` during construction of this instance.
     *
     * @return A `SparseRange` object describing the contents of the desired dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     * The identity of the predicted element is stored in `i`.
     */
    virtual SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) = 0;

    /**
     * Calling this function will automatically advance the `SparseOracleAwareExtractor` instance along the oracle's prediction stream,
     * i.e., `used_predictions` is incremented and the next call to `fetch()` will extract data for the next prediction.
     *
     * @param[out] i Predicted index of the element on the iteration dimension.
     * @param[out] vbuffer Pointer to an array with enough space for at least `extracted_length()` values.
     * On output, this is filled with the values of the structural non-zeros. 
     *
     * Ignored if `Options::sparse_extract_value` was set to `false` during construction of this instance.
     * Also ignored if set to `NULL`, in which case the values are extracted but not copied to `vbuffer`.
     * @param[out] ibuffer Pointer to an array with enough space for at least `extracted_length()` indices.
     * On output, this is filled with the indices of the structural non-zeros. 
     *
     * Ignored if `Options::sparse_extract_index` was set to `false` during construction of this instance.
     * Also ignored if set to `NULL`, in which case the indices are extracted but not copied to `ibuffer`.
     *
     * @return A `SparseRange` object describing the contents of the desired dimension element.
     * Either or both of `value` or `index` is set to `NULL` if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     * The identity of the predicted element is stored in `i`.
     */
    SparseRange<Value_, Index_> fetch_copy(Index_& i, Value_* vbuffer, Index_* ibuffer) {
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
     * Calling this function will automatically advance the `SparseOracleAwareExtractor` instance along the oracle's prediction stream,
     * i.e., `used_predictions` is incremented and the next call to `fetch()` will extract data for the next prediction.
     *
     * @param[out] i Predicted index of the element on the iteration dimension.
     * @return A `SparseRangeCopy` object containing the contents of the desired dimension element.
     * Either or both of `value` or `index` is empty if extraction of that field is skipped, 
     * based on the setting of `Options::sparse_extract_mode` used to construct this object.
     * The identity of the predicted element is stored in `i`.
     */
    SparseRangeCopy<Value_, Index_> fetch(Index_& i) {
        SparseRangeCopy<Value_, Index_> output(extracted_length<selection_, Index_>(*this));
        auto range = fetch_copy(i, output.value.data(), output.index.data());
        output.number = range.number;
        output.value.resize(range.value != NULL ? range.number : 0);
        output.index.resize(range.index != NULL ? range.number : 0);
        return output;
    }

    /**
     * Whether this class enables sparse access.
     */
    static constexpr bool sparse = true;

    /**
     * Type of selection on the extraction dimension.
     */
    static constexpr DimensionSelectionType selection = selection_;
};

/**
 * @tparam selection_ Type of selection on the extraction dimension.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 *
 * Conditional interface for an oracle-aware extractor that depends on the format type.
 */
template<DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_> 
using OracleAwareExtractor = typename std::conditional<
        sparse_, 
        SparseOracleAwareExtractor<selection_, Value_, Index_>, 
        DenseOracleAwareExtractor<selection_, Value_, Index_> 
    >::type;

/**
 * Oracle-aware extractor for dense extraction of full rows.
 */
template<typename Value_, typename Index_>
using FullDenseOracleAwareExtractor = OracleAwareExtractor<DimensionSelectionType::FULL, false, Value_, Index_>;

/**
 * Oracle-aware extractor for dense extraction of a block of each row.
 */
template<typename Value_, typename Index_>
using BlockDenseOracleAwareExtractor = OracleAwareExtractor<DimensionSelectionType::BLOCK, false, Value_, Index_>;

/**
 * Oracle-aware extractor for dense extraction of an indexed subset of each row.
 */
template<typename Value_, typename Index_>
using IndexDenseOracleAwareExtractor = OracleAwareExtractor<DimensionSelectionType::INDEX, false, Value_, Index_>;

/**
 * Oracle-aware extractor for sparse extraction of full rows.
 */
template<typename Value_, typename Index_>
using FullSparseOracleAwareExtractor = OracleAwareExtractor<DimensionSelectionType::FULL, true, Value_, Index_>;

/**
 * Oracle-aware extractor for sparse extraction of a block of each column.
 */
template<typename Value_, typename Index_>
using BlockSparseOracleAwareExtractor = OracleAwareExtractor<DimensionSelectionType::BLOCK, true, Value_, Index_>;

/**
 * Oracle-aware extractor for sparse extraction of an indexed subset of each column.
 */
template<typename Value_, typename Index_>
using IndexSparseOracleAwareExtractor = OracleAwareExtractor<DimensionSelectionType::INDEX, true, Value_, Index_>;

/*********************************************************
 *********************************************************/

/**
 * @cond
 */
template<DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_> 
struct DummyOracleAwareExtractor : public OracleAwareExtractor<selection_, sparse_, Value_, Index_> {
    DummyOracleAwareExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > ex, std::shared_ptr<Oracle<Index_> > o) : 
        internal(std::move(ex)), oracle(std::move(o)) 
    {
        this->total_predictions = oracle->total();
        if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
            this->full_length = internal->full_length;
        } else if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
            this->block_start = internal->block_start;
            this->block_length = internal->block_length;
        } else {
            this->index_length = internal->index_length;
        }
    }

    const Value_* fetch(Index_& i, Value_* buffer) {
        i = oracle->get(this->used_predictions);
        ++(this->used_predictions);
        return internal->fetch(i, buffer);
    }

    SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
        i = oracle->get(this->used_predictions);
        ++(this->used_predictions);
        return internal->fetch(i, vbuffer, ibuffer);
    }

    const Index_* index_start() const {
        return internal->index_start();
    }

private:
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > internal;
    std::shared_ptr<Oracle<Index_> > oracle;
};
/**
 * @endcond
 */

}

#endif
