#ifndef TATAMI_FRAGMENTED_SPARSE_MATRIX_H
#define TATAMI_FRAGMENTED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "../utils/ElementType.hpp"
#include "../utils/PseudoOracularExtractor.hpp"

#include "primary_extraction.hpp"
#include "secondary_extraction.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

/**
 * @file FragmentedSparseMatrix.hpp
 *
 * @brief Fragmented sparse matrix representation. 
 */

namespace tatami {

/**
 * @cond
 */
namespace FragmentedSparseMatrix_internal {

/********************
 *** Primary full ***
 ********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class PrimaryMyopicFullDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    PrimaryMyopicFullDense(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary) :
        my_values(values), my_indices(indices), my_secondary(secondary) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        const auto& curv = my_values[i];
        const auto& curi = my_indices[i];

        std::fill_n(buffer, my_secondary, static_cast<Value_>(0));
        for (decltype(curv.size()) x = 0, end = curv.size(); x < end; ++x) {
            buffer[curi[x]] = curv[x];
        }
        return buffer;
    }

private:
    const ValueVectorStorage_& my_values;
    const IndexVectorStorage_& my_indices;
    Index_ my_secondary;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class PrimaryMyopicFullSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    PrimaryMyopicFullSparse(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, [[maybe_unused]] Index_ secondary /* for consistency only */, const Options& opt) :
        my_values(values), my_indices(indices), my_needs_value(opt.sparse_extract_value), my_needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* index_buffer) {
        const auto& curv = my_values[i];
        const auto& curi = my_indices[i];

        SparseRange<Value_, Index_> output(curv.size(), NULL, NULL);
        if (my_needs_value) {
            output.value = sparse_utils::extract_primary_vector(curv, static_cast<decltype(curv.size())>(0), curv.size(), vbuffer);
        }
        if (my_needs_index) {
            output.index = sparse_utils::extract_primary_vector(curi, static_cast<decltype(curi.size())>(0), curi.size(), index_buffer);
        }
        return output;
    }

private:
    const ValueVectorStorage_& my_values;
    const IndexVectorStorage_& my_indices;
    bool my_needs_value, my_needs_index;
};

/*********************
 *** Primary block ***
 *********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class PrimaryMyopicBlockDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    PrimaryMyopicBlockDense(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, Index_ block_start, Index_ block_length) :
        my_values(values), my_indices(indices), my_secondary(secondary), my_block_start(block_start), my_block_length(block_length) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        const auto& curi = my_indices[i];
        const auto& curv = my_values[i];

        auto iStart = curi.begin();
        auto iEnd = curi.end();
        sparse_utils::refine_primary_block_limits(iStart, iEnd, my_secondary, my_block_start, my_block_length);
        auto start_pos = (iStart - curi.begin());
        auto end_pos = (iEnd - curi.begin());

        std::fill_n(buffer, my_block_length, static_cast<Value_>(0));
        for (auto x = start_pos; x < end_pos; ++x) {
            buffer[curi[x] - my_block_start] = curv[x];
        }
        return buffer;
    }

private:
    const ValueVectorStorage_& my_values;
    const IndexVectorStorage_& my_indices;
    Index_ my_secondary;
    Index_ my_block_start, my_block_length;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class PrimaryMyopicBlockSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    PrimaryMyopicBlockSparse(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, Index_ block_start, Index_ block_length, const Options& opt) :
        my_values(values),
        my_indices(indices), 
        my_secondary(secondary),
        my_block_start(block_start),
        my_block_length(block_length),
        my_needs_value(opt.sparse_extract_value), 
        my_needs_index(opt.sparse_extract_index)
    {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* index_buffer) {
        const auto& curi = my_indices[i];
        auto iStart = curi.begin();
        auto iEnd = curi.end();
        sparse_utils::refine_primary_block_limits(iStart, iEnd, my_secondary, my_block_start, my_block_length);
        auto offset = iStart - curi.begin();
        auto delta = iEnd - iStart;

        SparseRange<Value_, Index_> output(delta, NULL, NULL);
        if (my_needs_value) {
            output.value = sparse_utils::extract_primary_vector(my_values[i], offset, delta, vbuffer);
        }
        if (my_needs_index) {
            output.index = sparse_utils::extract_primary_vector(curi, offset, delta, index_buffer);
        }
        return output;
    }

private:
    const ValueVectorStorage_& my_values;
    const IndexVectorStorage_& my_indices;
    Index_ my_secondary;
    Index_ my_block_start, my_block_length;
    bool my_needs_value, my_needs_index;
};

/***********************
 *** Primary indexed ***
 ***********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class PrimaryMyopicIndexDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    PrimaryMyopicIndexDense(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, VectorPtr<Index_> indices_ptr) :
        my_values(values), my_indices(indices), my_retriever(*indices_ptr, secondary), my_num_indices(indices_ptr->size()) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        const auto& curi = my_indices[i];
        const auto& curv = my_values[i];
        std::fill_n(buffer, my_num_indices, static_cast<Value_>(0));
        my_retriever.populate(
            curi.begin(),
            curi.end(),
            [&](auto s, auto offset) -> void {
                buffer[s] = curv[offset];
            }
        );
        return buffer;
    }

private:
    const ValueVectorStorage_& my_values;
    const IndexVectorStorage_& my_indices;
    sparse_utils::RetrievePrimarySubsetDense<Index_> my_retriever;
    std::size_t my_num_indices;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class PrimaryMyopicIndexSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    PrimaryMyopicIndexSparse(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, VectorPtr<Index_> indices_ptr, const Options& opt) :
        my_values(values), my_indices(indices), my_retriever(*indices_ptr, secondary), my_needs_value(opt.sparse_extract_value), my_needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* index_buffer) {
        const auto& curi = my_indices[i];
        const auto& curv = my_values[i];
        Index_ count = 0;
        auto vcopy = vbuffer;
        auto icopy = index_buffer;

        my_retriever.populate(
            curi.begin(),
            curi.end(),
            [&](auto offset, auto ix) -> void {
                ++count;
                if (my_needs_value) {
                    *vcopy = curv[offset];
                    ++vcopy;
                }
                if (my_needs_index) {
                    *icopy = ix;
                    ++icopy;
                }
            }
        );

        return SparseRange<Value_, Index_>(count, my_needs_value ? vbuffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueVectorStorage_& my_values;
    const IndexVectorStorage_& my_indices;
    sparse_utils::RetrievePrimarySubsetSparse<Index_> my_retriever;
    bool my_needs_value, my_needs_index;
};

/**********************
 *** Secondary full ***
 **********************/

template<typename Index_, class IndexVectorStorage_>
class ServeIndices {
public:
    ServeIndices(const IndexVectorStorage_& indices) : my_indices(indices) {}

private:
    const IndexVectorStorage_& my_indices;

public:
    typedef decltype(std::declval<IndexVectorStorage_>()[0].size()) pointer_type;

    pointer_type start_offset(Index_) const {
        return 0;
    }

    pointer_type end_offset(Index_ primary) const {
        return my_indices[primary].size();
    }

    auto raw(Index_ primary) const {
        return my_indices[primary].begin();
    }
};

template<typename Index_, class IndexVectorStorage_>
auto make_ServeIndices(const IndexVectorStorage_& i) {
    return ServeIndices<Index_, IndexVectorStorage_>(i);
}

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_> 
class SecondaryMyopicFullDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    SecondaryMyopicFullDense(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices), secondary, indices.size()) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_cache.size(), static_cast<Value_>(0));
        my_cache.search(
            i, 
            [&](Index_ primary, Index_ index_primary, auto ptr) -> void {
                buffer[index_primary] = my_values[primary][ptr];
            }
        );
        return buffer;
    }

private:
    const ValueVectorStorage_& my_values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > my_cache;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class SecondaryMyopicFullSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    SecondaryMyopicFullSparse(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, const Options& opt) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices), secondary, indices.size()), my_needs_value(opt.sparse_extract_value), my_needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        Index_ count = 0;
        my_cache.search(
            i,
            [&](Index_ primary, Index_, auto ptr) -> void {
                if (my_needs_value) {
                    value_buffer[count] = my_values[primary][ptr];
                }
                if (my_needs_index) {
                    index_buffer[count] = primary;
                }
                ++count;
            }
        );
        return SparseRange<Value_, Index_>(count, my_needs_value ? value_buffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueVectorStorage_& my_values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > my_cache;
    bool my_needs_value, my_needs_index;
};

/***********************
 *** Secondary block ***
 ***********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class SecondaryMyopicBlockDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    SecondaryMyopicBlockDense(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, Index_ block_start, Index_ block_length) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices), secondary, block_start, block_length) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_cache.size(), static_cast<Value_>(0));
        my_cache.search(
            i,
            [&](Index_ primary, Index_ index_primary, auto ptr) -> void {
                buffer[index_primary] = my_values[primary][ptr];
            }
        );
        return buffer;
    }

private:
    const ValueVectorStorage_& my_values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > my_cache;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class SecondaryMyopicBlockSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    SecondaryMyopicBlockSparse(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, Index_ block_start, Index_ block_length, const Options& opt) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices), secondary, block_start, block_length), my_needs_value(opt.sparse_extract_value), my_needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        Index_ count = 0;
        my_cache.search(
            i, 
            [&](Index_ primary, Index_, auto ptr) -> void {
                if (my_needs_value) {
                    value_buffer[count] = my_values[primary][ptr];
                }
                if (my_needs_index) {
                    index_buffer[count] = primary;
                }
                ++count;
            }
        );
        return SparseRange<Value_, Index_>(count, my_needs_value ? value_buffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueVectorStorage_& my_values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > my_cache;
    bool my_needs_value, my_needs_index;
};

/***********************
 *** Secondary index ***
 ***********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class SecondaryMyopicIndexDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    SecondaryMyopicIndexDense(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, VectorPtr<Index_> indices_ptr) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices), secondary, std::move(indices_ptr)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill_n(buffer, my_cache.size(), static_cast<Value_>(0));
        my_cache.search(
            i,
            [&](Index_ primary, Index_ index_primary, auto ptr) -> void {
                buffer[index_primary] = my_values[primary][ptr];
            }
        );
        return buffer;
    }

private:
    const ValueVectorStorage_& my_values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > my_cache;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class SecondaryMyopicIndexSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    SecondaryMyopicIndexSparse(const ValueVectorStorage_& values, const IndexVectorStorage_& indices, Index_ secondary, VectorPtr<Index_> indices_ptr, const Options& opt) :
        my_values(values), my_cache(make_ServeIndices<Index_>(indices), secondary, std::move(indices_ptr)), my_needs_value(opt.sparse_extract_value), my_needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* index_buffer) {
        Index_ count = 0;
        my_cache.search(
            i,
            [&](Index_ primary, Index_, auto ptr) -> void {
                if (my_needs_value) {
                    vbuffer[count] = my_values[primary][ptr];
                }
                if (my_needs_index) {
                    index_buffer[count] = primary;
                }
                ++count;
            }
        );
        return SparseRange<Value_, Index_>(count, my_needs_value ? vbuffer : NULL, my_needs_index ? index_buffer : NULL);
    }

private:
    const ValueVectorStorage_& my_values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > my_cache;
    bool my_needs_value, my_needs_index;
};

}
/**
 * @endcond
 */

/**
 * @brief Options for `FragmentedSparseMatrix()`.
 */
struct FragmentedSparseMatrixOptions {
    /**
     * Should the input vectors be checked for validity in the `FragmentedSparseMatrix` constructor?
     * If true, the constructor will check that:
     *
     * - `values` and `indices` have the same length that is equal to the number of rows (for `row_sparse = true`) or columns (otherwise).
     * - corresponding elements of `values` and `indices` have the same length.
     * - each element of `indices` is ordered and contains non-negative values less than `ncol` (for `row_sparse = true`) or `nrow` (otherwise).
     *
     * This can be disabled for faster construction if the caller knows that the input vectors are valid.
     */
    bool check = true;
};

/**
 * @brief Fragmented sparse matrix representation.
 *
 * In a fragmented sparse matrix, each element of the primary dimension has its own vector of indices and data values.
 * This differs from a compressed sparse matrix (see `CompressedSparseMatrix`) where the index/value vectors are concatenated across all elements.
 * For row sparse matrices, the rows are the primary dimension, while for column sparse matrices, the columns are the primary dimension.
 * This representation is equivalent to SciPy's list-of-lists sparse matrix (Python), or SparseArray's SVT_SparseMatrix class (R/Bioconductor).
 *
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam ValueVectorStorage_ Vector class used to store the matrix value vectors.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * Each inner vector should also have methods for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Value_*`, it will also be used.
 * The inner vector does not necessarily have to contain `Value_`, as long as the type is convertible to `Value_`.
 * @tparam IndexVectorStorage_ Vector class used to store the row/column indices internally.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * Each inner vector should also have methods for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Index*`, it will also be used.
 * The inner vector does not necessarily have to contain `Index_`, as long as the type is convertible to `Index_`.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
class FragmentedSparseMatrix : public Matrix<Value_, Index_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of vectors of non-zero elements.
     * @param indices Vector of vectors of row indices (if `row_sparse = false`) or column indices (if `row_sparse = true`) for the non-zero elements.
     * @param row_sparse Whether this is a row sparse representation.
     * If `false`, a column sparse representation is assumed instead.
     * @param options Further options.
     */
    FragmentedSparseMatrix(Index_ nrow, Index_ ncol, ValueVectorStorage_ values, IndexVectorStorage_ indices, bool row_sparse, const FragmentedSparseMatrixOptions& options) : 
        my_nrow(nrow), my_ncol(ncol), my_values(std::move(values)), my_indices(std::move(indices)), my_row_sparse(row_sparse)
    {
        if (options.check) {
            if (my_values.size() != my_indices.size()) {
                throw std::runtime_error("'values' and 'indices' should be of the same length");
            }

            if (my_row_sparse) {
                if (!safe_non_negative_equal(my_indices.size(), my_nrow)) {
                    throw std::runtime_error("length of 'indices' should be equal to number of rows'");
                }
            } else {
                if (!safe_non_negative_equal(my_indices.size(), my_ncol)) {
                    throw std::runtime_error("length of 'indices' should be equal to number of columns");
                }
            }

            ElementType<ElementType<IndexVectorStorage_> > max_index = (my_row_sparse ? my_ncol : my_nrow);
            for (decltype(my_indices.size()) i = 0, end = my_indices.size(); i < end; ++i) {
                const auto& curv = my_values[i];
                const auto& curi = my_indices[i];
                if (!safe_non_negative_equal(curv.size(), curi.size())) {
                    throw std::runtime_error("corresponding elements of 'values' and 'indices' should have the same length");
                }

                for (auto x : curi) {
                    if (x < 0 || x >= max_index) {
                        throw std::runtime_error("'indices' should contain non-negative integers less than the number of " + (my_row_sparse ? std::string("columns") : std::string("rows")));
                    }
                }

                for (decltype(curi.size()) j = 1, jend = curi.size(); j < jend; ++j) {
                    if (curi[j] <= curi[j - 1]) {
                        throw std::runtime_error("my_indices should be strictly increasing within each element of 'indices'");
                    }
                }

                // Check that iterator subtraction is safe.
                // Various functions in 'sparse_utils' will subtract iterators when converting the lower_bound return value to an index.
                // We cast to Index_ as it gives us a chance to skip the check at compile time, given that the size is no greater than the dimension extents.
                sanisizer::can_ptrdiff<decltype(curi.begin())>(static_cast<Index_>(curi.size()));
            }
        }
    }

    /**
     * @cond
     */
    // Back-compatibility only.
    FragmentedSparseMatrix(Index_ nrow, Index_ ncol, ValueVectorStorage_ values, IndexVectorStorage_ indices, bool row_sparse, bool check = true) :
        FragmentedSparseMatrix(
            nrow,
            ncol,
            std::move(values),
            std::move(indices),
            row_sparse,
            [&]{
                FragmentedSparseMatrixOptions fopt;
                fopt.check = check;
                return fopt;
            }()
        )
    {}
    /**
     * @endcond
     */

private:
    Index_ my_nrow, my_ncol;
    ValueVectorStorage_ my_values;
    IndexVectorStorage_ my_indices;
    bool my_row_sparse;

public:
    Index_ nrow() const { return my_nrow; }

    Index_ ncol() const { return my_ncol; }

    bool is_sparse() const { return true; }

    double is_sparse_proportion() const { return 1; }

    bool prefer_rows() const { return my_row_sparse; }

    double prefer_rows_proportion() const { return static_cast<double>(my_row_sparse); }

    bool uses_oracle(bool) const { return false; }

    using Matrix<Value_, Index_>::dense;

    using Matrix<Value_, Index_>::sparse;

private:
    Index_ secondary() const {
        if (my_row_sparse) {
            return my_ncol;
        } else {
            return my_nrow;
        }
    }

    /*****************************
     ******* Dense myopic ********
     *****************************/
private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options&) const {
        if (my_row_sparse == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicFullDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary()
            );
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicFullDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary()
            ); 
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_end, const Options&) const {
        if (my_row_sparse == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicBlockDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), block_start, block_end
            );
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicBlockDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), block_start, block_end
            );
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> subset_ptr, const Options&) const {
        if (my_row_sparse == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicIndexDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), std::move(subset_ptr)
            );
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicIndexDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), std::move(subset_ptr)
            );
        }
    }

    /******************************
     ******* Sparse myopic ********
     ******************************/
private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        if (my_row_sparse == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicFullSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), opt
            );
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicFullSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), opt
            ); 
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_end, const Options& opt) const {
        if (my_row_sparse == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicBlockSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), block_start, block_end, opt
            );
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicBlockSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), block_start, block_end, opt
            );
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> subset_ptr, const Options& opt) const {
        if (my_row_sparse == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicIndexSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), std::move(subset_ptr), opt
            );
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicIndexSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(
                my_values, my_indices, secondary(), std::move(subset_ptr), opt
            );
        }
    }

    /*******************************
     ******* Dense oracular ********
     *******************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, block_start, block_end, opt));
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> subset_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, std::move(subset_ptr), opt));
    }

    /********************************
     ******* Sparse oracular ********
     ********************************/
public:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, Index_ block_start, Index_ block_end, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, block_start, block_end, opt));
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> subset_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, std::move(subset_ptr), opt));
    }
};

/**
 * @brief Fragmented sparse column matrix.
 *
 * See `tatami::FragmentedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_ = std::vector<std::vector<Value_> >, class IndexVectorStorage_ = std::vector<std::vector<Index_> > >
class FragmentedSparseColumnMatrix final : public FragmentedSparseMatrix<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of vectors of non-zero elements.
     * @param indices Vector of vectors of row indices for the non-zero elements.
     * @param check Should the input vectors be checked for validity?
     */
    FragmentedSparseColumnMatrix(Index_ nrow, Index_ ncol, ValueVectorStorage_ values, IndexVectorStorage_ indices, bool check = true) : 
        FragmentedSparseMatrix<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_>(nrow, ncol, std::move(values), std::move(indices), false, check) {}
};

/**
 * @brief Fragmented sparse row matrix.
 *
 * See `tatami::FragmentedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_ = std::vector<std::vector<Value_> >, class IndexVectorStorage_ = std::vector<std::vector<Index_> > >
class FragmentedSparseRowMatrix final : public FragmentedSparseMatrix<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> {
public:
    /**
     * @param nrow Number of rows.
     * @param ncol Number of columns.
     * @param values Vector of vectors of non-zero elements.
     * @param indices Vector of vectors of column indices for the non-zero elements.
     * @param check Should the input vectors be checked for validity?
     */
    FragmentedSparseRowMatrix(Index_ nrow, Index_ ncol, ValueVectorStorage_ values, IndexVectorStorage_ indices, bool check = true) : 
        FragmentedSparseMatrix<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_>(nrow, ncol, std::move(values), std::move(indices), true, check) {}
};


}

#endif
