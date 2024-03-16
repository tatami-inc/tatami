#ifndef TATAMI_FRAGMENTED_SPARSE_MATRIX_H
#define TATAMI_FRAGMENTED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "primary_extraction.hpp"
#include "secondary_extraction.hpp"
#include "../utils/ElementType.hpp"
#include "../utils/PseudoOracularExtractor.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>

/**
 * @file FragmentedSparseMatrix.hpp
 *
 * @brief Fragmented sparse matrix representation. 
 *
 * `typedef`s are provided for the usual row and column formats. 
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
struct PrimaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicFullDense(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec) :
        values(vstore), indices(istore), secondary(sec) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        const auto& curv = values[i];
        const auto& curi = indices[i];

        std::fill(buffer, buffer + secondary, static_cast<Value_>(0));
        for (size_t x = 0, end = curv.size(); x < end; ++x) {
            buffer[curi[x]] = curv[x];
        }
        return buffer;
    }

private:
    const ValueVectorStorage_& values;
    const IndexVectorStorage_& indices;
    Index_ secondary;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct PrimaryMyopicFullSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicFullSparse(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, const Options& opt) :
        values(vstore), indices(istore), secondary(sec), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        const auto& curv = values[i];
        const auto& curi = indices[i];

        SparseRange<Value_, Index_> output(curv.size(), NULL, NULL);
        if (needs_value) {
            output.value = sparse_utils::extract_primary_vector(curv, static_cast<size_t>(0), curv.size(), vbuffer);
        }
        if (needs_index) {
            output.index = sparse_utils::extract_primary_vector(curi, static_cast<size_t>(0), curi.size(), ibuffer);
        }
        return output;
    }

private:
    const ValueVectorStorage_& values;
    const IndexVectorStorage_& indices;
    Index_ secondary;
    bool needs_value, needs_index;
};

/*********************
 *** Primary block ***
 *********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct PrimaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicBlockDense(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, Index_ bs, Index_ bl) :
        values(vstore), indices(istore), secondary(sec), block_start(bs), block_length(bl) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        const auto& curi = indices[i];
        auto iStart = curi.begin();
        auto iEnd = curi.end();
        sparse_utils::refine_primary_block_limits(iStart, iEnd, secondary, block_start, block_length);

        std::fill(buffer, buffer + block_length, static_cast<Value_>(0));
        auto vIt = values[i].begin() + (iStart - curi.begin());
        for (; iStart != iEnd; ++iStart, ++vIt) {
            buffer[*iStart - block_start] = *vIt;
        }
        return buffer;
    }

private:
    const ValueVectorStorage_& values;
    const IndexVectorStorage_& indices;
    Index_ secondary;
    Index_ block_start, block_length;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct PrimaryMyopicBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicBlockSparse(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, Index_ bs, Index_ bl, const Options& opt) :
        values(vstore), indices(istore), secondary(sec), block_start(bs), block_length(bl), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        const auto& curi = indices[i];
        auto iStart = curi.begin();
        auto iEnd = curi.end();
        sparse_utils::refine_primary_block_limits(iStart, iEnd, secondary, block_start, block_length);
        size_t offset = iStart - curi.begin();
        size_t delta = iEnd - iStart;

        SparseRange<Value_, Index_> output(delta, NULL, NULL);
        if (needs_value) {
            output.value = sparse_utils::extract_primary_vector(values[i], offset, delta, vbuffer);
        }
        if (needs_index) {
            output.index = sparse_utils::extract_primary_vector(curi, offset, delta, ibuffer);
        }
        return output;
    }

private:
    const ValueVectorStorage_& values;
    const IndexVectorStorage_& indices;
    Index_ secondary;
    Index_ block_start, block_length;
    bool needs_value, needs_index;
};

/***********************
 *** Primary indexed ***
 ***********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct PrimaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicIndexDense(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, VectorPtr<Index_> sub_ptr) :
        values(vstore), indices(istore), subset_ptr(std::move(sub_ptr)) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        const auto& curi = indices[i];
        const auto& curv = values[i];
        const auto& subset = *subset_ptr;
        std::fill(buffer, buffer + subset.size(), static_cast<Value_>(0));
        sparse_utils::retrieve_primary_subset(
            curi.begin(),
            curi.end(),
            subset,
            [&](size_t s, size_t offset, Index_) {
                buffer[s] = curv[offset];
            }
        );
        return buffer;
    }

private:
    const ValueVectorStorage_& values;
    const IndexVectorStorage_& indices;
    VectorPtr<Index_> subset_ptr;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct PrimaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicIndexSparse(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, VectorPtr<Index_> sub_ptr, const Options& opt) :
        values(vstore), indices(istore), subset_ptr(std::move(sub_ptr)), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        const auto& curi = indices[i];
        const auto& curv = values[i];
        Index_ count = 0;
        auto vcopy = vbuffer;
        auto icopy = ibuffer;

        sparse_utils::retrieve_primary_subset(
            curi.begin(),
            curi.end(),
            *subset_ptr,
            [&](size_t, size_t offset, Index_ ix) {
                ++count;
                if (needs_value) {
                    *vcopy = curv[offset];
                    ++vcopy;
                }
                if (needs_index) {
                    *icopy = ix;
                    ++icopy;
                }
            }
        );

        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueVectorStorage_& values;
    const IndexVectorStorage_& indices;
    VectorPtr<Index_> subset_ptr;
    bool needs_value, needs_index;
};

/**********************
 *** Secondary full ***
 **********************/

template<typename Index_, class IndexVectorStorage_>
struct ServeIndices {
    ServeIndices(const IndexVectorStorage_& i) : indices(i) {}
    const IndexVectorStorage_& indices;

public:
    typedef size_t pointer_type;

    pointer_type start_offset(Index_) const {
        return 0;
    }

    pointer_type end_offset(Index_ primary) const {
        return indices[primary].size();
    }

    auto raw(Index_ primary) const {
        return indices[primary].begin();
    }
};

template<typename Index_, class IndexVectorStorage_>
auto make_ServeIndices(const IndexVectorStorage_& i) {
    return ServeIndices<Index_, IndexVectorStorage_>(i);
}

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_> 
struct SecondaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicFullDense(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec) :
        values(vstore), cache(make_ServeIndices<Index_>(istore), sec, istore.size()) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + cache.size(), static_cast<Value_>(0));
        cache.search(i, [&](Index_ primary, Index_ index_primary, size_t ptr) {
            buffer[index_primary] = values[primary][ptr];
        });
        return buffer;
    }

private:
    const ValueVectorStorage_& values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > cache;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct SecondaryMyopicFullSparse : public MyopicSparseExtractor<Value_, Index_> {
    SecondaryMyopicFullSparse(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, const Options& opt) :
        values(vstore), cache(make_ServeIndices<Index_>(istore), sec, istore.size()), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        cache.search(i, [&](Index_ primary, Index_, size_t ptr) {
            if (needs_value) {
                vbuffer[count] = values[primary][ptr];
            }
            if (needs_index) {
                ibuffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueVectorStorage_& values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > cache;
    bool needs_value, needs_index;
};

/***********************
 *** Secondary block ***
 ***********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct SecondaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicBlockDense(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, Index_ bs, Index_ bl) :
        values(vstore), cache(make_ServeIndices<Index_>(istore), sec, bs, bl) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + cache.size(), static_cast<Value_>(0));
        cache.search(i, [&](Index_ primary, Index_ index_primary, size_t ptr) {
            buffer[index_primary] = values[primary][ptr];
        });
        return buffer;
    }

private:
    const ValueVectorStorage_& values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > cache;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct SecondaryMyopicBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
    SecondaryMyopicBlockSparse(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, Index_ bs, Index_ bl, const Options& opt) :
        values(vstore), cache(make_ServeIndices<Index_>(istore), sec, bs, bl), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        cache.search(i, [&](Index_ primary, Index_, size_t ptr) {
            if (needs_value) {
                vbuffer[count] = values[primary][ptr];
            }
            if (needs_index) {
                ibuffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueVectorStorage_& values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > cache;
    bool needs_value, needs_index;
};

/***********************
 *** Secondary index ***
 ***********************/

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct SecondaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicIndexDense(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, VectorPtr<Index_> sub_ptr) :
        values(vstore), cache(make_ServeIndices<Index_>(istore), sec, std::move(sub_ptr)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + cache.size(), static_cast<Value_>(0));
        cache.search(i, [&](Index_ primary, Index_ index_primary, size_t ptr) {
            buffer[index_primary] = values[primary][ptr];
        });
        return buffer;
    }

private:
    const ValueVectorStorage_& values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > cache;
};

template<typename Value_, typename Index_, class ValueVectorStorage_, class IndexVectorStorage_>
struct SecondaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
    SecondaryMyopicIndexSparse(const ValueVectorStorage_& vstore, const IndexVectorStorage_& istore, Index_ sec, VectorPtr<Index_> sub_ptr, const Options& opt) :
        values(vstore), cache(make_ServeIndices<Index_>(istore), sec, std::move(sub_ptr)), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        cache.search(i, [&](Index_ primary, Index_, size_t ptr) {
            if (needs_value) {
                vbuffer[count] = values[primary][ptr];
            }
            if (needs_index) {
                ibuffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueVectorStorage_& values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexVectorStorage_> > cache;
    bool needs_value, needs_index;
};

}
/**
 * @endcond
 */


/**
 * @brief Fragmented sparse matrix representation.
 *
 * In a fragmented sparse matrix, each element of the primary dimension has its own vector of indices and data values.
 * This differs from a compressed sparse matrix (see `CompressedSparseMatrix`) where the index/value vectors are concatenated across all elements.
 * For row sparse matrices, the rows are the primary dimension, while for column sparse matrices, the columns are the primary dimension.
 *
 * @tparam row_ Whether this is a row sparse representation.
 * If `false`, a column sparse representation is assumed instead.
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
template<
    bool row_, 
    typename Value_, 
    typename Index_ = int, 
    class ValueVectorStorage_ = std::vector<std::vector<Value_> >,
    class IndexVectorStorage_ = std::vector<std::vector<Index_> >
>
class FragmentedSparseMatrix : public Matrix<Value_, Index_> {
public:
    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param vals Vector of vectors of non-zero elements.
     * @param idx Vector of vectors of row indices (if `ROW=false`) or column indices (if `ROW=true`) for the non-zero elements.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `vals` and `idx` have the same length that is equal to the number of rows (for `row_ = true`) or columns (otherwise);
     * that corresponding elements of `vals` and `idx` also have the same length;
     * and that each `idx` is ordered and contains non-negative values less than `nc` (for `row_ = true`) or `nr` (for `row_ = false`).
     */
    FragmentedSparseMatrix(Index_ nr, Index_ nc, ValueVectorStorage_ vals, IndexVectorStorage_ idx, bool check=true) : 
        nrows(nr), ncols(nc), values(std::move(vals)), indices(std::move(idx)) 
    {
        if (check) {
            if (values.size() != indices.size()) {
                throw std::runtime_error("'values' and 'indices' should be of the same length");
            }

            if (row_) {
                if (indices.size() != static_cast<size_t>(nrows)) {
                    throw std::runtime_error("length of 'indices' should be equal to number of rows'");
                }
            } else {
                if (indices.size() != static_cast<size_t>(ncols)) {
                    throw std::runtime_error("length of 'indices' should be equal to number of columns");
                }
            }

            ElementType<ElementType<IndexVectorStorage_> > max_index = (row_ ? ncols : nrows);
            for (size_t i = 0, end = indices.size(); i < end; ++i) {
                const auto& curv = values[i];
                const auto& curi = indices[i];
                if (curv.size() != curi.size()) {
                    throw std::runtime_error("corresponding elements of 'values' and 'indices' should have the same length");
                }

                for (auto x : curi) {
                    if (x < 0 || x >= max_index) {
                        if constexpr(row_) {
                            throw std::runtime_error("'indices' should contain non-negative integers less than the number of rows");
                        } else {
                            throw std::runtime_error("'indices' should contain non-negative integers less than the number of columns");
                        }
                    }
                }

                for (size_t j = 1, jend = curi.size(); j < jend; ++j) {
                    if (curi[j] <= curi[j - 1]) {
                        throw std::runtime_error("indices should be strictly increasing within each element of 'indices'");
                    }
                }
            }
        }
    }

private:
    Index_ nrows, ncols;
    ValueVectorStorage_ values;
    IndexVectorStorage_ indices;

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool sparse() const { return true; }

    double sparse_proportion() const { return 1; }

    /**
     * @return `true` if `row_ = true` (for `FragmentedSparseRowMatrix` objects), otherwise returns `false` (for `FragmentedSparseColumnMatrix` objects).
     */
    bool prefer_rows() const { return row_; }

    double prefer_rows_proportion() const { return static_cast<double>(row_); }

    bool uses_oracle(bool) const { return false; }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    Index_ secondary() const {
        if constexpr(row_) {
            return ncols;
        } else {
            return nrows;
        }
    }

    /*****************************
     ******* Dense myopic ********
     *****************************/
private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options&) const {
        if (row_ == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicFullDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary());
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicFullDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary()); 
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_end, const Options&) const {
        if (row_ == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicBlockDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), block_start, block_end);
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicBlockDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), block_start, block_end);
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> subset_ptr, const Options&) const {
        if (row_ == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicIndexDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, std::move(subset_ptr));
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicIndexDense<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), std::move(subset_ptr));
        }
    }

    /******************************
     ******* Sparse myopic ********
     ******************************/
private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        if (row_ == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicFullSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), opt);
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicFullSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), opt); 
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_end, const Options& opt) const {
        if (row_ == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicBlockSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), block_start, block_end, opt);
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicBlockSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), block_start, block_end, opt);
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> subset_ptr, const Options& opt) const {
        if (row_ == row) {
            return std::make_unique<FragmentedSparseMatrix_internal::PrimaryMyopicIndexSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, std::move(subset_ptr), opt);
        } else {
            return std::make_unique<FragmentedSparseMatrix_internal::SecondaryMyopicIndexSparse<Value_, Index_, ValueVectorStorage_, IndexVectorStorage_> >(values, indices, secondary(), std::move(subset_ptr), opt);
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
 * Fragmented sparse column matrix.
 * See `tatami::FragmentedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_ = std::vector<std::vector<Value_> >, class IndexVectorStorage_ = std::vector<std::vector<Index_> > >
using FragmentedSparseColumnMatrix = FragmentedSparseMatrix<false, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_>;

/**
 * Fragmented sparse row matrix.
 * See `tatami::FragmentedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_ = std::vector<std::vector<Value_> >, class IndexVectorStorage_ = std::vector<std::vector<Index_> > >
using FragmentedSparseRowMatrix = FragmentedSparseMatrix<true, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_>;

}

#endif
