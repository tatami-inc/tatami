#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "primary_extraction.hpp"
#include "secondary_extraction.hpp"

#include "../utils/ElementType.hpp"
#include "../utils/has_data.hpp"
#include "../utils/PseudoOracularExtractor.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>

/**
 * @file CompressedSparseMatrix.hpp
 *
 * @brief Compressed sparse matrix representation. 
 *
 * `typedef`s are provided for the usual row and column formats. 
 */

namespace tatami {

/**
 * @cond
 */
namespace CompressedSparseMatrix_internal {

/********************
 *** Primary full ***
 ********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicFullDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto offset = indptr[i];
        auto vIt = values.begin() + offset;
        auto iIt = indices.begin() + offset;
        size_t delta = indptr[i+1] - indptr[i];

        std::fill(buffer, buffer + secondary, static_cast<Value_>(0));
        for (size_t x = 0; x < delta; ++x, ++vIt, ++iIt) {
            buffer[*iIt] = *vIt;
        }
        return buffer;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicFullSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicFullSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, const Options& opt) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto offset = indptr[i];
        auto delta = indptr[i+1] - indptr[i];

        SparseRange<Value_, Index_> output(delta, NULL, NULL);
        if (needs_value) {
            output.value = sparse_utils::extract_primary_vector(values, offset, delta, vbuffer);
        }
        if (needs_index) {
            output.index = sparse_utils::extract_primary_vector(indices, offset, delta, ibuffer);
        }
        return output;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
    bool needs_value, needs_index;
};

/*********************
 *** Primary block ***
 *********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicBlockDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, Index_ bs, Index_ bl) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec), block_start(bs), block_length(bl) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto iStart = indices.begin() + indptr[i];
        auto iEnd = indices.begin() + indptr[i + 1];
        sparse_utils::refine_primary_block_limits(iStart, iEnd, secondary, block_start, block_length);

        std::fill(buffer, buffer + block_length, static_cast<Value_>(0));
        auto vIt = values.begin() + (iStart - indices.begin());
        for (; iStart != iEnd; ++iStart, ++vIt) {
            buffer[*iStart - block_start] = *vIt;
        }
        return buffer;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
    Index_ block_start, block_length;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicBlockSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, Index_ bs, Index_ bl, const Options& opt) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec), block_start(bs), block_length(bl), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto iStart = indices.begin() + indptr[i];
        auto iEnd = indices.begin() + indptr[i + 1];
        sparse_utils::refine_primary_block_limits(iStart, iEnd, secondary, block_start, block_length);
        size_t offset = iStart - indices.begin();
        size_t delta = iEnd - iStart;

        SparseRange<Value_, Index_> output(delta, NULL, NULL);
        if (needs_value) {
            output.value = sparse_utils::extract_primary_vector(values, offset, delta, vbuffer);
        }
        if (needs_index) {
            output.index = sparse_utils::extract_primary_vector(indices, offset, delta, ibuffer);
        }
        return output;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
    Index_ block_start, block_length;
    bool needs_value, needs_index;
};

/***********************
 *** Primary indexed ***
 ***********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicIndexDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, const VectorPtr<Index_>& sub_ptr) :
        values(vstore), indices(istore), indptr(pstore), retriever(*sub_ptr), extent(sub_ptr->size()) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + extent, static_cast<Value_>(0));
        auto vIt = values.begin() + indptr[i];
        retriever.populate(
            indices.begin() + indptr[i], 
            indices.begin() + indptr[i+1],
            [&](size_t s, size_t offset) {
                buffer[s] = *(vIt + offset);
            }
        );
        return buffer;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    sparse_utils::RetrievePrimarySubsetDense<Index_> retriever;
    size_t extent;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicIndexSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, const VectorPtr<Index_>& sub_ptr, const Options& opt) :
        values(vstore), indices(istore), indptr(pstore), retriever(*sub_ptr), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        auto vcopy = vbuffer;
        auto icopy = ibuffer;

        auto vIt = values.begin() + indptr[i];
        retriever.populate(
            indices.begin() + indptr[i], 
            indices.begin() + indptr[i+1],
            [&](size_t offset, Index_ ix) {
                ++count;
                if (needs_value) {
                    *vcopy = *(vIt + offset);
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
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    sparse_utils::RetrievePrimarySubsetSparse retriever;
    bool needs_value, needs_index;
};

/**********************
 *** Secondary full ***
 **********************/

template<typename Index_, class IndexStorage_, class PointerStorage_>
struct ServeIndices {
    ServeIndices(const IndexStorage_& i, const PointerStorage_& p) : indices(i), indptr(p) {}
    const IndexStorage_& indices;
    const PointerStorage_& indptr;

public:
    typedef ElementType<PointerStorage_> pointer_type;

    pointer_type start_offset(Index_ primary) const {
        return indptr[primary];
    }

    pointer_type end_offset(Index_ primary) const {
        return indptr[primary + 1];
    }

    auto raw(Index_) const {
        return indices.begin();
    }
};

template<typename Index_, class IndexStorage_, class PointerStorage_>
auto make_ServeIndices(const IndexStorage_& i, const PointerStorage_& p) {
    return ServeIndices<Index_, IndexStorage_, PointerStorage_>(i, p);
}

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicFullDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec) :
        values(vstore), cache(make_ServeIndices<Index_>(istore, pstore), sec, pstore.size() - 1) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + cache.size(), static_cast<Value_>(0));
        cache.search(i, [&](Index_, Index_ index_primary, ElementType<PointerStorage_> ptr) {
            buffer[index_primary] = values[ptr];
        });
        return buffer;
    }

private:
    const ValueStorage_& values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > cache;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicFullSparse : public MyopicSparseExtractor<Value_, Index_> {
    SecondaryMyopicFullSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, const Options& opt) :
        values(vstore), cache(make_ServeIndices<Index_>(istore, pstore), sec, pstore.size() - 1), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        cache.search(i, [&](Index_ primary, Index_, ElementType<PointerStorage_> ptr) {
            if (needs_value) {
                vbuffer[count] = values[ptr];
            }
            if (needs_index) {
                ibuffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueStorage_& values;
    sparse_utils::FullSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > cache;
    bool needs_value, needs_index;
};

/***********************
 *** Secondary block ***
 ***********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicBlockDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, Index_ bs, Index_ bl) :
        values(vstore), cache(make_ServeIndices<Index_>(istore, pstore), sec, bs, bl) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + cache.size(), static_cast<Value_>(0));
        cache.search(i, [&](Index_, Index_ index_primary, ElementType<PointerStorage_> ptr) {
            buffer[index_primary] = values[ptr];
        });
        return buffer;
    }

private:
    const ValueStorage_& values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > cache;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicBlockSparse : public MyopicSparseExtractor<Value_, Index_> {
    SecondaryMyopicBlockSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, Index_ bs, Index_ bl, const Options& opt) :
        values(vstore), cache(make_ServeIndices<Index_>(istore, pstore), sec, bs, bl), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        cache.search(i, [&](Index_ primary, Index_, ElementType<PointerStorage_> ptr) {
            if (needs_value) {
                vbuffer[count] = values[ptr];
            }
            if (needs_index) {
                ibuffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueStorage_& values;
    sparse_utils::BlockSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > cache;
    bool needs_value, needs_index;
};

/***********************
 *** Secondary index ***
 ***********************/

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicIndexDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, VectorPtr<Index_> sub_ptr) :
        values(vstore), cache(make_ServeIndices<Index_>(istore, pstore), sec, std::move(sub_ptr)) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + cache.size(), static_cast<Value_>(0));
        cache.search(i, [&](Index_, Index_ index_primary, ElementType<PointerStorage_> ptr) {
            buffer[index_primary] = values[ptr];
        });
        return buffer;
    }

private:
    const ValueStorage_& values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > cache;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
    SecondaryMyopicIndexSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, VectorPtr<Index_> sub_ptr, const Options& opt) :
        values(vstore), cache(make_ServeIndices<Index_>(istore, pstore), sec, std::move(sub_ptr)), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {} 

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        cache.search(i, [&](Index_ primary, Index_, ElementType<PointerStorage_> ptr) {
            if (needs_value) {
                vbuffer[count] = values[ptr];
            }
            if (needs_index) {
                ibuffer[count] = primary;
            }
            ++count;
        });
        return SparseRange<Value_, Index_>(count, needs_value ? vbuffer : NULL, needs_index ? ibuffer : NULL);
    }

private:
    const ValueStorage_& values;
    sparse_utils::IndexSecondaryExtractionCache<Index_, ServeIndices<Index_, IndexStorage_, PointerStorage_> > cache;
    bool needs_value, needs_index;
};

}
/**
 * @endcond
 */

/**
 * @brief Compressed sparse matrix representation.
 *
 * @tparam row_ Whether this is a compressed sparse row representation.
 * If `false`, a compressed sparse column representation is expected instead.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam ValueStorage_ Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `Value_`, as long as the type is convertible to `Value_`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Value_*`, it will also be used.
 * @tparam IndexStorage_ Vector class used to store the row/column indices internally.
 * This does not necessarily have to contain `Index_`, as long as the type is convertible to `Index_`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const Index_*`, it will also be used.
 * @tparam PointerStorage_ Vector class used to store the column/row index pointers.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 */
template<
    bool row_, 
    typename Value_, 
    typename Index_,
    class ValueStorage_ = std::vector<Value_>, 
    class IndexStorage_ = std::vector<Index_>, 
    class PointerStorage_ = std::vector<size_t> 
>
class CompressedSparseMatrix : public Matrix<Value_, Index_> {
public:
    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param vals Vector of non-zero elements.
     * @param idx Vector of row indices (if `row_ = false`) or column indices (if `row_ = true`) for the non-zero elements.
     * @param ptr Vector of index pointers.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `vals` and `idx` have the same length, equal to the number of structural non-zero elements;
     * `ptr` has length equal to the number of rows (if `row_ = true`) or columns (otherwise) plus one;
     * `ptr` is non-decreasing with first and last values set to 0 and the number of structural non-zeroes, respectively;
     * `idx` is strictly increasing within each interval defined by successive elements of `ptr`;
     * and all values of `idx` are non-negative and less than the number of columns (if `row_ = true`) or rows (otherwise).
     */
    CompressedSparseMatrix(Index_ nr, Index_ nc, ValueStorage_ vals, IndexStorage_ idx, PointerStorage_ ptr, bool check=true) : 
        nrows(nr), ncols(nc), values(std::move(vals)), indices(std::move(idx)), indptrs(std::move(ptr)) 
    {
        if (check) {
            if (values.size() != indices.size()) {
                throw std::runtime_error("'values' and 'indices' should be of the same length");
            }

            if constexpr(row_) {
                if (indptrs.size() != static_cast<size_t>(nrows) + 1){
                    throw std::runtime_error("length of 'indptrs' should be equal to 'nrows + 1'");
                }
            } else {
                if (indptrs.size() != static_cast<size_t>(ncols) + 1){
                    throw std::runtime_error("length of 'indptrs' should be equal to 'ncols + 1'");
                }
            }

            if (indptrs[0] != 0) {
                throw std::runtime_error("first element of 'indptrs' should be zero");
            }

            auto last = indptrs[indptrs.size() - 1]; // don't use back() as this is not guaranteed to be available for arbitrary PointerStorage_.
            if (static_cast<size_t>(last) != indices.size()) {
                throw std::runtime_error("last element of 'indptrs' should be equal to length of 'indices'");
            }

            ElementType<IndexStorage_> max_index = (row_ ? ncols : nrows);
            for (size_t i = 1; i < indptrs.size(); ++i) {
                auto start = indptrs[i- 1], end = indptrs[i];
                if (end < start || end > last) {
                    throw std::runtime_error("'indptrs' should be in non-decreasing order");
                }

                for (auto x = start; x < end; ++x) {
                    if (indices[x] < 0 || indices[x] >= max_index) {
                        if constexpr(row_) {
                            throw std::runtime_error("'indices' should contain non-negative integers less than the number of rows");
                        } else {
                            throw std::runtime_error("'indices' should contain non-negative integers less than the number of columns");
                        }
                    }
                }

                for (size_t j = start + 1; j < end; ++j) {
                    if (indices[j] <= indices[j - 1]) {
                        if constexpr(row_) {
                            throw std::runtime_error("'indices' should be strictly increasing within each row");
                        } else {
                            throw std::runtime_error("'indices' should be strictly increasing within each column");
                        }
                    }
                }
            }
        }
    }

private:
    Index_ nrows, ncols;
    ValueStorage_ values;
    IndexStorage_ indices;
    PointerStorage_ indptrs;

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool sparse() const { return true; }

    double sparse_proportion() const { return 1; }

    /**
     * @return `true` if `row_ = true` (for `CompressedSparseRowMatrix` objects), otherwise returns `false` (for `CompressedSparseColumnMatrix` objects).
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
public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, const Options&) const {
        if (row_ == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicFullDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary());
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicFullDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary()); 
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_end, const Options&) const {
        if (row_ == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicBlockDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), block_start, block_end);
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicBlockDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), block_start, block_end);
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(bool row, VectorPtr<Index_> subset_ptr, const Options&) const {
        if (row_ == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicIndexDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, std::move(subset_ptr));
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicIndexDense<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), std::move(subset_ptr));
        }
    }

    /******************************
     ******* Sparse myopic ********
     ******************************/
public:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const Options& opt) const {
        if (row_ == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicFullSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), opt);
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicFullSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), opt); 
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_end, const Options& opt) const {
        if (row_ == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicBlockSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), block_start, block_end, opt);
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicBlockSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), block_start, block_end, opt);
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(bool row, VectorPtr<Index_> subset_ptr, const Options& opt) const {
        if (row_ == row) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicIndexSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, std::move(subset_ptr), opt);
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::SecondaryMyopicIndexSparse<Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_> >(values, indices, indptrs, secondary(), std::move(subset_ptr), opt);
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

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(std::move(oracle), dense(row, std::move(indices_ptr), opt));
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

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const Oracle<Index_> > oracle, VectorPtr<Index_> indices_ptr, const Options& opt) const {
        return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(std::move(oracle), sparse(row, std::move(indices_ptr), opt));
    }
};

/**
 * Compressed sparse column matrix.
 * See `tatami::CompressedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueStorage_ = std::vector<Value_>, class IndexStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
using CompressedSparseColumnMatrix = CompressedSparseMatrix<false, Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_>;

/**
 * Compressed sparse row matrix.
 * See `tatami::CompressedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueStorage_ = std::vector<Value_>, class IndexStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
using CompressedSparseRowMatrix = CompressedSparseMatrix<true, Value_, Index_, ValueStorage_, IndexStorage_, PointerStorage_>;

}

#endif
