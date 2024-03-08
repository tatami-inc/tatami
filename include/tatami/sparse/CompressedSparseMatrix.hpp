#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "primary_extraction.hpp"
#include "SparseSecondaryExtractorCore.hpp"
#include "../utils/ElementType.hpp"
#include "../utils/OracleUnawareMatrix.hpp"

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

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicFullDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto offset = indptr[i];
        auto vIt = values.begin() + offset;
        auto iIt = indices.begin() + offset;
        auto delta = indptr[i+1] - indptr[i];

        std::fill(buffer, buffer + secondary, static_cast<Value_>(0));
        for (size_t x = 0; x < delta; ++x, ++vIt, ++iIt) {
            buffer[*iIt] = *vIt;
        }
        return buffer;
    }

    Index_ number() const {
        return secondary;
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

    const Value_* fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
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

    Index_ number() const {
        return secondary;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
    bool needs_value, needs_index;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicBlockDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicBlockDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec, Index_ bs, Index_ bl) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec), block_start(bs), block_length(bl) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto iStart = indices.begin() + indptr[i];
        auto iEnd = indices.begin() + indptr[i + 1];
        sparse_utils::refine_primary_block_limits(iStart, iEnd, secondary, block_start, block_end);

        std::fill(buffer, buffer + block_length, static_cast<Value_>(0));
        auto vIt = values.begin() + (iStart - indices.begin());
        for (; iStart != iEnd; ++iStart, ++vIt) {
            buffer[*iStart - block_start] = *vIt;
        }
        return buffer;
    }

    Index_ number() const {
        return block_length;
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

    const Value_* fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto iStart = indices.begin() + indptr[i];
        auto iEnd = indices.begin() + indptr[i + 1];
        sparse_utils::refine_primary_block_limits(iStart, iEnd, secondary, block_start, block_end);
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

    Index_ number() const {
        return block_length;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
    Index_ block_start, block_length;
    bool needs_value, needs_index;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicIndexDense : public MyopicDenseExtractor<Value_, Index_> {
    PrimaryMyopicIndexDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, std::vector<Index_> sub) :
        values(vstore), indices(istore), indptr(pstore), subset(std::move(sub)) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        std::fill(buffer, buffer + subset.size(), static_cast<Value_>(0));
        auto vIt = values.begin() + indptr[i];
        retrieve_primary_subset(
            indices.begin() + indptr[i], 
            indices.begin() + indptr[i+1],
            subset,
            [&](size_t s, size_t offset, Index_) {
                buffer[s] = *(vIt + offset);
            }
        );
        return buffer;
    }

    Index_ number() const {
        return indices.size();
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    std::vector<Index_> subset;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct PrimaryMyopicIndexSparse : public MyopicSparseExtractor<Value_, Index_> {
    PrimaryMyopicIndexSparse(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, std::vector<Index_> sub, const Options& opt) :
        values(vstore), indices(istore), indptr(pstore), subset(std::move(sub)), needs_value(opts.sparse_extract_value), needs_index(opts.sparse_extract_index) {} 

    const Value_* fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ count = 0;
        auto vcopy = vbuffer;
        auto icopy = ibuffer;

        auto vIt = values.begin() + indptr[i];
        retrieve_primary_subset(
            indices.begin() + indptr[i], 
            indices.begin() + indptr[i+1],
            subset,
            [&](size_t, size_t offset, Index_ ix) {
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

    Index_ number() const {
        return subset.size();
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    std::vector<Index_> subset;
    bool needs_value, needs_index;
};

template<typename Value_, typename Index_, class ValueStorage_, class IndexStorage_, class PointerStorage_>
struct SecondaryMyopicFullDense : public MyopicDenseExtractor<Value_, Index_> {
    SecondaryMyopicFullDense(const ValueStorage_& vstore, const IndexStorage_& istore, const PointerStorage_& pstore, Index_ sec) :
        values(vstore), indices(istore), indptr(pstore), secondary(sec) {} 

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto offset = indptr[i];
        auto vIt = values.begin() + offset;
        auto iIt = indices.begin() + offset;
        auto delta = indptr[i+1] - indptr[i];

        std::fill(buffer, buffer + secondary, static_cast<Value_>(0));
        for (size_t x = 0; x < delta; ++x, ++vIt, ++iIt) {
            buffer[*iIt] = *vIt;
        }
        return buffer;
    }

    Index_ number() const {
        return indptr.size() - 1;
    }

private:
    const ValueStorage_& values;
    const IndexStorage_& indices;
    const PointerStorage_& indptr;
    Index_ secondary;
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
    typename Index_ = int, 
    class ValueStorage_ = std::vector<Value_>, 
    class IndexStorage_ = std::vector<Index_>, 
    class PointerStorage_ = std::vector<size_t> 
>
class CompressedSparseMatrix : public OracleUnawareMatrix<Value_, Index_> {
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

    /*************************************
     ******* Secondary extraction ********
     *************************************/
private:
    typedef ElementType<IndexStorage_> StoredIndex;
    typedef ElementType<PointerStorage_> StoredPointer;

    struct SecondaryModifier {
        static void increment(StoredPointer& ptr, const IndexStorage_&, StoredPointer) { ++ptr; }
        static void decrement(StoredPointer& ptr, const IndexStorage_&, StoredPointer) { --ptr; }
        static StoredPointer get(StoredPointer ptr) { return ptr; }
        static void set(StoredPointer& ptr, StoredPointer val) { ptr = val; }
    };

    struct SecondaryCore : public SparseSecondaryExtractorCore<Index_, StoredIndex, StoredPointer, SecondaryModifier> {
        SecondaryCore() = default;

        SecondaryCore(StoredIndex max_index, const IndexStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
            SparseSecondaryExtractorCore<Index_, StoredIndex, StoredPointer, SecondaryModifier>(max_index, length)
        {
            auto idpIt = idp.begin() + start;
            for (Index_ i = 0; i < length; ++i, ++idpIt) {
                this->current_indptrs[i] = *idpIt;
                this->current_indices[i] = (*idpIt < *(idpIt + 1) ? idx[*idpIt] : max_index);
            }
            this->closest_current_index = (length ? *std::min_element(this->current_indices.begin(), this->current_indices.end()) : max_index);
            return;
        } 

        SecondaryCore(StoredIndex max_index, const IndexStorage_& idx, const PointerStorage_& idp) :
            SecondaryCore(max_index, idx, idp, static_cast<Index_>(0), static_cast<Index_>(idp.size() - 1)) {}

        SecondaryCore(StoredIndex max_index, const IndexStorage_& idx, const PointerStorage_& idp, const Index_* subset, Index_ length) :
            SparseSecondaryExtractorCore<Index_, StoredIndex, StoredPointer, SecondaryModifier>(max_index, length)
        {
            for (Index_ i0 = 0; i0 < length; ++i0) {
                auto i = subset[i0];
                this->current_indptrs[i0] = idp[i];
                this->current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
            }
            this->closest_current_index = (length ? *std::min_element(this->current_indices.begin(), this->current_indices.end()) : max_index);
            return;
        }

        template<class PrimaryFunction_, class StoreFunction_, class SkipFunction_>
        bool search(StoredIndex secondary, Index_ primary_length, PrimaryFunction_&& to_primary, const IndexStorage_& indices, const PointerStorage_& indptrs, StoreFunction_&& store, SkipFunction_&& skip) {
            return this->search_base(
                secondary, 
                primary_length, 
                std::forward<PrimaryFunction_>(to_primary), 
                indices, 
                indptrs, 
                std::forward<StoreFunction_>(store), 
                std::forward<SkipFunction_>(skip)
            );
        }
    };

    template<DimensionSelectionType selection_, bool sparse_>
    struct SecondaryExtractorBase : public CompressedExtractorBase<!row_, selection_, sparse_> {
        template<typename ...Args_>
        SecondaryExtractorBase(const CompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : CompressedExtractorBase<!row_, selection_, sparse_>(p, opt, std::forward<Args_>(args)...) {
            auto max_index = (row_ ? this->parent->ncols : this->parent->nrows);

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                state = SecondaryCore(max_index, this->parent->indices, this->parent->indptrs);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                state = SecondaryCore(max_index, this->parent->indices, this->parent->indptrs, this->block_start, this->block_length);
            } else {
                state = SecondaryCore(max_index, this->parent->indices, this->parent->indptrs, this->subset_indices.data(), this->index_length);
            }
        }

    private:
        SecondaryCore state;

    protected:
        template<class Store_>
        void secondary_dimension_loop(Index_ i, Index_ start, Index_ length, Store_& store) {
            state.search(
                i, 
                length, 
                [&](Index_ p) -> Index_ { 
                    return p + start; 
                },
                this->parent->indices,
                this->parent->indptrs,
                [&](Index_ primary, StoredPointer curptr) -> void {
                    store.add(primary, curptr);
                },
                [&](Index_ primary) -> void {
                    store.skip(primary);
                }
            );
        }

        template<class Store_>
        void secondary_dimension_loop(Index_ i, const Index_* subset, Index_ length, Store_& output) {
            state.search(
                i, 
                length, 
                [&](Index_ p) -> Index_ { 
                    return subset[p];
                },
                this->parent->indices,
                this->parent->indptrs,
                [&](Index_ primary, StoredPointer curptr) -> void {
                    output.add(primary, curptr);
                },
                [&](Index_ primary) -> void {
                    output.skip(primary);
                }
            );
            return;
        }
    };

    template<DimensionSelectionType selection_>
    struct DenseSecondaryExtractor : public SecondaryExtractorBase<selection_, false> {
        template<typename ...Args_>
        DenseSecondaryExtractor(const CompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : SecondaryExtractorBase<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

    private:
        struct ExpandedStoreBlock {
            ExpandedStoreBlock(const ValueStorage_& iv, Value_* ov) : in_values(iv), out_values(ov) {}
            const ValueStorage_& in_values;
            Value_* out_values;
            Index_ first;

            void add(Index_ i, StoredPointer ptr) {
                out_values[i - first] = in_values[ptr];
                return;
            }

            void skip(Index_) {} 
        };

        struct ExpandedStoreIndexed {
            ExpandedStoreIndexed(const ValueStorage_& iv, Value_* ov) : in_values(iv), out_values(ov) {}
            const ValueStorage_& in_values;
            Value_* out_values;

            void add(Index_, StoredPointer ptr) {
                *out_values = in_values[ptr];
                ++out_values;
                return;
            }

            void skip(Index_) {
                ++out_values;
            } 
        };

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            typename std::conditional<selection_ == DimensionSelectionType::INDEX, ExpandedStoreIndexed, ExpandedStoreBlock>::type store(this->parent->values, buffer);
            std::fill(buffer, buffer + extracted_length<selection_, Index_>(*this), static_cast<Value_>(0));

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                store.first = 0;
                this->secondary_dimension_loop(i, static_cast<Index_>(0), this->full_length, store);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                store.first = this->block_start;
                this->secondary_dimension_loop(i, this->block_start, this->block_length, store);
            } else {
                this->secondary_dimension_loop(i, this->subset_indices.data(), this->index_length, store);
            }

            return buffer;
        }
    };

    template<DimensionSelectionType selection_>
    struct SparseSecondaryExtractor : public SecondaryExtractorBase<selection_, true> {
        template<typename ...Args_>
        SparseSecondaryExtractor(const CompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : SecondaryExtractorBase<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            if (!this->needs_value) {
                vbuffer = NULL;
            }
            if (!this->needs_index) {
                ibuffer = NULL;
            }

            sparse_utils::SimpleRawStore<Value_, Index_, ValueStorage_> store(this->parent->values, vbuffer, ibuffer);
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->secondary_dimension_loop(i, static_cast<Index_>(0), this->full_length, store);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->secondary_dimension_loop(i, this->block_start, this->block_length, store);
            } else {
                this->secondary_dimension_loop(i, this->subset_indices.data(), this->index_length, store);
            }

            return SparseRange<Value_, Index_>(store.n, vbuffer, ibuffer);
        }
    };

    /*************************************
     ******* Extraction overrides ********
     *************************************/
private:
    Index_ secondary() const {
        if constexpr(row_) {
            return ncols;
        } else {
            return nrows;
        }
    }

    Index_ primary() const {
        if constexpr(row_) {
            return nrows;
        } else {
            return ncols;
        }
    }

public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(const Options&) const {
        if constexpr(row_) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicFullDense<Value_, Index_> >(values, indices, indptrs, secondary());
        } else {
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options&) const {
        if constexpr(row_) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicBlockDense<Value_, Index_> >(values, indices, indptrs, secondary(), block_start, block_end);
        } else {

        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> subset, const Options& opt) const {
        if constexpr(row_) {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicIndexDense<Value_, Index_> >(values, indices, indptrs, std::move(subset));
        } else {
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        if constexpr(row_) {
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicFullDense<Value_, Index_> >(values, indices, indptrs, secondary());
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        if constexpr(row_) {
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicBlockDense<Value_, Index_> >(values, indices, indptrs, secondary(), block_start, block_end);
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        if constexpr(row_) {
        } else {
            return std::make_unique<CompressedSparseMatrix_internal::PrimaryMyopicIndexDense<Value_, Index_> >(values, indices, indptrs, std::move(subset));
        }
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
