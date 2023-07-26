#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "primary_extraction.hpp"
#include "SparseSecondaryExtractorCore.hpp"

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

            Stored<IndexStorage_> max_index = (row_ ? ncols : nrows);
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
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_>
    struct CompressedExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        CompressedExtractorBase(const CompressedSparseMatrix* p, const Options& opt) : parent(p), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? parent->ncols : parent->nrows);
            }
        }

        CompressedExtractorBase(const CompressedSparseMatrix* p, const Options& opt, Index_ bs, Index_ bl) : CompressedExtractorBase(p, opt) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
        }

        CompressedExtractorBase(const CompressedSparseMatrix* p, const Options& opt, std::vector<Index_> i) : CompressedExtractorBase(p, opt) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                subset_indices = std::move(i);
                this->index_length = subset_indices.size();
            }
        }

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return subset_indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> >) {
            return;
        }

    protected:
        const CompressedSparseMatrix* parent;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type subset_indices;
        bool needs_value = false;
        bool needs_index = false;
    };

    /***********************************
     ******* Primary extraction ********
     ***********************************/
private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct PrimaryExtractorBase : public CompressedExtractorBase<row_, selection_, sparse_> {
        template<typename ...Args_>
        PrimaryExtractorBase(const CompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : CompressedExtractorBase<row_, selection_, sparse_>(p, opt, std::forward<Args_>(args)...) {
            bool spawn_cache = false;

            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                // only need to create a cache if block does not start at 0, see primary_dimension for details.
                spawn_cache = (opt.cache_for_reuse && this->block_start); 
            } else if constexpr(selection_ == DimensionSelectionType::INDEX) {
                // only need to create a cache if indices are non-empty and the first is not 0, see primary_dimension for details.
                spawn_cache = (opt.cache_for_reuse && this->index_length && this->subset_indices[0]);
            }

            if (spawn_cache) {
                auto len = row_ ? this->parent->nrows : this->parent->ncols;
                if constexpr(selection_ == DimensionSelectionType::INDEX) {
                    cached.resize(len, -1);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    cached.resize(len, std::pair<size_t, size_t>(-1, 0));
                }
            }
        }

    protected:
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, 
             std::vector<size_t>, 
             typename std::conditional<selection_ == DimensionSelectionType::BLOCK,
                std::vector<std::pair<size_t, size_t> >,
                bool
            >::type
        >::type cached;
    };

    template<DimensionSelectionType selection_>
    struct DensePrimaryExtractor : public PrimaryExtractorBase<selection_, false> {
        template<typename ...Args_>
        DensePrimaryExtractor(const CompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : PrimaryExtractorBase<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                auto obtained = sparse_utils::extract_primary_dimension(i, this->parent->indices, this->parent->indptrs);
                sparse_utils::transplant_primary_expanded(this->parent->values, this->parent->indices, obtained, buffer, static_cast<Index_>(0), this->full_length);

            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                auto obtained = sparse_utils::extract_primary_dimension(i, this->block_start, this->block_length, this->parent->indices, this->parent->indptrs, this->cached);
                sparse_utils::transplant_primary_expanded(this->parent->values, this->parent->indices, obtained, buffer, this->block_start, this->block_length);

            } else {
                std::fill(buffer, buffer + this->index_length, static_cast<Value_>(0));
                sparse_utils::SimpleExpandedStore<Value_, Index_, ValueStorage_> store(this->parent->values, buffer);
                sparse_utils::primary_dimension(i, this->subset_indices.data(), this->index_length, this->parent->indices, this->parent->indptrs, this->cached, store);
            }

            return buffer;
        }

    };

    template<DimensionSelectionType selection_>
    struct SparsePrimaryExtractor : public PrimaryExtractorBase<selection_, true> {
        template<typename ...Args_>
        SparsePrimaryExtractor(const CompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : PrimaryExtractorBase<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            if (!this->needs_value) {
                vbuffer = NULL;
            }
            if (!this->needs_index) {
                ibuffer = NULL;
            }

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                auto obtained = sparse_utils::extract_primary_dimension(i, this->parent->indices, this->parent->indptrs);
                SparseRange<Value_, Index_> output(obtained.second);
                sparse_utils::transplant_primary_values(this->parent->values, obtained, output, vbuffer);
                sparse_utils::transplant_primary_indices(this->parent->indices, obtained, output, ibuffer);
                return output;

            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                auto obtained = sparse_utils::extract_primary_dimension(i, this->block_start, this->block_length, this->parent->indices, this->parent->indptrs, this->cached);
                SparseRange<Value_, Index_> output(obtained.second);
                sparse_utils::transplant_primary_values(this->parent->values, obtained, output, vbuffer);
                sparse_utils::transplant_primary_indices(this->parent->indices, obtained, output, ibuffer);
                return output;

            } else {
                sparse_utils::SimpleRawStore<Value_, Index_, ValueStorage_> store(this->parent->values, vbuffer, ibuffer);
                sparse_utils::primary_dimension(i, this->subset_indices.data(), this->index_length, this->parent->indices, this->parent->indptrs, this->cached, store);
                return SparseRange<Value_, Index_>(store.n, vbuffer, ibuffer);
            }
        }
    };

    /*************************************
     ******* Secondary extraction ********
     *************************************/
private:
    typedef Stored<IndexStorage_> StoredIndex;
    typedef Stored<PointerStorage_> StoredPointer;

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
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_> 
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options& opt, Args_&& ... args) const { 
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(accrow_ == row_) {
            if constexpr(sparse_) {
                output.reset(new SparsePrimaryExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            } else {
                output.reset(new DensePrimaryExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            }
        } else {
            if constexpr(sparse_) {
                output.reset(new SparseSecondaryExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            } else {
                output.reset(new DenseSecondaryExtractor<selection_>(this, opt, std::forward<Args_>(args)...));
            }
        }

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return populate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return populate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
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
