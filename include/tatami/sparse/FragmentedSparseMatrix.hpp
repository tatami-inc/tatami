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
#include <string>

/**
 * @file FragmentedSparseMatrix.hpp
 *
 * @brief Fragmented sparse matrix representation. 
 *
 * `typedef`s are provided for the usual row and column formats. 
 */

namespace tatami {

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
    class ValueVectorStorage_ = std::vector<std::vector<Value_> >
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
        check_values(check); 
        return;
    }

private:
    Index_ nrows, ncols;
    ValueVectorStorage_ values;
    IndexVectorStorage_ indices;

    void check_values(bool check) {
        if (!check) {
            return;
        }

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

        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            auto start = indptrs[i- 1], end = indptrs[i];
            if (end < start || end > last) {
                throw std::runtime_error("'indptrs' should be in non-decreasing order");
            }

            for (size_t x = 1, xend = currenti.size(); x < xend; ++x) {
                if (previous >= indices[x]) {
                    std::string msg = "'indices' should be strictly increasing within each ";
                    if constexpr(row_) {
                        msg += "row";
                    } else {
                        msg += "column";
                    }
                    throw std::runtime_error(msg);
                }
            }
        }

        return;
    }

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool sparse() const { return true; }

    double sparse_proportion() const { return 1; }

    /**
     * @return `true` if `ROW = true` (for `FragmentedSparseRowMatrix` objects), otherwise returns `false` (for `FragmentedSparseColumnMatrix` objects).
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
    struct FragmentedExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        FragmentedExtractorBase(const FragmentedSparseMatrix* p, const Options& opt) : parent(p), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? parent->ncols : parent->nrows);
            }
        }

        FragmentedExtractorBase(const FragmentedSparseMatrix* p, const Options& opt, Index_ bs, Index_ bl) : FragmentedExtractorBase(p, opt) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
        }

        FragmentedExtractorBase(const FragmentedSparseMatrix* p, const Options& opt, std::vector<Index_> i) : FragmentedExtractorBase(p, opt) {
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
        const FragmentedSparseMatrix* parent;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type subset_indices;
        bool needs_value = false;
        bool needs_index = false;
    };

    /***********************************
     ******* Primary extraction ********
     ***********************************/
private:
    struct RawStore {
        RawStore(const ValueVectorStorage_& iv, Value_* ov, Index_* oi) : in_values(iv), out_values(ov), out_indices(oi) {}

    private:
        const ValueVectorStorage_& in_values;
        Value_* out_values;
        Index_* out_indices;

    public:
        Index_ n = 0;

        void add(Index_ i, size_t ptr) {
            ++n;
            if (out_indices) {
                *out_indices = i;
                ++out_indices;
            }
            if (out_values) {
                *out_values = in_values[ptr];
                ++out_values;
            }
            return;
        }

        void skip(Index_) {} 
    };

private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct PrimaryExtractorBase : public FragmentedExtractorBase<row_, selection_, sparse_> {
        template<typename ...Args_>
        PrimaryExtractorBase(const FragmentedSparseMatrix* p, const Options& opt, Args_&& ... args) : FragmentedExtractorBase<row_, selection_, sparse_>(p, opt, std::forward<Args_>(args)...) {
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
        DensePrimaryExtractor(const FragmentedSparseMatrix* p, const Options& opt, Args_&& ... args) : PrimaryExtractorBase<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

    private:
        struct ExpandedStore {
            ExpandedStore(const ValueVectorStorage_& iv, Value_* ov) : in_values(iv), out_values(ov) {}

        private:
            const ValueVectorStorage_& in_values;
            Value_* out_values;

        public:
            void add(Index_, size_t ptr) {
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
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                auto obtained = sparse_utils::extract_primary_dimension(i, this->parent->indices, this->parent->indptrs);
                sparse_utils::transplant_primary_expanded(this->parent->values, this->parent->indices, obtained, buffer, static_cast<Index_>(0), this->full_length);

            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                auto obtained = sparse_utils::extract_primary_dimension(i, this->block_start, this->block_length, this->parent->indices, this->parent->indptrs, this->cached);
                sparse_utils::transplant_primary_expanded(this->parent->values, this->parent->indices, obtained, buffer, this->block_start, this->block_length);

            } else {
                std::fill(buffer, buffer + this->index_length, static_cast<Value_>(0));
                ExpandedStore store(this->parent->values, buffer);
                sparse_utils::primary_dimension(i, this->subset_indices.data(), this->index_length, this->parent->indices, this->parent->indptrs, this->cached, store);
            }

            return buffer;
        }

    };

    template<DimensionSelectionType selection_>
    struct SparsePrimaryExtractor : public PrimaryExtractorBase<selection_, true> {
        template<typename ...Args_>
        SparsePrimaryExtractor(const FragmentedSparseMatrix* p, const Options& opt, Args_&& ... args) : PrimaryExtractorBase<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

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
                RawStore store(this->parent->values, vbuffer, ibuffer);
                sparse_utils::primary_dimension(i, this->subset_indices.data(), this->index_length, this->parent->indices, this->parent->indptrs, this->cached, store);
                return SparseRange<Value_, Index_>(store.n, vbuffer, ibuffer);
            }
        }
    };

    /*************************************
     ******* Secondary extraction ********
     *************************************/
private:
    typedef Stored<IndexVectorStorage_> StoredIndex;
    typedef Stored<PointerStorage_> StoredPointer;

    struct SecondaryModifier {
        static void increment(StoredPointer& ptr, const IndexVectorStorage_&, StoredPointer) { ++ptr; }
        static void decrement(StoredPointer& ptr, const IndexVectorStorage_&, StoredPointer) { --ptr; }
        static StoredPointer get(StoredPointer ptr) { return ptr; }
        static void set(StoredPointer& ptr, StoredPointer val) { ptr = val; }
    };

    struct SecondaryCore : public SparseSecondaryExtractorCore<Index_, StoredIndex, StoredPointer, SecondaryModifier> {
        SecondaryCore() = default;

        SecondaryCore(StoredIndex max_index, const IndexVectorStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
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

        SecondaryCore(StoredIndex max_index, const IndexVectorStorage_& idx, const PointerStorage_& idp) :
            SecondaryCore(max_index, idx, idp, static_cast<Index_>(0), static_cast<Index_>(idp.size() - 1)) {}

        SecondaryCore(StoredIndex max_index, const IndexVectorStorage_& idx, const PointerStorage_& idp, const Index_* subset, Index_ length) :
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
        bool search(StoredIndex secondary, Index_ primary_length, PrimaryFunction_&& to_primary, const IndexVectorStorage_& indices, const PointerStorage_& indptrs, StoreFunction_&& store, SkipFunction_&& skip) {
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
    struct SecondaryExtractorBase : public FragmentedExtractorBase<!row_, selection_, sparse_> {
        template<typename ...Args_>
        SecondaryExtractorBase(const FragmentedSparseMatrix* p, const Options& opt, Args_&& ... args) : FragmentedExtractorBase<!row_, selection_, sparse_>(p, opt, std::forward<Args_>(args)...) {
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
        DenseSecondaryExtractor(const FragmentedSparseMatrix* p, const Options& opt, Args_&& ... args) : SecondaryExtractorBase<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

    private:
        struct ExpandedStoreBlock {
            const ValueVectorStorage_* in_values;
            Value_* out_values;
            Index_ first;

            void add(Index_ i, StoredPointer ptr) {
                out_values[i - first] = (*in_values)[ptr];
                return;
            }

            void skip(Index_) {} 
        };

        struct ExpandedStoreIndexed {
            const ValueVectorStorage_* in_values;
            Value_* out_values;

            void add(Index_, StoredPointer ptr) {
                *out_values = (*in_values)[ptr];
                ++out_values;
                return;
            }

            void skip(Index_) {
                ++out_values;
            } 
        };

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            typename std::conditional<selection_ == DimensionSelectionType::INDEX, ExpandedStoreIndexed, ExpandedStoreBlock>::type store;
            store.in_values = &(this->parent->values);
            store.out_values = buffer;
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
        SparseSecondaryExtractor(const FragmentedSparseMatrix* p, const Options& opt, Args_&& ... args) : SecondaryExtractorBase<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            if (!this->needs_value) {
                vbuffer = NULL;
            }
            if (!this->needs_index) {
                ibuffer = NULL;
            }

            RawStore store(this->parent->values, vbuffer, ibuffer);
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
 * Fragmented sparse column matrix.
 * See `tatami::FragmentedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_ = std::vector<Value_>, class IndexVectorStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
using FragmentedSparseColumnMatrix = FragmentedSparseMatrix<false, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_, PointerStorage_>;

/**
 * Fragmented sparse row matrix.
 * See `tatami::FragmentedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class ValueVectorStorage_ = std::vector<Value_>, class IndexVectorStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
using FragmentedSparseRowMatrix = FragmentedSparseMatrix<true, Value_, Index_, ValueVectorStorage_, IndexVectorStorage_, PointerStorage_>;

}

#endif
