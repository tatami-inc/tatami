#ifndef TATAMI_SEMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_SEMI_COMPRESSED_SPARSE_MATRIX_H

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "utils.hpp"
#include "SparseSecondaryExtractorCore.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>
#include <string>

/**
 * @file SemiCompressedSparseMatrix.hpp
 *
 * @brief Semi-compressed sparse matrix representation. 
 *
 * `typedef`s are provided for the usual row and column formats. 
 */

namespace tatami {

/**
 * @brief Semi-compressed sparse matrix representation.
 *
 * This refers to a sparse matrix of non-negative integers, where counts greater than 1 are represented by duplicated indices.
 * For data where the vast majority of counts are 1, the semi-compressed representation reduces memory usage by eliminating the storage of the values.
 * Nonetheless, it is still capable of handling larger counts on the rare occasions that they do occur.
 *
 * @tparam row_ Whether this is a compressed sparse row representation.
 * If `false`, a compressed sparse column representation is used instead.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
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
    class IndexStorage_ = std::vector<Index_>, 
    class PointerStorage_ = std::vector<size_t> 
>
class SemiCompressedSparseMatrix : public Matrix<Value_, Index_> {
public:
    /**
     * @param nr Number of rows.
     * @param nc Number of columns.
     * @param idx Vector of row indices (if `row_ = false`) or column indices (if `row_ = true`) for the non-zero elements.
     * @param ptr Vector of index pointers.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `ptr` has length equal to the number of rows (if `row_ = true`) or columns (otherwise) plus one;
     * `ptr` is ordered with first and last values set to 0 and the length of `idx`, respectively;
     * `idx` is non-decreasing within each interval defined by successive elements of `ptr`;
     * and all values of `idx` are non-negative and less than the number of columns (if `row_ = true`) or rows (otherwise).
     */
    SemiCompressedSparseMatrix(Index_ nr, Index_ nc, IndexStorage_ idx, PointerStorage_ ptr, bool check=true) : 
        nrows(nr), ncols(nc), indices(std::move(idx)), indptrs(std::move(ptr)) 
    {
        if (check) {
            if (row_) {
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

            for (size_t i = 1; i < indptrs.size(); ++i) {
                auto start = indptrs[i- 1], end = indptrs[i];
                if (end < start || end > last) {
                    throw std::runtime_error("'indptrs' should be in non-decreasing order");
                }

                for (auto x = start + 1; x < end; ++x) {
                    if (indices[x - 1] > indices[x]) {
                        std::string msg = "'indices' should be non-decreasing within each ";
                        if constexpr(row_) {
                            msg += "row";
                        } else {
                            msg += "column";
                        }
                        throw std::runtime_error(msg);
                    }
                }
            }
        }
    }

private:
    Index_ nrows, ncols;
    IndexStorage_ indices;
    PointerStorage_ indptrs;

    typedef Stored<IndexStorage_> StoredIndex;
    typedef Stored<PointerStorage_> StoredPointer;

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool sparse() const { return true; }

    double sparse_proportion() const { return 1; }

    /**
     * @return `true` if `ROW = true` (for `CompressedSparseRowMatrix` objects), otherwise returns `false` (for `CompressedSparseColumnMatrix` objects).
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
    struct SemiCompressedExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        SemiCompressedExtractorBase(const SemiCompressedSparseMatrix* p, const Options& opt) : parent(p), needs_value(opt.sparse_extract_value), needs_index(opt.sparse_extract_index) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = (accrow_ ? parent->ncols : parent->nrows);
            }
        }

        SemiCompressedExtractorBase(const SemiCompressedSparseMatrix* p, const Options& opt, Index_ bs, Index_ bl) : SemiCompressedExtractorBase(p, opt) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
        }

        SemiCompressedExtractorBase(const SemiCompressedSparseMatrix* p, const Options& opt, std::vector<Index_> i) : SemiCompressedExtractorBase(p, opt) {
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
        const SemiCompressedSparseMatrix* parent;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type subset_indices;
        bool needs_value = false;
        bool needs_index = false;
    };

    /***********************************
     ******* Primary extraction ********
     ***********************************/
private:
    struct RawStore {
        Value_* out_values;
        Index_* out_indices;
        Index_ n = 0;

        void add(Index_ i, Value_ val) {
            ++n;
            if (out_indices) {
                *out_indices = i;
                ++out_indices;
            }
            if (out_values) {
                *out_values = val;
                ++out_values;
            }
            return;
        }

        void skip(Index_) {} 
    };

private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct PrimaryExtractorBase : public SemiCompressedExtractorBase<row_, selection_, sparse_> {
        template<typename ...Args_>
        PrimaryExtractorBase(const SemiCompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : SemiCompressedExtractorBase<row_, selection_, sparse_>(p, opt, std::forward<Args_>(args)...) {
            bool spawn_cache = false;

            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                // only need to create a cache if block does not start at 0, see primary_dimension for details.
                spawn_cache = (opt.cache_for_reuse && this->block_start); 
            } else if constexpr(selection_ == DimensionSelectionType::INDEX) {
                // only need to create a cache if indices are non-empty and the first is not 0, see primary_dimension for details.
                spawn_cache = (opt.cache_for_reuse && this->index_length && this->subset_indices[0]);
            }

            if (spawn_cache) {
                cached.resize(row_ ? this->parent->nrows : this->parent->ncols, -1);
            }
        }

    private:
        std::vector<size_t> cached;

    protected:
        StoredPointer primary_start(Index_ i, Index_ start) {
            bool do_cache = !cached.empty();
            if (do_cache) {
                auto val = cached[i];
                if (val != -1) {
                    return val;
                }
            }

            auto outstart = this->parent->indptrs[i];
            if (start) { // Jumping ahead if non-zero.
                auto iIt = this->parent->indices.begin() + outstart;
                auto eIt = this->parent->indices.begin() + this->parent->indptrs[i + 1]; 
                outstart = std::lower_bound(iIt, eIt, start) - this->parent->indices.begin();
            } 

            if (do_cache) {
                cached[i] = outstart;
            }

            return outstart;
        }

        template<class Store_>
        void primary_dimension(Index_ i, const Index_* subset, Index_ length, Store_& store) {
            if (!length) {
                return;
            }

            auto iIt = this->parent->indices.begin() + this->parent->indptrs[i];
            auto eIt = this->parent->indices.begin() + this->parent->indptrs[i + 1]; 

            if (subset[0]) { // Only jumping ahead if the start is non-zero.
                bool do_cache = !cached.empty();
                if (do_cache) {
                    if (cached[i] != -1) { // retrieving the jump from cache, if we came here before.
                        iIt += cached[i];
                    } else {
                        auto iIt2 = std::lower_bound(iIt, eIt, subset[0]);
                        cached[i] = iIt2 - iIt;
                        iIt = iIt2;
                    }
                } else {
                    iIt = std::lower_bound(iIt, eIt, subset[0]);
                }
            } 

            if (iIt == eIt) {
                return;
            }

            for (Index_ counter = 0; counter < length; ++counter) {
                auto current = subset[counter];

                while (iIt != eIt && current > *iIt) {
                    ++iIt;
                }
                if (iIt == eIt) {
                    break;
                }

                if (current == *iIt) {
                    int count = 0;
                    do {
                        ++iIt;
                        ++count;
                    } while (iIt != eIt && current == *iIt);
                    store.add(current, count);
                } else {
                    store.skip(current);
                }
            }

            return;
        }
    };

private:
    template<DimensionSelectionType selection_>
    struct DensePrimaryExtractor : public PrimaryExtractorBase<selection_, false> {
        template<typename ...Args_>
        DensePrimaryExtractor(const SemiCompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : PrimaryExtractorBase<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

    private:
        struct ExpandedStore {
            Value_* out_values;

            void add(Index_, Value_ val) {
                *out_values = val;
                ++out_values;
                return;
            }

            void skip(Index_) {
                ++out_values;
            }
        };

        void primary_dimension_expanded(Index_ i, Index_ start, Index_ length, Value_* out_values) {
            std::fill(out_values, out_values + length, static_cast<Value_>(0));
            auto position = this->primary_start(i, start);
            auto limit = this->parent->indptrs[i + 1];
            auto end = start + length;

            const auto& indices = this->parent->indices;
            while (position < limit && indices[position] < end) {
                auto curdex = indices[position];
                auto copy = position;
                ++copy;
                for (; copy < limit && indices[copy] == curdex; ++copy) {}
                out_values[curdex - start] = copy - position;
                position = copy;
            }

            return;
        }

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                primary_dimension_expanded(i, static_cast<Index_>(0), this->full_length, buffer);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                primary_dimension_expanded(i, this->block_start, this->block_length, buffer);
            } else {
                std::fill(buffer, buffer + this->index_length, static_cast<Value_>(0));
                ExpandedStore store;
                store.out_values = buffer;
                this->primary_dimension(i, this->subset_indices.data(), this->index_length, store);
            }

            return buffer;
        }
    };

private:
    template<DimensionSelectionType selection_>
    struct SparsePrimaryExtractor : public PrimaryExtractorBase<selection_, true> {
        template<typename ...Args_>
        SparsePrimaryExtractor(const SemiCompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : PrimaryExtractorBase<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

    private:
        SparseRange<Value_, Index_> primary_dimension_raw(Index_ i, Index_ start, Index_ length, Value_* out_values, Index_* out_indices) {
            auto position = this->primary_start(i, start);
            auto limit = this->parent->indptrs[i + 1];
            auto end = start + length;

            SparseRange<Value_, Index_> output;
            output.number = 0;
            output.value = out_values;
            output.index = out_indices;
            auto& count = output.number;

            const auto& indices = this->parent->indices;
            while (position < limit && indices[position] < end) {
                auto curdex = indices[position];
                auto copy = position;
                ++copy;
                for (; copy < limit && indices[copy] == curdex; ++copy) {}

                if (out_values) {
                    *out_values = copy - position;
                    ++out_values;
                }
                if (out_indices) {
                    *out_indices = curdex;
                    ++out_indices;
                }

                ++count;
                position = copy;
            }

            return output;
        }

    public:
        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            if (!this->needs_value) {
                vbuffer = NULL;
            }
            if (!this->needs_index) {
                ibuffer = NULL;
            }

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                return primary_dimension_raw(i, static_cast<Index_>(0), this->full_length, vbuffer, ibuffer);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                return primary_dimension_raw(i, this->block_start, this->block_length, vbuffer, ibuffer);
            } else {
                RawStore store;
                store.out_values = vbuffer;
                store.out_indices = ibuffer;
                this->primary_dimension(i, this->subset_indices.data(), this->index_length, store);
                return SparseRange<Value_, Index_>(store.n, vbuffer, ibuffer);
            }
        }
    };

    /*************************************
     ******* Secondary extraction ********
     *************************************/
private:
    struct Position {
        Position() = default;
        Position(StoredPointer p) : pointer(p) {}
        StoredPointer pointer;
        Index_ count = 0;
        bool scanned = false;
    };

    static void update_secondary_position(Position& current, const IndexStorage_& indices, StoredPointer limit) {
        if (current.scanned) {
            return;
        }

        auto curdex = indices[current.pointer];
        auto copy = current.pointer;
        do {
            ++copy;
        } while (copy < limit && indices[copy] == curdex);
        current.count = copy - current.pointer;
        current.scanned = true;
        return;
    }

    struct SecondaryModifier {
        static void increment(Position& ptr, const IndexStorage_& indices, StoredPointer limit) {
            update_secondary_position(ptr, indices, limit);
            ptr.pointer += ptr.count;
            ptr.scanned = false;
            ptr.count = 0;
        }

        static void decrement(Position& ptr, const IndexStorage_& indices, StoredPointer limit) {
            if (ptr.pointer == limit) {
                return;
            }

            auto copy = ptr.pointer;
            --copy;
            auto curdex = indices[copy];
            while (copy > limit && indices[copy - 1] == curdex) {
                --copy;
            }

            ptr.count = ptr.pointer - copy;
            ptr.scanned = true;
            ptr.pointer = copy;
        }

        static StoredPointer get(const Position& ptr) {
            return ptr.pointer;
        }

        static void set(Position& ptr, StoredPointer val) {
            ptr.pointer = val;
            ptr.scanned = false;
            ptr.count = 0;
        }
    };

    struct SecondaryCore : public SparseSecondaryExtractorCore<Index_, StoredIndex, Position, SecondaryModifier> {
    public:
        SecondaryCore() = default;

        SecondaryCore(StoredIndex max_index, const IndexStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
            SparseSecondaryExtractorCore<Index_, StoredIndex, Position, SecondaryModifier>(max_index, length)
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
            SparseSecondaryExtractorCore<Index_, StoredIndex, Position, SecondaryModifier>(max_index, length)
        {
            for (Index_ i0 = 0; i0 < length; ++i0) {
                auto i = subset[i0];
                this->current_indptrs[i0] = idp[i];
                this->current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
            }
            this->closest_current_index = (length ? *std::min_element(this->current_indices.begin(), this->current_indices.end()) : max_index);
            return;
        }

    public:
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

private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct SecondaryExtractorBase : public SemiCompressedExtractorBase<!row_, selection_, sparse_> {
        template<typename ...Args_>
        SecondaryExtractorBase(const SemiCompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : SemiCompressedExtractorBase<!row_, selection_, sparse_>(p, opt, std::forward<Args_>(args)...) {
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
        void secondary_dimension_loop(Index_ i, Index_ start, Index_ length, Store_& output) {
            state.search(
                i, 
                length, 
                [&](Index_ p) -> Index_ { 
                    return p + start; 
                },
                this->parent->indices,
                this->parent->indptrs,
                [&](Index_ primary, Position& curptr) -> void {
                    update_secondary_position(curptr, this->parent->indices, this->parent->indptrs[primary + 1]);
                    output.add(primary, curptr.count);
                },
                [&](Index_ primary) -> void {
                    output.skip(primary);
                }
            );
            return;
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
                [&](Index_ primary, Position& curptr) -> void {
                    update_secondary_position(curptr, this->parent->indices, this->parent->indptrs[primary + 1]);
                    output.add(primary, curptr.count);
                },
                [&](Index_ primary) -> void {
                    output.skip(primary);
                }
            );
            return;
        }
    };

private:
    template<DimensionSelectionType selection_>
    struct DenseSecondaryExtractor : public SecondaryExtractorBase<selection_, false> {
        template<typename ...Args_>
        DenseSecondaryExtractor(const SemiCompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : SecondaryExtractorBase<selection_, false>(p, opt, std::forward<Args_>(args)...) {}

    private:
        struct ExpandedStoreIndexed {
            Value_* out_values;

            void add(Index_, Value_ val) {
                *out_values = val;
                ++out_values;
                return;
            }

            void skip(Index_) {
                ++out_values;
            }
        };

        struct ExpandedStoreBlock {
            Value_* out_values;
            Index_ first;

            static constexpr bool can_ignore_values = false;

            void add(Index_ i, Value_ val) {
                out_values[i - first] = val;
                return;
            }

            void skip(Index_) {} 
        };

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            typename std::conditional<selection_ == DimensionSelectionType::INDEX, ExpandedStoreIndexed, ExpandedStoreBlock>::type store;
            store.out_values = buffer;
            std::fill(buffer, buffer + extracted_length<selection_, Index_>(*this), static_cast<Value_>(0));

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                store.first = 0;
                this->secondary_dimension_loop(i, static_cast<Index_>(0), this->full_length, store);

            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                store.first = this->block_start;
                this->secondary_dimension_loop(i, this->block_start, this->block_length, store);

            } else {
                ExpandedStoreIndexed store;
                store.out_values = buffer;
                this->secondary_dimension_loop(i, this->subset_indices.data(), this->index_length, store);
            }

            return buffer;
        }
    };

private:
    template<DimensionSelectionType selection_>
    struct SparseSecondaryExtractor : public SecondaryExtractorBase<selection_, true> {
        template<typename ...Args_>
        SparseSecondaryExtractor(const SemiCompressedSparseMatrix* p, const Options& opt, Args_&& ... args) : SecondaryExtractorBase<selection_, true>(p, opt, std::forward<Args_>(args)...) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            if (!this->needs_value) {
                vbuffer = NULL;
            }
            if (!this->needs_index) {
                ibuffer = NULL;
            }

            RawStore store;
            store.out_values = vbuffer;
            store.out_indices = ibuffer;

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
 * Semi-compressed sparse column matrix.
 * See `tatami::SemiCompressedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class IndexStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
using SemiCompressedSparseColumnMatrix = SemiCompressedSparseMatrix<false, Value_, Index_, IndexStorage_, PointerStorage_>;

/**
 * Semi-compressed sparse row matrix.
 * See `tatami::SemiCompressedSparseMatrix` for details on the template parameters.
 */
template<typename Value_, typename Index_, class IndexStorage_ = std::vector<Index_>, class PointerStorage_ = std::vector<size_t> >
using SemiCompressedSparseRowMatrix = SemiCompressedSparseMatrix<true, Value_, Index_, IndexStorage_, PointerStorage_>;

}

#endif
