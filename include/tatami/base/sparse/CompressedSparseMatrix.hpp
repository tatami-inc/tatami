#ifndef TATAMI_COMPRESSED_SPARSE_MATRIX_H
#define TATAMI_COMPRESSED_SPARSE_MATRIX_H

#include "../Matrix.hpp"
#include "../utils.hpp"

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
 * If `false`, a compressed sparse column representation is used instead.
 * @tparam Value_ Type of the matrix values.
 * @tparam Index_ Type of the row/column indices.
 * @tparam ValueStorage_ Vector class used to store the matrix values internally.
 * This does not necessarily have to contain `T`, as long as the type is convertible to `T`.
 * Methods should be available for `size()`, `begin()`, `end()` and `[]`.
 * If a method is available for `data()` that returns a `const T*`, it will also be used.
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
     * @param idx Vector of row indices (if `ROW=false`) or column indices (if `ROW=true`) for the non-zero elements.
     * @param ptr Vector of index pointers.
     * @param check Should the input vectors be checked for validity?
     *
     * If `check=true`, the constructor will check that `vals` and `idx` have the same length;
     * `ptr` is ordered with first and last values set to 0 and the number of non-zero elements, respectively;
     * and `idx` is ordered within each interval defined by successive elements of `ptr`.
     */
    CompressedSparseMatrix(Index_ nr, Index_ nc, ValueStorage_ vals, IndexStorage_ idx, PointerStorage_ ptr, bool check=true) : 
        nrows(nr), ncols(nc), values(std::move(vals)), indices(std::move(idx)), indptrs(std::move(ptr)) 
    {
        check_values(check); 
        return;
    }

private:
    Index_ nrows, ncols;
    ValueStorage_ values;
    IndexStorage_ indices;
    PointerStorage_ indptrs;

    void check_values(bool check) {
        if (!check) {
            return;
        }

        if (values.size() != indices.size()) {
            throw std::runtime_error("'values' and 'indices' should be of the same length");
        }

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
        if (static_cast<size_t>(indptrs[indptrs.size() - 1]) != indices.size()) {
            throw std::runtime_error("last element of 'indptrs' should be equal to length of 'indices'");
        }

        size_t counter = 0;
        for (size_t i = 1; i < indptrs.size(); ++i) {
            if (indptrs[i] < indptrs[i-1]) {
                throw std::runtime_error("'indptrs' should be in increasing order");
            }

            if (counter < indices.size()) {
                auto previous = indices[counter];
                ++counter;
                while (counter < static_cast<size_t>(indptrs[i])) {
                    if (previous >= indices[counter]) {
                        if (row_) {
                            throw std::runtime_error("'indices' should be strictly increasing within each row");
                        } else {
                            throw std::runtime_error("'indices' should be strictly increasing within each column");
                        }
                    }
                    ++counter;
                }
            }
        }

        return;
    }

public:
    Index_ nrow() const { return nrows; }

    Index_ ncol() const { return ncols; }

    bool sparse() const { return true; }

    /**
     * @return `true` if `ROW = true` (for `CompressedSparseRowMatrix` objects), otherwise returns `false` (for `CompressedSparseColumnMatrix` objects).
     */
    bool prefer_rows() const { return row_; }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_, bool sparse_>
    struct CompressedExtractorBase : public StandardExtractor<accrow_, sparse_, Value_, Index_> {
        CompressedExtractorBase(const CompressedSparseMatrix* p, ExtractionOptions<Index_>& opt) :
            StandardExtractor<accrow_, sparse_, Value_, Index_>(p, opt), // (note: potentially invalidates non_target.selections.indices)
            parent(p) 
        {}

    protected:
        const CompressedSparseMatrix* parent;
    };

    /***********************************
     ******* Primary extraction ********
     ***********************************/
private:
    struct PrimaryWorkspace {
        std::vector<std::pair<size_t, size_t> > cached;
    };

    Index_ max_secondary_index() const {
        if constexpr(row_) {
            return ncols;
        } else {
            return nrows;
        }
    }

    std::pair<size_t, size_t> primary_dimension(Index_ i, Index_ start, Index_ length, PrimaryWorkspace& work) const {
        bool do_cache = !(work.cached.empty());
        if (do_cache) {
            auto val = work.cached[i];
            if (val.first != -1) {
                return val;
            }
        }

        auto iIt = indices.begin() + indptrs[i], eIt = indices.begin() + indptrs[i + 1]; 

        if (start) { // Jumping ahead if non-zero.
            iIt = std::lower_bound(iIt, eIt, start);
        } 

        auto last = start + length;
        if (last != max_secondary_index()) { // Jumping to last element.
            eIt = std::lower_bound(iIt, eIt, last);
        }

        size_t outstart = iIt - indices.begin();
        size_t outlength = eIt - iIt;
        if (do_cache) {
            work.cached[i].first = outstart;
            work.cached[i].second = outlength;
        }

        return std::make_pair(outstart, outlength);
    }

    SparseRange<Value_, Index_> primary_dimension_raw(Index_ i, Index_ start, Index_ length, PrimaryWorkspace& work, Value_* out_values, Index_* out_indices) const {
        auto obtained = primary_dimension(i, start, length, work);
        SparseRange<Value_, Index_> output(obtained.second);

        if (out_values) {
            if constexpr(has_data<Value_, ValueStorage_>::value) {
                output.value = values.data() + obtained.first;
            } else {
                auto vIt = values.begin() + obtained.first;
                std::copy(vIt, vIt + obtained.second, out_values);
                output.value = out_values;
            }
        } else {
            output.value = NULL;
        }

        if (out_indices) {
            if constexpr(has_data<Index_, IndexStorage_>::value) {
                output.index = indices.data() + obtained.first;
            } else {
                auto iIt = indices.begin() + obtained.first;
                std::copy(iIt, iIt + obtained.second, out_indices);
                output.index = out_indices;
            }
        } else {
            output.index = NULL;
        }

        return output;
    }

    void primary_dimension_expanded(Index_ i, Index_ start, Index_ length, PrimaryWorkspace& work, Value_* out_values) const {
        std::fill(out_values, out_values + length, static_cast<Value_>(0));
        auto obtained = primary_dimension(i, start, length, work);
        auto vIt = values.begin() + obtained.first;
        auto iIt = indices.begin() + obtained.first;
        for (size_t x = 0; x < obtained.second; ++x, ++vIt, ++iIt) {
            out_values[*iIt - start] = *vIt;
        }
        return;
    }

private:
    template<class Store_>
    void add_to_store(Index_ primary, size_t curptr, Store_& output) const {
        // If the value isn't actually needed, we don't want to do a look-up on 'values'.
        if constexpr(Store_::can_ignore_values) {
            if (output.ignore_values()) {
                output.add(primary);
                return;
            }
        }
        output.add(primary, values[curptr]);
    }

    template<class Store_>
    void primary_dimension(Index_ i, const Index_* subset, Index_ length, PrimaryWorkspace& work, Store_& store) const {
        if (!length) {
            return;
        }

        auto iIt = indices.begin() + indptrs[i], eIt = indices.begin() + indptrs[i + 1]; 

        if (indices[0]) { // Only jumping ahead if the start is non-zero.
            bool do_cache = !work.cached.empty();
            if (do_cache) {
                if (work.cached[i].first != -1) { // retrieving the jump from cache, if we came here before.
                    iIt += work.cached[i].first;
                } else {
                    auto iIt2 = std::lower_bound(iIt, eIt, subset[0]);
                    work.cached[i].first = iIt2 - iIt;
                    iIt = iIt2;
                }
            } else {
                iIt = std::lower_bound(iIt, eIt, subset[0]);
            }
        } 

        if (iIt == eIt) {
            return;
        }

        Index_ counter = 0;
        while (counter < length) {
            auto current = subset[counter];

            while (iIt != eIt && current > *iIt) {
                ++iIt;
            }
            if (iIt == eIt) {
                break;
            }

            if (current == *iIt) {
                add_to_store(current, iIt - indices.begin(), store);
            } else {
                store.skip(current);
            }
            ++counter;
        }

        return;
    }

    // Overview of the Store contract:
    struct RawStore {
        Value_* out_values;
        Index_* out_indices;
        Index_ n = 0;

        // Whether or not 'values' might be ignored, typically
        // because a NULL pointer is supplied in the 'vbuffer'
        static constexpr bool can_ignore_values = true;

        // If 'can_ignore_values = true', we need this method
        // to define whether the values are ignored for this instance.
        bool ignore_values() const {
            return !out_values;
        }

        // If 'can_ignore_values = true', we need this method to handle
        // hits, but _without_ a costly fetch of the value from memory.
        void add(Index_ i) {
            ++n;
            if (out_indices) {
                *out_indices = i;
                ++out_indices;
            }
        }

        // Storing a hit in terms of its index and its value. If 'can_ignore_values = true',
        // it can be assumed that 'skip_values()' is false if this method is called. Note 
        // that add() and skip() is always called on indices in increasing order.
        void add(Index_ i, Value_ val) {
            add(i);
            *out_values = val;
            ++out_values;
            return;
        }

        // Skip an index if it isn't found in the data. 
        void skip(Index_) {} 
    };

    SparseRange<Value_, Index_> primary_dimension_raw(Index_ i, const Index_* indices, Index_ length, PrimaryWorkspace& work, Value_* out_values, Index_* out_indices) const {
        RawStore store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        primary_dimension(i, indices, length, work, store);
        return SparseRange<Value_, Index_>(store.n, out_values, out_indices);
    } 

    struct ExpandedStoreIndexed {
        Value_* out_values;

        static constexpr bool can_ignore_values = false;

        void add(Index_, Value_ val) {
            *out_values = val;
            ++out_values;
            return;
        }

        void skip(Index_) {
            ++out_values;
        }
    };

    void primary_dimension_expanded(Index_ i, const Index_* indices, Index_ length, PrimaryWorkspace& work, Value_* out_values) const {
        std::fill(out_values, out_values + length, static_cast<Value_>(0));
        ExpandedStoreIndexed store;
        store.out_values = out_values;
        primary_dimension(i, indices, length, work, store);
        return;
    }

private:
    template<bool sparse_>
    struct PrimaryExtractorBase : public CompressedExtractorBase<row_, sparse_> {
        PrimaryExtractorBase(const CompressedSparseMatrix* p, IterationOptions<Index_>& target, ExtractionOptions<Index_>& non_target) : 
            CompressedExtractorBase<row_, sparse_>(p, non_target) 
        {
            bool spawn_cache = false;

            switch (non_target.selection.type) {
                case DimensionSelectionType::BLOCK:
                    // only need to create a cache if block does not start at 0, see primary_dimension for details.
                    spawn_cache = (target.cache_for_reuse && this->extracted_block); 
                    break;
                case DimensionSelectionType::INDEX:
                    // only need to create a cache if indices are non-empty and the first is not 0, see primary_dimension for details.
                    spawn_cache = (target.cache_for_reuse && this->extracted_length && this->quick_extracted_index()[0]);
                    break;
            }

            if (spawn_cache) {
                work.cached.resize(row_ ? this->parent->nrows : this->parent->ncols, std::pair<size_t, size_t>(-1, 0));
            }
        }

    protected:
        PrimaryWorkspace work;
    };

    struct DensePrimaryExtractor : public PrimaryExtractorBase<false> {
        DensePrimaryExtractor(const CompressedSparseMatrix* p, IterationOptions<Index_>& target, ExtractionOptions<Index_>& non_target) : 
            PrimaryExtractorBase<false>(p, target, non_target) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            switch (this->extracted_selection) {
                case DimensionSelectionType::FULL:
                    this->parent->primary_dimension_expanded(i, 0, this->extracted_length, this->work, buffer);
                    break;
                case DimensionSelectionType::BLOCK:
                    this->parent->primary_dimension_expanded(i, this->extracted_block, this->extracted_length, this->work, buffer);
                    break;
                case DimensionSelectionType::INDEX:
                    this->parent->primary_dimension_expanded(i, this->quick_extracted_index(), this->extracted_length, this->work, buffer);
                    break;
            }
            return buffer;
        }
    };

    struct SparsePrimaryExtractor : public PrimaryExtractorBase<true> {
        SparsePrimaryExtractor(const CompressedSparseMatrix* p, IterationOptions<Index_>& target, ExtractionOptions<Index_>& non_target) : 
            PrimaryExtractorBase<true>(p, target, non_target) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            switch (this->extracted_selection) {
                case DimensionSelectionType::FULL:
                    return this->parent->primary_dimension_raw(i, 0, this->extracted_length, this->work, vbuffer, ibuffer);
                case DimensionSelectionType::BLOCK:
                    return this->parent->primary_dimension_raw(i, this->extracted_block, this->extracted_length, this->work, vbuffer, ibuffer);
                case DimensionSelectionType::INDEX:
                    return this->parent->primary_dimension_raw(i, this->quick_extracted_index(), this->extracted_length, this->work, vbuffer, ibuffer);
            }
        }
    };

    /*************************************
     ******* Secondary extraction ********
     *************************************/
private:
    struct SecondaryWorkspace {
        typedef typename std::remove_reference<decltype(std::declval<IndexStorage_>()[0])>::type index_type;
        typedef typename std::remove_reference<decltype(std::declval<PointerStorage_>()[0])>::type indptr_type;

        SecondaryWorkspace() = default;

        SecondaryWorkspace(Index_ max_index, const IndexStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
            current_indptrs(idp.begin() + start, idp.begin() + start + length), current_indices(length)
        {
            /* Here, the general idea is to store a local copy of the actual
             * row indices (for CSC matrices; column indices, for CSR matrices)
             * so that we don't have to keep on doing cache-unfriendly look-ups
             * for the indices based on the pointers that we do have. This assumes
             * that the density is so low that updates to the local indices are
             * rare relative to the number of comparisons to those same indices.
             * Check out the `secondary_dimension()` function for how this is used.
             */
            auto idpIt = idp.begin() + start;
            for (Index_ i = 0; i < length; ++i, ++idpIt) {
                current_indices[i] = (*idpIt < *(idpIt + 1) ? idx[*idpIt] : max_index);
            }
            return;
        } 

        SecondaryWorkspace(Index_ max_index, const IndexStorage_& idx, const PointerStorage_& idp) :
            SecondaryWorkspace(max_index, idx, idp, 0, idp.size() - 1) {}

        SecondaryWorkspace(Index_ max_index, const IndexStorage_& idx, const PointerStorage_& idp, const Index_* subset, Index_ length) :
            current_indptrs(length), current_indices(length)
        {
            for (Index_ i0 = 0; i0 < length; ++i0) {
                auto i = subset[i0];
                current_indptrs[i0] = idp[i];
                current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
            }
            return;
        } 

        Index_ previous_request = 0;
        std::vector<indptr_type> current_indptrs; // the current position of the pointer

        // The current index being pointed to, i.e., current_indices[0] <= indices[current_ptrs[0]]. 
        // If current_ptrs[0] is out of range, the current index is instead set to the maximum index
        // (e.g., max rows for CSC matrices).
        std::vector<index_type> current_indices;
    };

private:
    template<class Store_>
    void secondary_dimension_above(Index_ secondary, Index_ primary, Index_ index_primary, SecondaryWorkspace& work, Store_& output) const {
        auto& curdex = work.current_indices[index_primary];
        auto& curptr = work.current_indptrs[index_primary];

        // Remember that we already store the index corresponding to the current indptr.
        // So, we only need to do more work if the request is greater than that index.
        if (secondary > curdex) {
            auto max_index = max_secondary_index();
            auto limit = indptrs[primary + 1];

            // Having a peek at the index of the next non-zero
            // element; maybe we're lucky enough that the requested
            // index is below this, as would be the case for
            // consecutive or near-consecutive accesses.
            ++curptr;
            if (curptr < limit) {
                auto candidate = indices[curptr];
                if (candidate >= secondary) {
                    curdex = candidate;

                } else if (secondary + 1 == max_index) {
                    // Special case if the requested index is at the end of the matrix,
                    // in which case we can just jump there directly rather than 
                    // doing an unnecessary binary search.
                    curptr = limit - 1;
                    if (indices[curptr] == secondary) {
                        curdex = indices[curptr];
                    } else {
                        ++curptr;
                        curdex = max_index;
                    }

                } else {
                    // Otherwise we need to search indices above the existing position.
                    curptr = std::lower_bound(indices.begin() + curptr, indices.begin() + limit, secondary) - indices.begin();
                    curdex = (curptr < limit ? indices[curptr] : max_index);
                }
            } else {
                curdex = max_index;
            }
        } 

        if (secondary == curdex) { // assuming secondary < max_index, of course.
            add_to_store(primary, curptr, output);
        } else {
            output.skip(primary);
        }
    }

    template<class Store_>
    void secondary_dimension_below(Index_ secondary, Index_ primary, Index_ index_primary, SecondaryWorkspace& work, Store_& output) const {
        auto& curdex = work.current_indices[index_primary];
        auto& curptr = work.current_indptrs[index_primary];

        if (secondary < curdex) {
            auto limit = indptrs[primary + 1];
            if (indptrs[primary] < limit) {
                if (secondary == 0) {
                    // Special case if the requested index is at the end of the matrix,
                    // in which case we can just jump there directly rather than 
                    // doing an unnecessary binary search.
                    curptr = indptrs[primary];
                    curdex = indices[curptr];
                } else {
                    // Otherwise, searching indices below the existing position.
                    curptr = std::lower_bound(indices.begin() + indptrs[primary], indices.begin() + curptr, secondary) - indices.begin();
                    curdex = (curptr < limit ? indices[curptr] : max_secondary_index());
                }
            }
        }

        if (secondary == curdex) { // assuming secondary < max_index, of course.
            add_to_store(primary, curptr, output);
        } else {
            output.skip(primary);
        }
    }

private:
    template<class Store_>
    void secondary_dimension_loop(Index_ i, Index_ start, Index_ length, SecondaryWorkspace& work, Store_& store) const {
        Index_ prev_i = work.previous_request;

        if (i >= prev_i) {
            for (Index_ current = 0; current < length; ++current) {
                secondary_dimension_above(i, current + start, current, work, store);
            }
        } else {
            for (Index_ current = 0; current < length; ++current) {
                secondary_dimension_below(i, current + start, current, work, store);
            }
        }

        work.previous_request = i;
        return;
    }

    SparseRange<Value_, Index_> secondary_dimension_raw(Index_ i, Index_ start, Index_ length, SecondaryWorkspace& work, Value_* out_values, Index_* out_indices) const {
        RawStore store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        secondary_dimension_loop(i, start, length, work, store);
        return SparseRange<Value_, Index_>(store.n, out_values, out_indices);
    }

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

    void secondary_dimension_expanded(Index_ i, Index_ start, Index_ length, SecondaryWorkspace& work, Value_* out_values) const {
        std::fill(out_values, out_values + length, static_cast<Value_>(0));
        ExpandedStoreBlock store;
        store.out_values = out_values;
        store.first = start;
        secondary_dimension_loop(i, start, length, work, store);
        return;
    }

    template<class Store_>
    void secondary_dimension_loop(Index_ i, const Index_* subset, Index_ length, SecondaryWorkspace& work, Store_& output) const {
        Index_ prev_i = work.previous_request;

        if (i >= prev_i) {
            for (Index_ current = 0; current < length; ++current) {
                secondary_dimension_above(i, subset[current], current, work, output);
            }
        } else {
            for (Index_ current = 0; current < length; ++current) {
                secondary_dimension_below(i, subset[current], current, work, output);
            }
        }

        work.previous_request = i;
        return;
    }

    SparseRange<Value_, Index_> secondary_dimension_raw(Index_ i, const Index_* subset, Index_ length, SecondaryWorkspace& work, Value_* out_values, Index_* out_indices) const {
        RawStore store;
        store.out_values = out_values;
        store.out_indices = out_indices;
        secondary_dimension_loop(i, subset, length, work, store);
        return SparseRange<Value_, Index_>(store.n, out_values, out_indices);
    }

    void secondary_dimension_expanded(Index_ i, const Index_* subset, Index_ length, SecondaryWorkspace& work, Value_* out_values) const {
        std::fill(out_values, out_values + length, static_cast<Value_>(0));
        ExpandedStoreIndexed store;
        store.out_values = out_values;
        secondary_dimension_loop(i, subset, length, work, store);
        return;
    }

private:
    template<bool sparse_>
    struct SecondaryExtractorBase : public CompressedExtractorBase<!row_, sparse_> {
        SecondaryExtractorBase(const CompressedSparseMatrix* p, ExtractionOptions<Index_>& non_target) : CompressedExtractorBase<!row_, sparse_>(p, non_target) {
            auto max_index = row_ ? this->parent->nrows : this->parent->ncols; // maximum possible index of the primary dimension, as we're iterating over the secondary dimension.

            switch (non_target.selection.type) {
                case DimensionSelectionType::FULL:
                    work = SecondaryWorkspace(max_index, this->parent->indices, this->parent->indptrs);
                    break;
                case DimensionSelectionType::BLOCK:
                    work = SecondaryWorkspace(max_index, this->parent->indices, this->parent->indptrs, this->extracted_block, this->extracted_length);
                    break;
                case DimensionSelectionType::INDEX:
                    work = SecondaryWorkspace(max_index, this->parent->indices, this->parent->indptrs, this->quick_extracted_index(), this->extracted_length);
                    break;
            }
        }

        // Keep public for testing.
        SecondaryWorkspace work;
    };

    struct DenseSecondaryExtractor : public SecondaryExtractorBase<false> {
        DenseSecondaryExtractor(const CompressedSparseMatrix* p, ExtractionOptions<Index_>& non_target) : SecondaryExtractorBase<false>(p, non_target) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            switch (this->extracted_selection) {
                case DimensionSelectionType::FULL:
                    this->parent->secondary_dimension_expanded(i, 0, this->extracted_length, this->work, buffer);
                    break;
                case DimensionSelectionType::BLOCK:
                    this->parent->secondary_dimension_expanded(i, this->extracted_block, this->extracted_length, this->work, buffer);
                    break;
                case DimensionSelectionType::INDEX:
                    this->parent->secondary_dimension_expanded(i, this->quick_extracted_index(), this->extracted_length, this->work, buffer);
                    break;
            }
            return buffer;
        }
    };

    struct SparseSecondaryExtractor : public SecondaryExtractorBase<true> {
        SparseSecondaryExtractor(const CompressedSparseMatrix* p, ExtractionOptions<Index_>& non_target) : SecondaryExtractorBase<true>(p, non_target) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            switch (this->extracted_selection) {
                case DimensionSelectionType::FULL:
                    return this->parent->secondary_dimension_raw(i, 0, this->extracted_length, this->work, vbuffer, ibuffer);
                case DimensionSelectionType::BLOCK:
                    return this->parent->secondary_dimension_raw(i, this->extracted_block, this->extracted_length, this->work, vbuffer, ibuffer);
                case DimensionSelectionType::INDEX:
                    return this->parent->secondary_dimension_raw(i, this->quick_extracted_index(), this->extracted_length, this->work, vbuffer, ibuffer);
            }
        }
    };

    /*************************************
     ******* Extraction overrides ********
     *************************************/
public:
    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_row(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        if constexpr(row_) {
            return std::unique_ptr<DenseExtractor<Value_, Index_> >(new DensePrimaryExtractor(this, iopt, eopt));
        } else {
            return std::unique_ptr<DenseExtractor<Value_, Index_> >(new DenseSecondaryExtractor(this, eopt));
        }
    }

    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_column(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        if constexpr(row_) {
            return std::unique_ptr<DenseExtractor<Value_, Index_> >(new DenseSecondaryExtractor(this, eopt));
        } else {
            return std::unique_ptr<DenseExtractor<Value_, Index_> >(new DensePrimaryExtractor(this, iopt, eopt));
        }
    }

    std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_row(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        if constexpr(row_) {
            return std::unique_ptr<SparseExtractor<Value_, Index_> >(new SparsePrimaryExtractor(this, iopt, eopt));
        } else {
            return std::unique_ptr<SparseExtractor<Value_, Index_> >(new SparseSecondaryExtractor(this, eopt));
        }
    }

    std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_column(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        if constexpr(row_) {
            return std::unique_ptr<SparseExtractor<Value_, Index_> >(new SparseSecondaryExtractor(this, eopt));
        } else {
            return std::unique_ptr<SparseExtractor<Value_, Index_> >(new SparsePrimaryExtractor(this, iopt, eopt));
        }
    }

public:
    // For testing only.
    typedef SparseSecondaryExtractor SSE;
    typedef SecondaryWorkspace SW;
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
