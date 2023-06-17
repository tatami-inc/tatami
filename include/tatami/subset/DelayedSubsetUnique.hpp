#ifndef TATAMI_DELAYED_SUBSET_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_UNIQUE_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetUnique.hpp
 *
 * @brief Delayed subsetting by unique row/column indices.
 */

namespace tatami {

/**
 * @brief Delayed subsetting of a matrix with unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of unique indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetUnique : public Matrix<Value_, Index_> {
private:
    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<IndexStorage_>()[0])>::type>::type storage_type;

    static void finish_assembly(
        const std::vector<std::pair<storage_type, Index_> >& collected,
        const IndexStorage_& indices, 
        std::vector<Index_>& reverse_mapping,
        std::vector<Index_>& sorted,
        Index_ mapping_dim,
        std::vector<Index_>& mapping_single
    ) {
        sorted.reserve(indices.size());
        reverse_mapping.resize(indices.size());

        for (Index_ i = 0, end = collected.size(); i < end; ++i) {
            const auto& current = collected[i];
            sorted.push_back(current.first);
            reverse_mapping[current.second] = sorted.size() - 1;
        }

        mapping_single.resize(mapping_dim);
        for (Index_  i = 0, end = indices.size(); i < end; ++i) {
            mapping_single[indices[i]] = i;
        }
    }

public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * This should be unique, but may be unsorted.
     * @param check Whether to check `idx` for unique values.
     */
    DelayedSubsetUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        std::vector<std::pair<storage_type, Index_> > collected;
        collected.reserve(indices.size());
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            collected.emplace_back(indices[i], i);
        }
        std::sort(collected.begin(), collected.end());

        if (check) {
            for (Index_ i = 0, end = collected.size(); i < end; ++i) {
                if (i) {
                    if (collected[i].first == collected[i-1].first) {
                        throw std::runtime_error("indices should be unique");
                        break;
                    }
                }
            }
        }

        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            sorted,
            margin_ == 0 ? mat->nrow() : mat->ncol(),
            mapping_single
        );
    }

    /**
     * @cond
     */
    DelayedSubsetUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, const std::vector<std::pair<storage_type, Index_> >& collected, IndexStorage_ idx) : 
        mat(std::move(p)), indices(std::move(idx)) 
    {
        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            sorted,
            margin_ == 0 ? mat->nrow() : mat->ncol(),
            mapping_single
        );
    }
    /**
     * @endcond
     */

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

    std::vector<Index_> reverse_mapping;
    std::vector<Index_> sorted;
    std::vector<Index_> mapping_single;

public:
    Index_ nrow() const {
        if constexpr(margin_==0) {
            return indices.size();
        } else {
            return mat->nrow();
        }
    }

    Index_ ncol() const {
        if constexpr(margin_==0) {
            return mat->ncol();
        } else {
            return indices.size();
        }
    }

    bool sparse() const {
        return mat->sparse();
    }

    double sparse_proportion() const {
        return mat->sparse_proportion();
    }

    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    double prefer_rows_proportion() const {
        return mat->prefer_rows_proportion();
    }

    bool uses_oracle(bool row) const {
        return mat->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::sparse_column;

    using Matrix<Value_, Index_>::sparse_row;

private:
    struct DenseBase {
        DenseBase(size_t bufsize) : temp(bufsize) {}
        std::vector<Value_> temp;
    };

    struct SparseBase {
        /*
         * There's no need to have a vbuffer here, because we can directly extract into the user-supplied vbuffer if we want the values. 
         * The inner extraction is guaranteed to be no greater than the vbuffer length, as it's a 1:1 mapping between unique indices in the row/column call and the unique internal subset indices.
         * Temporaries are not a concern because we end up copying everything into sortspace for re-sorting anyway.
         *
         * Similarly, an ibuffer is _usually_ unnecessary as we can write directly into the user-supplied ibuffer.
         * There's only a need to allocate our own ibuffer if (i) the user wants the values but not the indices, and (ii) they want the output to be sorted. 
         * In this case, the user-supplied ibuffer (i) may not be valid but (ii) we still need the indices, so we need to create our ibuffer.
         */ 

        SparseBase(const Options& opt, size_t bufsize) :
            report_index(opt.sparse_extract_index),
            needs_sort(opt.sparse_ordered_index),
            itemp(opt.sparse_extract_value && !opt.sparse_extract_index && needs_sort ? bufsize : 0) 
        {
            if (needs_sort) {
                sortspace.reserve(bufsize);
            }
        }

        bool report_index;
        bool needs_sort;
        std::vector<Index_> itemp;
        std::vector<std::pair<Index_, Value_> > sortspace;
    };

private:
    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > create_inner_extractor(const Options& opt, std::vector<Index_> idx) const {
        if constexpr(sparse_) {
            if (opt.sparse_ordered_index) {
                auto copy = opt;

                // Turning off the sorting to enable possible optimizations in the underlying matrix.
                // We don't need sorted output as we'll be resorting ourselves later.
                copy.sparse_ordered_index = false;

                if (!opt.sparse_extract_index && opt.sparse_extract_value) {
                    // Need to force extraction of indices for correct resorting.
                    copy.sparse_extract_index = true;
                }

                return new_extractor<margin_ != 0, sparse_>(mat.get(), std::move(idx), copy);
            }
        }

        return new_extractor<margin_ != 0, sparse_>(mat.get(), std::move(idx), opt);
    }

    template<class Extractor_>
    SparseRange<Value_, Index_> reorganize_sparse(
        Index_ i, 
        Value_* vbuffer, 
        Index_* ibuffer, 
        Extractor_* work,
        std::vector<Index_>& itemp, 
        std::vector<std::pair<Index_, Value_> >& sortspace, 
        bool report_index, 
        bool needs_sort) const 
    {
        if (!needs_sort) {
            // When we don't need sorting, validity of ibuffer and vbuffer should be consistent 
            // with the extraction mode used to construct 'work'.
            auto raw = work->fetch(i, vbuffer, ibuffer);
            if (raw.index) {
                auto icopy = ibuffer;
                for (Index_ i = 0; i < raw.number; ++i, ++icopy) {
                    *icopy = mapping_single[raw.index[i]];
                }
                raw.index = ibuffer;
            }
            return raw;
        }

        // If the workspace's ibuffer is empty, we're either extracting the indices
        // directly into the user's ibuffer, or we don't need the indices at all.
        // Either way, it doesn't hurt to use the user's ibuffer.
        Index_* iin = (itemp.empty() ? ibuffer : itemp.data()); 

        auto raw = work->fetch(i, vbuffer, iin);
        if (!raw.value && !raw.index) {
            return raw;
        }

        // raw.index is guaranteed to be valid as this point, as we forced the
        // extraction of indices if the user requests sorted output.
        sortspace.clear();
        if (raw.value) {
            for (Index_ i = 0; i < raw.number; ++i) {
                sortspace.emplace_back(mapping_single[raw.index[i]], raw.value[i]);
            }
        } else {
            for (Index_ i = 0; i < raw.number; ++i) {
                sortspace.emplace_back(mapping_single[raw.index[i]], 0); // no-op value.
            }
        }

        std::sort(sortspace.begin(), sortspace.end());

        if (raw.value) {
            auto vcopy = vbuffer;
            for (const auto& x : sortspace) {
                *vcopy = x.second;
                ++vcopy;
            }
        } else {
            vbuffer = NULL;
        }

        if (report_index) {
            auto icopy = ibuffer;
            for (const auto& x : sortspace) {
                *icopy = x.first;
                ++icopy;
            }
        } else {
            ibuffer = NULL;
        }

        return SparseRange<Value_, Index_>(raw.number, vbuffer, ibuffer);
    }

    /*****************************************
     ************ Full extraction ************
     *****************************************/
private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct ParallelExtractor : public Extractor<selection_, sparse_, Value_, Index_> {
        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            internal->set_oracle(std::move(o));
        }
    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
    };

    template<bool sparse_>
    struct FullParallelExtractor : public ParallelExtractor<DimensionSelectionType::FULL, sparse_> {
        FullParallelExtractor(const DelayedSubsetUnique* p, const Options& opt) : parent(p) {
            this->full_length = parent->indices.size();
            this->internal = parent->create_inner_extractor<sparse_>(opt, parent->sorted); // copy is deliberate.
        }
    protected:
        const DelayedSubsetUnique* parent;
    };

    struct FullDenseParallelExtractor : public FullParallelExtractor<false>, public DenseBase {
        FullDenseParallelExtractor(const DelayedSubsetUnique* p, const Options& opt) : 
            FullParallelExtractor<false>(p, opt), DenseBase(this->internal->index_length) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, this->parent->reverse_mapping);
        }
    };

    struct FullSparseParallelExtractor : public FullParallelExtractor<true>, public SparseBase {
        FullSparseParallelExtractor(const DelayedSubsetUnique* p, const Options& opt) : 
            FullParallelExtractor<true>(p, opt), SparseBase(opt, this->internal->index_length) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return this->parent->reorganize_sparse(
                i, vbuffer, ibuffer,
                this->internal.get(),
                this->itemp, this->sortspace, this->report_index, this->needs_sort
            );
        }
    };

    /******************************************
     ************ Block extraction ************
     ******************************************/
private:
    template<class Function_> 
    void transplant_indices(std::vector<Index_>& local, Index_ len, Function_ get_index) const {
        // Avoid re-sorting by going through the indices
        // and filtering the existing 'sorted'.
        local.resize(sorted.size());

        for (Index_ b = 0; b < len; ++b) {
            local[reverse_mapping[get_index(b)]] = 1;
        }

        Index_ counter = 0;
        for (Index_ b = 0, end = local.size(); b < end; ++b) {
            if (local[b]) {
                local[counter] = sorted[b];
                ++counter;
            }
        }

        local.resize(counter);
    }

    template<class Function>
    void transplant_indices(std::vector<Index_>& local, Index_ len, Function get_index, std::vector<Index_>& my_reverse_mapping) const {
        // Same principle, but we need to keep track of the indices of the hits
        // to create a reverse mapping vector for the dense extractors.
        std::vector<unsigned char> hits;
        hits.resize(sorted.size());
        local.resize(sorted.size());

        for (Index_ b = 0; b < len; ++b) {
            auto r = reverse_mapping[get_index(b)];
            hits[r] = 1;
            local[r] = b;
        }

        Index_ counter = 0;
        my_reverse_mapping.resize(len);
        for (Index_ b = 0, end = local.size(); b < end; ++b) {
            if (hits[b]) {
                my_reverse_mapping[local[b]] = counter;
                local[counter] = sorted[b];
                ++counter;
            }
        }

        local.resize(counter);
    }

private:
    template<bool sparse_>
    struct BlockParallelExtractor : public ParallelExtractor<DimensionSelectionType::BLOCK, sparse_> {
        BlockParallelExtractor(const DelayedSubsetUnique* p, const Options& opt, Index_ bs, Index_ bl) : parent(p) {
            this->block_start = bs;
            this->block_length = bl;

            std::vector<Index_> local;
            auto fun = [&](Index_ i) -> Index_ { return i  + bs; };
            if constexpr(sparse_) {
                parent->transplant_indices(local, bl, std::move(fun));
            } else {
                parent->transplant_indices(local, bl, std::move(fun), reverse_mapping);
            }

            this->internal = parent->create_inner_extractor<sparse_>(opt, std::move(local));
        }
    protected:
        const DelayedSubsetUnique* parent;
        typename std::conditional<!sparse_, std::vector<Index_>, bool>::type reverse_mapping;
    };

    struct BlockDenseParallelExtractor : public BlockParallelExtractor<false>, public DenseBase {
        BlockDenseParallelExtractor(const DelayedSubsetUnique* p, const Options& opt, Index_ bs, Index_ bl) : 
            BlockParallelExtractor<false>(p, opt, bs, bl), DenseBase(this->internal->index_length) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, this->reverse_mapping);
        }
    };

    struct BlockSparseParallelExtractor : public BlockParallelExtractor<true>, public SparseBase {
        BlockSparseParallelExtractor(const DelayedSubsetUnique* p, const Options& opt, Index_ bs, Index_ bl) : 
            BlockParallelExtractor<true>(p, opt, bs, bl), SparseBase(opt, this->internal->index_length) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return this->parent->reorganize_sparse(
                i, vbuffer, ibuffer,
                this->internal.get(),
                this->itemp, this->sortspace, this->report_index, this->needs_sort
            );
        }
    };

    /*****************************************
     *********** Index extraction ************
     *****************************************/
private:
    template<bool sparse_>
    struct IndexParallelExtractor : public ParallelExtractor<DimensionSelectionType::INDEX, sparse_> {
        IndexParallelExtractor(const DelayedSubsetUnique* p, const Options& opt, std::vector<Index_> idx) : parent(p) {
            Index_ il = idx.size();
            this->index_length = il;
            indices = std::move(idx);

            std::vector<Index_> local;
            auto fun = [&](Index_ i) -> Index_ { return indices[i]; };
            if constexpr(sparse_) {
                parent->transplant_indices(local, il, std::move(fun));
            } else {
                parent->transplant_indices(local, il, std::move(fun), reverse_mapping);
            }

            this->internal = parent->create_inner_extractor<sparse_>(opt, std::move(local));
        }

        const Index_* index_start() const {
            return indices.data();
        }

    protected:
        const DelayedSubsetUnique* parent;
        std::vector<Index_> indices;
        typename std::conditional<!sparse_, std::vector<Index_>, bool>::type reverse_mapping;
    };

    struct IndexDenseParallelExtractor : public IndexParallelExtractor<false>, public DenseBase {
        IndexDenseParallelExtractor(const DelayedSubsetUnique* p, const Options& opt, std::vector<Index_> idx) : 
            IndexParallelExtractor<false>(p, opt, std::move(idx)), DenseBase(this->internal->index_length) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, this->reverse_mapping);
        }
    };

    struct IndexSparseParallelExtractor : public IndexParallelExtractor<true>, public SparseBase {
        IndexSparseParallelExtractor(const DelayedSubsetUnique* p, const Options& opt, std::vector<Index_> idx) : 
            IndexParallelExtractor<true>(p, opt, std::move(idx)), SparseBase(opt, this->internal->index_length) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return this->parent->reorganize_sparse(
                i, vbuffer, ibuffer,
                this->internal.get(),
                this->itemp, this->sortspace, this->report_index, this->needs_sort
            );
        }
    };

    /**************************************************
     ************ Public virtual overrides ************
     **************************************************/
private:
    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::FULL, sparse_, Value_, Index_> > populate_parallel(const Options& options) const {
        std::unique_ptr<Extractor<DimensionSelectionType::FULL, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new FullSparseParallelExtractor(this, options));
        } else {
            output.reset(new FullDenseParallelExtractor(this, options));
        }
        return output;
    }

    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > populate_parallel(const Options& options, Index_ bs, Index_ bl) const {
        std::unique_ptr<Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new BlockSparseParallelExtractor(this, options, bs, bl));
        } else {
            output.reset(new BlockDenseParallelExtractor(this, options, bs, bl));
        }
        return output;
    }

    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > populate_parallel(const Options& options, std::vector<Index_> idx) const {
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new IndexSparseParallelExtractor(this, options, std::move(idx)));
        } else {
            output.reset(new IndexDenseParallelExtractor(this, options, std::move(idx)));
        }
        return output;
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options& options, Args_&& ... args) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            return subset_utils::populate_perpendicular<accrow_, selection_, sparse_>(mat.get(), indices, options, std::forward<Args_>(args)...);
        } else {
            return populate_parallel<sparse_>(options, std::forward<Args_>(args)...);
        }
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& options) const {
        return populate<true, DimensionSelectionType::FULL, false>(options);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& options) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(options, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& options) const {
        return populate<true, DimensionSelectionType::INDEX, false>(options, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& options) const {
        return populate<false, DimensionSelectionType::FULL, false>(options);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& options) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(options, block_start, block_length);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& options) const {
        return populate<false, DimensionSelectionType::INDEX, false>(options, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& options) const {
        return populate<true, DimensionSelectionType::FULL, true>(options);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& options) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(options, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& options) const {
        return populate<true, DimensionSelectionType::INDEX, true>(options, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& options) const {
        return populate<false, DimensionSelectionType::FULL, true>(options);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& options) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(options, block_start, block_length);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& options) const {
        return populate<false, DimensionSelectionType::INDEX, true>(options, std::move(indices));
    }
};

}

#endif
