#ifndef TATAMI_DELAYED_SUBSET_HPP
#define TATAMI_DELAYED_SUBSET_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubset.hpp
 *
 * @brief Delayed subsetting by rows or columns.
 *
 * This is equivalent to the `DelayedSubset` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed subsetting of a matrix with general indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of arbitrary indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubset : public Matrix<Value_, Index_> {
private:
    typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<IndexStorage_>()[0])>::type>::type storage_type;

    static void finish_assembly(
        const std::vector<std::pair<storage_type, Index_> >& collected,
        const IndexStorage_& indices, 
        std::vector<Index_>& reverse_mapping,
        std::vector<Index_>& unique_and_sorted,
        Index_ mapping_dim,
        std::vector<std::pair<Index_, Index_> >& mapping_duplicates,
        std::vector<Index_>& mapping_duplicates_pool
    ) {
        unique_and_sorted.reserve(indices.size());
        reverse_mapping.resize(indices.size());

        mapping_duplicates.resize(mapping_dim);
        mapping_duplicates_pool.reserve(indices.size());

        for (Index_ i = 0, end = collected.size(); i < end; ++i) {
            const auto& current = collected[i];
            auto& range = mapping_duplicates[current.first];
            if (unique_and_sorted.empty() || current.first != unique_and_sorted.back()) {
                unique_and_sorted.push_back(current.first);
                range.first = mapping_duplicates_pool.size();
            }

            mapping_duplicates_pool.push_back(current.second);
            reverse_mapping[current.second] = unique_and_sorted.size() - 1;
            ++range.second;
        }
    }

public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * These may be duplicated and/or unsorted.
     */
    DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx) : mat(std::move(p)), indices(std::move(idx)) {
        std::vector<std::pair<storage_type, Index_> > collected;
        collected.reserve(indices.size());
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            collected.emplace_back(indices[i], i);
        }
        std::sort(collected.begin(), collected.end());

        finish_assembly(
            collected,
            indices, 
            reverse_mapping,
            unique_and_sorted,
            margin_ == 0 ? mat->nrow() : mat->ncol(),
            mapping_duplicates,
            mapping_duplicates_pool
        );

        return;
    }

    /**
     * @cond
     */
    DelayedSubset(std::shared_ptr<const Matrix<Value_, Index_> > p, const std::vector<std::pair<storage_type, Index_> >& collected, IndexStorage_ idx) : 
        mat(std::move(p)), indices(std::move(idx)) 
    {
        finish_assembly(
            collected,
            indices, 
            reverse_mapping, 
            unique_and_sorted,
            margin_ == 0 ? mat->nrow() : mat->ncol(),
            mapping_duplicates,
            mapping_duplicates_pool
        );
    }
    /**
     * @endcond
     */

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

    std::vector<Index_> reverse_mapping;
    std::vector<Index_> unique_and_sorted;
    std::vector<std::pair<Index_, Index_> > mapping_duplicates; // holds (position, size) in the pool.
    std::vector<Index_> mapping_duplicates_pool; 

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
         * The behavior of the temp buffers depends on whether we want sorted output.
         *
         * If we want sorted output, then the logic roughly follows that of DelayedSubsetUnique::SparseBase.
         * We can directly extract into the user-supplied ibuffer/vbuffer and then copy everything into sortspace for resorting. 
         * The only exception is if the doesn't want the indices, in which case we need to do some allocation.
         *
         * If we don't want sorted output, then the logic follows that of DelayedSubsetSorted::SparseBase.
         * Here, we can't extract into the user-supplied ibuffer/vbuffer, because we need some
         * place to hold the to-be-processed values while we expand the duplicates.
         * Of course, we can omit the vtemp allocation if the user doesn't want the values.
         */ 

        SparseBase(const Options& opt, size_t bufsize) : report_index(opt.sparse_extract_index), needs_sort(opt.sparse_ordered_index) {
            if (needs_sort) {
                if (!opt.sparse_extract_index) {
                    itemp.resize(bufsize);
                }
                sortspace.reserve(bufsize);
            } else {
                if (opt.sparse_extract_value) {
                    vtemp.resize(bufsize);
                }
                itemp.resize(bufsize);
            }
        }

    protected:
        bool report_index;
        bool needs_sort;
        std::vector<Value_> vtemp;
        std::vector<Index_> itemp;
        std::vector<std::pair<Index_, Value_> > sortspace;
    };

private:
    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > create_inner_extractor(const Options& opt, std::vector<Index_> idx) const {
        if constexpr(sparse_) {
            if (opt.sparse_ordered_index || !opt.sparse_extract_index) {
                auto copy = opt;

                // Turning off the sorting to enable possible optimizations in the underlying matrix.
                // We don't need sorted output as we'll be resorting ourselves later.
                copy.sparse_ordered_index = false; 

                // Need to force extraction of indices for correct duplication.
                copy.sparse_extract_index = true;

                return new_extractor<margin_ != 0, sparse_>(mat.get(), std::move(idx), copy);
            }
        }

        return new_extractor<margin_ != 0, sparse_>(mat.get(), std::move(idx), opt);
    }

    template<class Extractor_>
    static SparseRange<Value_, Index_> reorganize_sparse_unsorted(
        Index_ i, 
        Value_* vbuffer, 
        Index_* ibuffer, 
        std::vector<Value_>& vtemp,
        std::vector<Index_>& itemp,
        bool report_index,
        Extractor_* work, 
        const std::vector<std::pair<Index_, Index_> >& dups, 
        const std::vector<Index_>& pool)
    {
        // Allocation status of vtemp depends on the extraction mode used to construct 'work'.
        Value_* vin = vtemp.data();

        // itemp should always be allocated, as we need this to get the expanded counts.
        Index_* iin = itemp.data();

        auto raw = work->fetch(i, vin, iin);
        if (!raw.value) {
            vbuffer = NULL;
        }
        if (!report_index) {
            ibuffer = NULL;
        }

        auto vcopy = vbuffer;
        auto icopy = ibuffer;
        Index_ counter = 0;

        for (Index_ i = 0; i < raw.number; ++i) {
            const auto& pool_pos = dups[raw.index[i]];
            counter += pool_pos.second;

            if (vcopy) {
                std::fill(vcopy, vcopy + pool_pos.second, raw.value[i]);
                vcopy += pool_pos.second;
            }

            if (icopy) {
                auto istart = pool.begin() + pool_pos.first;
                std::copy(istart, istart + pool_pos.second, icopy);
                icopy += pool_pos.second;
            }
        }

        return SparseRange<Value_, Index_>(counter, vbuffer, ibuffer);
    }

    template<class Extractor_>
    static SparseRange<Value_, Index_> reorganize_sparse_sorted(
        Index_ i, 
        Value_* vbuffer, 
        Index_* ibuffer, 
        std::vector<std::pair<Index_, Value_> >& sortspace,
        std::vector<Index_>& itemp,
        bool report_index,
        Extractor_* work, 
        const std::vector<std::pair<Index_, Index_> >& dups, 
        const std::vector<Index_>& pool)
    {
        // itemp may not be allocated, if indices are to be extracted into ibuffer.
        Index_* iin = (itemp.empty() ? ibuffer : itemp.data());
        auto raw = work->fetch(i, vbuffer, iin);

        sortspace.clear();
        for (Index_ i = 0; i < raw.number; ++i) {
            const auto& pool_pos = dups[raw.index[i]];
            Index_ pool_end = pool_pos.first + pool_pos.second;

            if (raw.value) {
                for (Index_ j = pool_pos.first; j < pool_end; ++j) {
                    sortspace.emplace_back(pool[j], raw.value[i]);
                }
            } else {
                for (Index_ j = pool_pos.first; j < pool_end; ++j) {
                    sortspace.emplace_back(pool[j], 0);
                }
            }
        }

        std::sort(sortspace.begin(), sortspace.end());

        if (report_index) {
            auto icopy = ibuffer;
            for (const auto& x : sortspace) {
                *icopy= x.first;
                ++icopy;
            }
        } else {
            ibuffer = NULL;
        }

        if (raw.value) {
            auto vcopy = vbuffer;
            for (const auto& x : sortspace) {
                *vcopy = x.second;
                ++vcopy;
            }
        } else {
            vbuffer = NULL;
        }

        return SparseRange<Value_, Index_>(sortspace.size(), vbuffer, ibuffer);
    }

    /**************************************************
     ************ Full parallel extraction ************
     **************************************************/
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
        FullParallelExtractor(const DelayedSubset* p, const Options& opt) : parent(p) {
            this->full_length = parent->indices.size();
            this->internal = parent->create_inner_extractor<sparse_>(opt, parent->unique_and_sorted); // copy is deliberate.
        }

    protected:
        const DelayedSubset* parent;
    };

    struct DenseFullParallelExtractor : public FullParallelExtractor<false>, public DenseBase {
        DenseFullParallelExtractor(const DelayedSubset* p, const Options& opt) : 
            FullParallelExtractor<false>(p, opt), DenseBase(this->internal->index_length)
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, this->parent->reverse_mapping);
        }
    };

    struct SparseFullParallelExtractor : public FullParallelExtractor<true>, public SparseBase {
        SparseFullParallelExtractor(const DelayedSubset* p, const Options& opt) : 
            FullParallelExtractor<true>(p, opt), SparseBase(opt, this->internal->index_length)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            const auto& dups = this->parent->mapping_duplicates;
            const auto& pool = this->parent->mapping_duplicates_pool;
            if (this->needs_sort) {
                return DelayedSubset::reorganize_sparse_sorted(i, vbuffer, ibuffer, this->sortspace, this->itemp, this->report_index, this->internal.get(), dups, pool);
            } else {
                return DelayedSubset::reorganize_sparse_unsorted(i, vbuffer, ibuffer, this->vtemp, this->itemp, this->report_index, this->internal.get(), dups, pool);
            }
        }
    };

    /***************************************************
     ************ Block parallel extraction ************
     ***************************************************/
private:
    void transplant_indices(
        std::vector<Index_>& local, 
        std::vector<std::pair<storage_type, Index_> >& collected, 
        std::vector<Index_>& reverse_mapping) 
    const {
        std::sort(collected.begin(), collected.end());

        reverse_mapping.resize(collected.size());
        local.reserve(collected.size());

        for (const auto& current : collected) {
            if (local.empty() || current.first != local.back()) {
                local.push_back(current.first);
            }
            reverse_mapping[current.second] = local.size() - 1;
        }
    }

    void transplant_indices(
        std::vector<Index_>& local,
        std::vector<std::pair<storage_type, Index_> >& collected, 
        std::vector<std::pair<Index_, Index_> >& dups,
        std::vector<Index_>& pool) 
    const {
        std::sort(collected.begin(), collected.end());

        Index_ mapping_dim = margin_ == 0 ? mat->nrow() : mat->ncol();
        dups.resize(mapping_dim);
        pool.reserve(collected.size());
        local.reserve(collected.size());

        for (const auto& current : collected) {
            auto& range = dups[current.first];
            if (local.empty() || current.first != local.back()) {
                local.push_back(current.first);
                range.first = pool.size();
            }
            ++range.second;
            pool.push_back(current.second);
        }
    }

private:
    template<bool sparse_>
    struct BlockParallelExtractor : public ParallelExtractor<DimensionSelectionType::BLOCK, sparse_> {
        BlockParallelExtractor(const DelayedSubset* parent, const Options& opt, Index_ bs, Index_ bl) {
            this->block_start = bs;
            this->block_length = bl;

            const auto& parent_indices = parent->indices;
            std::vector<std::pair<storage_type, Index_> > collected;
            collected.reserve(bl);

            auto block_end = bs + bl;
            for (Index_ i = bs; i < block_end; ++i) {
                if constexpr(sparse_) {
                    collected.emplace_back(parent_indices[i], i);
                } else {
                    collected.emplace_back(parent_indices[i], i - bs);
                }
            }

            std::vector<Index_> local;
            if constexpr(sparse_) {
                parent->transplant_indices(local, collected, mapping_duplicates, mapping_duplicates_pool);
            } else {
                parent->transplant_indices(local, collected, reverse_mapping);
            }

            this->internal = parent->create_inner_extractor<sparse_>(opt, std::move(local));
        }

    protected:
        typename std::conditional<!sparse_, std::vector<Index_>, bool>::type reverse_mapping;
        typename std::conditional<sparse_, std::vector<std::pair<Index_, Index_> >, bool>::type mapping_duplicates;
        typename std::conditional<sparse_, std::vector<Index_>, bool>::type mapping_duplicates_pool;
    };

    struct DenseBlockParallelExtractor : public BlockParallelExtractor<false>, public DenseBase {
        DenseBlockParallelExtractor(const DelayedSubset* p, const Options& opt, Index_ bs, Index_ bl) :
            BlockParallelExtractor<false>(p, opt, bs, bl), DenseBase(this->internal->index_length) 
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto raw = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(raw, buffer, this->reverse_mapping);
        }
    };

    struct SparseBlockParallelExtractor : public BlockParallelExtractor<true>, public SparseBase {
        SparseBlockParallelExtractor(const DelayedSubset* p, const Options& opt, Index_ bs, Index_ bl) :
            BlockParallelExtractor<true>(p, opt, bs, bl), SparseBase(opt, this->internal->index_length)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            const auto& dups = this->mapping_duplicates;
            const auto& pool = this->mapping_duplicates_pool;
            if (this->needs_sort) {
                return DelayedSubset::reorganize_sparse_sorted(i, vbuffer, ibuffer, this->sortspace, this->itemp, this->report_index, this->internal.get(), dups, pool);
            } else {
                return DelayedSubset::reorganize_sparse_unsorted(i, vbuffer, ibuffer, this->vtemp, this->itemp, this->report_index, this->internal.get(), dups, pool);
            }
        }
    };

    /***************************************************
     ************ Index parallel extraction ************
     ***************************************************/
private:
    template<bool sparse_>
    struct IndexParallelExtractor : public ParallelExtractor<DimensionSelectionType::INDEX, sparse_> {
        IndexParallelExtractor(const DelayedSubset* parent, const Options& opt, std::vector<Index_> idx) {
            Index_ il = idx.size();
            this->index_length = il;
            indices = std::move(idx);

            const auto& parent_indices = parent->indices;
            std::vector<std::pair<storage_type, Index_> > collected;
            collected.reserve(il);
            for (Index_ i = 0; i < il; ++i) {
                if constexpr(sparse_) {
                    collected.emplace_back(parent_indices[indices[i]], indices[i]);
                } else {
                    collected.emplace_back(parent_indices[indices[i]], i);
                }
            }

            std::vector<Index_> local;
            if constexpr(sparse_) {
                parent->transplant_indices(local, collected, mapping_duplicates, mapping_duplicates_pool);
            } else {
                parent->transplant_indices(local, collected, reverse_mapping);
            }
            this->internal = parent->create_inner_extractor<sparse_>(opt, std::move(local));
        }

        const Index_* index_start() const {
            return indices.data();
        }

    protected:
        std::vector<Index_> indices;

        typename std::conditional<!sparse_, std::vector<Index_>, bool>::type reverse_mapping;
        typename std::conditional<sparse_, std::vector<std::pair<Index_, Index_> >, bool>::type mapping_duplicates;
        typename std::conditional<sparse_, std::vector<Index_>, bool>::type mapping_duplicates_pool;
    };

    struct DenseIndexParallelExtractor : public IndexParallelExtractor<false>, public DenseBase {
        DenseIndexParallelExtractor(const DelayedSubset* p, const Options& opt, std::vector<Index_> idx) : 
            IndexParallelExtractor<false>(p, opt, std::move(idx)), DenseBase(this->internal->index_length) 
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto raw = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(raw, buffer, this->reverse_mapping);
        }
    };

    struct SparseIndexParallelExtractor : public IndexParallelExtractor<true>, public SparseBase {
        SparseIndexParallelExtractor(const DelayedSubset* p, const Options& opt, std::vector<Index_> idx) : 
            IndexParallelExtractor<true>(p, opt, std::move(idx)), SparseBase(opt, this->internal->index_length)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            const auto& dups = this->mapping_duplicates;
            const auto& pool = this->mapping_duplicates_pool;
            if (this->needs_sort) {
                return DelayedSubset::reorganize_sparse_sorted(i, vbuffer, ibuffer, this->sortspace, this->itemp, this->report_index, this->internal.get(), dups, pool);
            } else {
                return DelayedSubset::reorganize_sparse_unsorted(i, vbuffer, ibuffer, this->vtemp, this->itemp, this->report_index, this->internal.get(), dups, pool);
            }
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
            output.reset(new SparseFullParallelExtractor(this, options));
        } else {
            output.reset(new DenseFullParallelExtractor(this, options));
        }
        return output;
    }

    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > populate_parallel(const Options& options, Index_ bs, Index_ bl) const {
        std::unique_ptr<Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new SparseBlockParallelExtractor(this, options, bs, bl));
        } else {
            output.reset(new DenseBlockParallelExtractor(this, options, bs, bl));
        }
        return output;
    }

    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > populate_parallel(const Options& options, std::vector<Index_> idx) const {
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new SparseIndexParallelExtractor(this, options, std::move(idx)));
        } else {
            output.reset(new DenseIndexParallelExtractor(this, options, std::move(idx)));
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
