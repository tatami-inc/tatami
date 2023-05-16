#ifndef TATAMI_DELAYED_SUBSET_SORTED_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetSorted.hpp
 *
 * @brief Delayed subsetting with sorted row/column indices.
 *
 * This is equivalent to the `DelayedSubset` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed subsetting of a matrix with sorted indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of sorted indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetSorted : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * This should be sorted, but may be duplicated.
     * @param check Whether to check `idx` for sorted values.
     */
    DelayedSubsetSorted(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : mat(std::move(p)), indices(std::move(idx)) {
        if (check) {
            for (Index_ i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] < indices[i-1]) {
                    throw std::runtime_error("indices should be sorted");
                }
            }
        }

        Index_ mapping_dim = get_mapping_dim();
        unique.reserve(indices.size());
        reverse_mapping.reserve(indices.size());
        duplicate_starts.resize(mapping_dim);
        duplicate_lengths.resize(mapping_dim);

        Index_ ucount = 0;
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            Index_ curdex = indices[i];
            auto& len = duplicate_lengths[curdex];
            if (len == 0) {
                unique.push_back(curdex);
                duplicate_starts[curdex] = i;
                ++ucount;
            }
            reverse_mapping.push_back(ucount - 1);
            ++len;
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    IndexStorage_ indices;

    std::vector<Index_> unique;
    std::vector<Index_> reverse_mapping;
    std::vector<Index_> duplicate_starts; // holds the start position of each duplicate stretch on 'indices'.
    std::vector<Index_> duplicate_lengths; // holds the length of each duplicate stretch on 'indices'.

    Index_ get_mapping_dim() const {
        if constexpr(margin_ == 0) {
            return mat->nrow();
        } else {
            return mat->ncol();
        }
    }

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
    template<class Extractor_, class Indexer_>
    static SparseRange<Value_, Index_> remap_sparse_duplicates(
        Index_ i, 
        Value_* vbuffer, 
        Index_* ibuffer, 
        std::vector<Value_>& vtemp,
        std::vector<Index_>& itemp,
        bool report_index,
        Extractor_* internal,
        const std::vector<Index_>& dup_starts,
        const std::vector<Index_>& dup_lengths,
        Indexer_ custom_index) 
    {
        // Allocation status depends on the extraction mode used to construct internal; see SparseBase.
        Value_* vin = vtemp.data();

        // This should always be allocated, as we need this to get the expanded counts; see SparseBase.
        Index_* iin = itemp.data();

        auto raw = internal->fetch(i, vin, iin);
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
            auto len = dup_lengths[raw.index[i]];
            counter += len;

            if (vcopy) {
                std::fill(vcopy, vcopy + len, raw.value[i]);
                vcopy += len;
            }

            if (icopy) {
                if constexpr(std::is_same<Indexer_, bool>::value) {
                    std::iota(icopy, icopy + len, dup_starts[raw.index[i]]);
                } else {
                    // For the indexed extraction case, see IndexSparseParallelExtractor::fetch().
                    auto custom_start = custom_index + dup_starts[raw.index[i]];
                    std::copy(custom_start, custom_start + len, icopy);
                }
                icopy += len;
            }
        }

        return SparseRange<Value_, Index_>(counter, vbuffer, ibuffer);
    }

private:
    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > create_inner_extractor(const Options& opt, std::vector<Index_> idx) const {
        if constexpr(sparse_) {
            if (!opt.sparse_extract_index) {
                // Need to force extraction of indices for correct duplication.
                auto copy = opt;
                copy.sparse_extract_index = true;
                return new_extractor<margin_ != 0, sparse_>(mat.get(), std::move(idx), copy);
            }
        }

        return new_extractor<margin_ != 0, sparse_>(mat.get(), std::move(idx), opt);
    }

    struct DenseBase {
        DenseBase(size_t n) : temp(n) {}
    protected:
        std::vector<Value_> temp;
    };

    struct SparseBase {
        SparseBase(const Options& opt, size_t n) : 
            // Only need vtemp if we want to extract values.
            vtemp(opt.sparse_extract_value ? n : 0), 

            // Always need itemp to get indices back for duplicate counting + expansion.
            itemp(n), 

            // Need to store a flag indicating whether the indices are to be reported, though.
            report_index(opt.sparse_extract_index) 
        {}
    protected:
        std::vector<Value_> vtemp;
        std::vector<Index_> itemp;
        bool report_index;
    };

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
        FullParallelExtractor(const DelayedSubsetSorted* p, const Options& opt) : parent(p) {
            this->full_length = parent->indices.size();
            this->internal = parent->create_inner_extractor<sparse_>(opt, parent->unique); // copy is deliberate here.
        }

    protected:
        const DelayedSubsetSorted* parent;
    };

    struct FullDenseParallelExtractor : public FullParallelExtractor<false>, public DenseBase {
        FullDenseParallelExtractor(const DelayedSubsetSorted* p, const Options& opt) : 
            FullParallelExtractor<false>(p, opt), DenseBase(this->internal->index_length)
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, this->parent->reverse_mapping);
        }
    };

    struct FullSparseParallelExtractor : public FullParallelExtractor<true>, public SparseBase {
        FullSparseParallelExtractor(const DelayedSubsetSorted* p, const Options& opt) : 
            FullParallelExtractor<true>(p, opt), SparseBase(opt, this->internal->index_length) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return remap_sparse_duplicates(
                i, vbuffer, ibuffer, 
                this->vtemp, this->itemp, this->report_index, this->internal.get(),
                this->parent->duplicate_starts, this->parent->duplicate_lengths,
                false
            );
        }
    };

    /***************************************************
     ************ Block parallel extraction ************
     ***************************************************/
private:
    template<bool sparse_>
    struct BlockParallelExtractor : public ParallelExtractor<DimensionSelectionType::BLOCK, sparse_> {
        BlockParallelExtractor(const DelayedSubsetSorted* parent, const Options& opt, Index_ bs, Index_ bl) {
            this->block_start = bs;
            this->block_length = bl;

            const auto& punique = parent->unique;
            auto pstart = punique.begin();
            auto pend = punique.end();

            Index_ to = 0;
            if (bl) {
                // Finding the edges of the 'unique' subset so that we
                // don't have to allocate another vector here.
                auto left = parent->indices[bs];
                from = std::lower_bound(pstart, pend, left) - pstart;

                auto right = parent->indices[bs + bl - 1];
                to = std::upper_bound(pstart + from, pend, right) - pstart;
            }

            this->internal = parent->create_inner_extractor<sparse_>(opt, std::vector<Index_>(pstart + from, pstart + to));
        }

    protected:
        Index_ from = 0; // needed for the dense constructor, oh well.
    };

    struct BlockDenseParallelExtractor : public BlockParallelExtractor<false>, public DenseBase {
        BlockDenseParallelExtractor(const DelayedSubsetSorted* parent, const Options& opt, Index_ bs, Index_ bl) : 
            BlockParallelExtractor<false>(parent, opt, bs, bl), DenseBase(this->internal->index_length)
        {
            if (bl) {
                // Adjusting the mappings to account for the truncation to the indexing subset used to construct 'inner'.
                const auto& revmap = parent->reverse_mapping;
                reverse_mapping.reserve(bl);
                for (Index_ b = 0; b < bl; ++b) {
                    reverse_mapping.push_back(revmap[b + bs] - this->from);
                }
            }
        }

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, reverse_mapping);
        }

    protected:
        std::vector<Index_> reverse_mapping;
    };

    struct BlockSparseParallelExtractor : public BlockParallelExtractor<true>, public SparseBase {
        BlockSparseParallelExtractor(const DelayedSubsetSorted* parent, const Options& opt, Index_ bs, Index_ bl) : 
            BlockParallelExtractor<true>(parent, opt, bs, bl), SparseBase(opt, this->internal->index_length)
        {
            if (bl) {
                const auto& parent_indices = parent->indices;

                auto left = parent_indices[bs];
                auto right = parent_indices[bs + bl - 1];
                auto mapping_dim = parent->get_mapping_dim();

                duplicate_starts.resize(mapping_dim);
                const auto& dups = parent->duplicate_starts;
                std::copy(dups.begin() + left, dups.begin() + right + 1, duplicate_starts.begin() + left);

                duplicate_lengths.resize(mapping_dim);
                const auto& dupl = parent->duplicate_lengths;
                std::copy(dupl.begin() + left, dupl.begin() + right + 1, duplicate_lengths.begin() + left);

                // Now we have to adjust for the loss of duplicates at the block boundaries.
                Index_ lx = bs;
                while (lx > 0) {
                    --lx;
                    if (parent_indices[lx] != left) {
                        break;
                    }
                    --duplicate_lengths[left];
                    ++duplicate_starts[left];
                }

                Index_ rx = bs + bl;
                Index_ rend = parent_indices.size();
                while (rx < rend) {
                    if (parent_indices[rx] != right) {
                        break;
                    }
                    --duplicate_lengths[right];
                    ++rx;
                }
            }
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return remap_sparse_duplicates(
                i, vbuffer, ibuffer, 
                this->vtemp, this->itemp, this->report_index, this->internal.get(),
                this->duplicate_starts, this->duplicate_lengths,
                false
            );
        }

    protected:
        std::vector<Index_> duplicate_starts;
        std::vector<Index_> duplicate_lengths;
    };

    /***************************************************
     ************ Index parallel extraction ************
     ***************************************************/
private:
    template<bool sparse_>
    struct IndexParallelExtractor : public ParallelExtractor<DimensionSelectionType::INDEX, sparse_> {
        IndexParallelExtractor(const DelayedSubsetSorted* parent, const Options& opt, std::vector<Index_> idx) {
            Index_ il = idx.size();
            this->index_length = il;
            indices = std::move(idx);

            const auto& parent_indices = parent->indices;
            std::vector<Index_> local;
            local.reserve(il);

            if constexpr(!sparse_) {
                reverse_mapping.reserve(il);
                Index_ ucount = 0;
                for (Index_ i = 0; i < il; ++i) {
                    Index_ curdex = parent_indices[indices[i]];
                    if (local.empty() || local.back() != curdex) {
                        local.push_back(curdex);
                        ++ucount;
                    }
                    reverse_mapping.push_back(ucount - 1);
                }

            } else {
                Index_ mapping_dim = parent->get_mapping_dim();
                duplicate_starts.resize(mapping_dim);
                duplicate_lengths.resize(mapping_dim);

                for (Index_ i = 0; i < il; ++i) {
                    Index_ curdex = parent_indices[indices[i]];
                    auto& len = duplicate_lengths[curdex];
                    if (len == 0) {
                        local.push_back(curdex);
                        duplicate_starts[curdex] = i; // references a range on the this->indices array, see the remap_sparse_duplicates() call below.
                    }
                    ++len;
                }
            }

            this->internal = parent->create_inner_extractor<sparse_>(opt, std::move(local));
        }

        const Index_* index_start() const {
            return indices.data();
        }

    protected:
        std::vector<Index_> indices;

        typename std::conditional<!sparse_, std::vector<Index_>, bool>::type reverse_mapping;
        typename std::conditional<sparse_, std::vector<Index_>, bool>::type duplicate_starts;
        typename std::conditional<sparse_, std::vector<Index_>, bool>::type duplicate_lengths;
    };

    struct IndexDenseParallelExtractor : public IndexParallelExtractor<false>, public DenseBase {
        IndexDenseParallelExtractor(const DelayedSubsetSorted* parent, const Options& opt, std::vector<Index_> idx) :
            IndexParallelExtractor<false>(parent, opt, std::move(idx)), DenseBase(this->internal->index_length) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return subset_utils::remap_dense(ref, buffer, this->reverse_mapping);
        }
    };

    struct IndexSparseParallelExtractor : public IndexParallelExtractor<true>, public SparseBase {
        IndexSparseParallelExtractor(const DelayedSubsetSorted* parent, const Options& opt, std::vector<Index_> idx) :
            IndexParallelExtractor<true>(parent, opt, std::move(idx)), SparseBase(opt, this->internal->index_length) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return remap_sparse_duplicates(
                i, vbuffer, ibuffer, 
                this->vtemp, this->itemp, this->report_index, this->internal.get(),
                this->duplicate_starts, this->duplicate_lengths,
                this->indices.data()
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
