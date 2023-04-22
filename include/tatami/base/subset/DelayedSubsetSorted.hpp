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
            for (size_t i = 1, end = indices.size(); i < end; ++i) {
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
            auto& len = duplicates_length[curdex];
            if (len == 0) {
                unique.push_back(curdex);
                duplicates_start[curdex] = i;
                ++ucount;
            }
            reverse_mapping.push_back(ucount);
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
            return nrow();
        } else {
            return ncol();
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

    bool prefer_rows() const {
        return mat->prefer_rows();
    }

    std::pair<double, double> dimension_preference() const {
        return mat->dimension_preference();
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
        size_t counter = 0;

        for (size_t i = 0; i < raw.number; ++i) {
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
                    // For the indexed extraction case, see SparseIndexParallelExtractor::fetch().
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
    void create_inner_extractor(const Options<Index_>& opt, const Index_* start, size_t length) {
        if constexpr(sparse_) {
            if (!opt.sparse.extract_index) {
                // Need to force extraction of indices for correct duplication.
                auto copy = opt;
                copy.sparse.extract_index = true;
                return new_extractor<margin_ != 0, sparse_>(mat.get(), start, length, copy);
            }
        }

        return new_extractor<margin_ != 0, sparse_>(mat.get(), start, length, opt);
    }

    struct DenseBase {
        DenseBase(Index_ n) : temp(n) {}
    protected:
        std::vector<Value_> temp;
    };

    struct SparseBase {
        SparseBase(const Options<Index_>& opt, Index_ n) : 
            // Only need vtemp if we want to extract values.
            vtemp(opt.sparse.extract_value ? n : 0), 

            // Always need itemp to get indices back for duplicate counting + expansion.
            itemp(n), 

            // Need to store a flag indicating whether the indices are to be reported, though.
            report_index(opt.sparse.extract_index) 
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
    template<bool sparse_>
    struct FullParallelExtractor : public Extractor<DimensionSelectionType::FULL, sparse_, Value_, Index_> {
        FullParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt) : parent(p) {
            this->full_length = parent->indices.size();
            internal = create_inner_extractor(opt, parent->unique.data(), parent->unique.size());
        }

    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        const DelayedSubsetSorted* parent;
    };

    struct DenseFullParallelExtractor : public FullParallelExtractor<false>, public DenseBase {
        DenseFullParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt) : 
            FullParallelExtractor<false>(p, opt), DenseBase(this->internal->index_length)
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = internal->fetch(i, this->temp.data());
            return remap_dense(ref, buffer, this->parent->reverse_mapping);
        }
    };

    struct SparseFullParallelExtractor : public FullParallelExtractor<true>, public SparseBase {
        SparseFullParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt) : 
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
    struct BlockParallelExtractor : public Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> {
        BlockExtractor(const DelayedSubsetSorted* parent, const Options<Index_>& opt, Index_ bs, Index_ bl) {
            this->block_start = bs;
            this->block_length = bl;

            if (bl) {
                const auto& punique = parent->unique;
                auto pstart = punique.begin();
                auto pend = punique.end();

                // Finding the edges of the 'unique' subset so that we
                // don't have to allocate another vector here.
                auto left = parent->indices[bs];
                from = std::lower_bound(pstart, pend, left) - pstart;

                auto right = parent->indices[bs + bl - 1];
                Index_ to = std::upper_bound(pstart + from, pend, right) - pstart;

                internal = create_inner_extractor(opt, punique.data() + from, to - from);
            }
        }

    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        Index_ from = 0; // needed for the dense constructor, oh well.
    };

    struct DenseBlockParallelExtractor : public BlockParallelExtractor<false>, public DenseBase {
        DenseBlockParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt, Index_ bs, Index_ bl) : 
            BlockParallelExtractor<false>(p, opt, bs, bl), DenseBase(this->internal->index_length)
        {
            if (bl) {
                // Adjusting the mappings to account for the truncation to the indexing subset used to construct 'inner'.
                const auto& revmap = this->parent->reverse_mapping;
                reverse_mapping.reserve(bl);
                for (Index_ b = 0; b < bl; ++b) {
                    reverse_mapping.push_back(revmap[b + bs] - this->from);
                }
            }
        }

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = internal->fetch(i, this->temp.data());
            return remap_dense(ref, buffer, reverse_mapping);
        }

    protected:
        std::vector<Index_> reverse_mapping;
    };

    struct SparseBlockParallelExtractor : public BlockParallelExtractor<true>, public SparseBase {
        SparseBlockParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt, Index_ bs, Index_ bl) : 
            BlockParallelExtractor<true>(p, opt, bs, bl), SparseBase(opt, this->internal->index_length)
        {
            if (bl) {
                auto left = parent->indices[bs];
                auto right = parent->indices[bs + bl - 1];
                auto mapping_dim = get_mapping_dim();

                duplicate_starts.resize(mapping_dim);
                const auto& dups = this->parent->duplicate_starts;
                std::copy(dups.begin() + left, dups.begin() + right + 1, duplicate_starts.begin());

                duplicate_lengths.resize(mapping_dim);
                const auto& dupl = this->parent->duplicate_lengths;
                std::copy(dupl.begin() + left, dupl.begin() + right + 1, duplicate_lengths.begin());

                // Now we have to adjust for the loss of duplicates at the block boundaries.
                Index_ lx = bs;
                while (lx > 0) {
                    --lx;
                    if (parent->indices[lx] != left) {
                        break;
                    }
                    --duplicate_lengths[left];
                    ++duplicate_starts[left];
                }

                Index_ rx = bs + bl;
                Index_ rend = parent->indices.size();
                while (rx < rend) {
                    if (parent->indices[rx] != right) {
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
    struct IndexParallelExtractor : public Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> {
        const Index_* index_start() const {
            return indices.data();
        }
    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        std::vector<Index_> indices;
    }

    struct DenseIndexParallelExtractorRaw : public IndexParallelExtractor<false> {
        DenseIndexParallelExtractorRaw(const DelayedSubsetSorted* p, const Options<Index_>& opt, const Index_* is, size_t il) {
            const auto& pindices = this->parent->indices;
            this->indices.reserve(il);
            reverse_mapping.reserve(il);

            Index_ ucount = 0;
            for (size_t i = 0; i < il; ++i) {
                Index_ curdex = is[i];
                if (indices.empty() || indices.back() == curdex) {
                    this->indices.push_back(curdex);
                    ++ucount;
                }
                reverse_mapping.push_back(ucount);
            }

            this->internal = create_inner_extractor(opt, this->indices.data(), this->indices.size());
            this->index_length = il;
            this->indices = std::vector<Index_>(is, is + il);
        }
    protected:
        std::vector<Index_> reverse_mapping;
    };

    struct DenseIndexParallelExtractor : public DenseIndexParallelExtractorRaw, public DenseBase {
        DenseIndexParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt, const Index_* is, size_t il) :
            DenseIndexParallelExtractorRaw(p, opt, is, il), DenseBase(this->internal->index_length) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ref = this->internal->fetch(i, this->temp.data());
            return remap_dense(ref, buffer, this->reverse_mapping);
        }
    };

    struct SparseIndexParallelExtractorRaw : public IndexParallelExtractor<true> {
        SparseIndexParallelExtractorRaw(const DelayedSubsetSorted* p, const Options<Index_>& opt, const Index_* is, Index_ il) {
            const auto& pindices = this->parent->indices;
            this->indices.reserve(il);
            Index_ mapping_dim = get_mapping_dim();
            duplicate_starts.resize(mapping_dim);
            duplicate_lengths.resize(mapping_dim);

            for (size_t i = 0; i < il; ++i) {
                Index_ curdex = pindices[is[i]];
                auto& len = duplicate_lengths[curdex];
                if (len == 0) {
                    this->indices.push_back(curdex);
                    duplicate_starts[curdex] = i; // references a range on the this->indices array, see the remap_spares_duplicates() call below.
                }
                ++len;
            }

            this->internal = create_inner_extractor(opt, this->indices.data(), this->indices.size());
            this->index_length = il;
            this->indices = std::vector<Index_>(is, is + il);
        }
    protected:
        std::vector<Index_> duplicate_starts;
        std::vector<Index_> duplicate_lengths;
    };

    struct SparseIndexParallelExtractor : public SparseIndexParallelExtractorRaw, public SparseBase {
        SparseIndexParallelExtractor(const DelayedSubsetSorted* p, const Options<Index_>& opt, const Index_* is, Index_ il) :
            SparseIndexParallelExtractorRaw(p, opt, is, il), SparseBase(opt, this->internal->index_length) {}

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
     ************ Perpendicular extraction ************
     **************************************************/
private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct PerpendicularExtractor : public Extractor<selection_, sparse_, Value_, Index_> {
        PerpendicularExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > i, const DelayedSubsetSorted* p) : 
            internal(std::move(i)), parent(p)
        {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = internal->full_length;
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = internal->block_start;
                this->block_length = internal->block_length;
            } else {
                this->index_length = internal->index_length;
            }
        }

        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return internal->index_start();
            } else {
                return NULL;
            }
        }

    protected:
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > internal;
        const DelayedSubsetSorted* parent;
    };

    template<DimensionSelectionType selection_>
    struct DensePerpendicularExtractor : public Extractor<selection_, false, Value_, Index_> {
        DensePerpendicularExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > i, const DelayedSubsetSorted* p) : 
            PerpendicularExtractor(std::move(i), p) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            return this->internal->fetch(this->parent->indices[i], buffer);
        }
    };

    template<DimensionSelectionType selection_>
    struct SparsePerpendicularExtractor : public Extractor<selection_, true, Value_, Index_> {
        SparsePerpendicularExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > i, const DelayedSubsetSorted* p) : 
            PerpendicularExtractor(std::move(i), p) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return this->internal->fetch(this->parent->indices[i], vbuffer, ibuffer);
        }
    };

    /**************************************************
     ************ Public virtual overrides ************
     **************************************************/
private:
    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::FULL, sparse_, Value_, Index_> > populate_parallel(const Options<Index_>& options) const {
        std::unique_ptr<Extractor<DimensionSelectionType::FULL, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new SparseFullParallelWorkspace(options));
        } else {
            output.reset(new DenseFullParallelWorkspace(options));
        }
        return output;
    }

    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > populate_parallel(const Options<Index_>& options, Index_ bs, Index_ bl) const {
        std::unique_ptr<Extractor<DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new SparseBlockParallelWorkspace(options, bs, bl));
        } else {
            output.reset(new DenseBlockParallelWorkspace(options, bs, bl));
        }
        return output;
    }

    template<bool sparse_>
    std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > populate_parallel(const Options<Index_>& options, const Index_* is, size_t il) const {
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > output;
        if constexpr(sparse_) {
            output.reset(new SparseIndexParallelWorkspace(options, bs, bl));
        } else {
            output.reset(new DenseIndexParallelWorkspace(options, bs, bl));
        }
        return output;
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options<Index_>& options, Args_... args) const {
        if constexpr(accrow_ == (margin_ == 0)) {
            // TODO: handle variable access patterns here.
            std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;
            if constexpr(sparse_) {
                output.reset(new SparsePerpendicularExtractor(new_extractor<accrow_, sparse_>(mat.get(), args...), this));
            } else {
                output.reset(new DensePerpendicularExtractor(new_extractor<accrow_, sparse_>(mat.get(), args...), this));
            }
            return output;
        } else {
            return populate_parallel(options, args...);
        }
    }

public:
    std::unique_ptr<DenseFullExtractor<Value_, Index_> > dense_row(const Options<Index_>& options) const {
        return populate<true, DimensionSelectionType::FULL, false>(options);
    }

    std::unique_ptr<DenseBlockExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options<Index_>& options) const {
        return populate<true, DimensionSelectionType::BLOCK, false>(options, block_start, block_length);
    }

    std::unique_ptr<DenseIndexExtractor<Value_, Index_> > dense_row(const Index_* index_start, size_t index_length, const Options<Index_>& options) const {
        return populate<true, DimensionSelectionType::INDEX, false>(options, index_start, index_length);
    }

    std::unique_ptr<DenseFullExtractor<Value_, Index_> > dense_column(const Options<Index_>& options) const {
        return populate<false, DimensionSelectionType::FULL, false>(options);
    }

    std::unique_ptr<DenseBlockExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options<Index_>& options) const {
        return populate<false, DimensionSelectionType::BLOCK, false>(options, block_start, block_length);
    }

    std::unique_ptr<DenseIndexExtractor<Value_, Index_> > dense_column(const Index_* index_start, size_t index_length, const Options<Index_>& options) const {
        return populate<false, DimensionSelectionType::INDEX, false>(options, index_start, index_length);
    }

public:
    std::unique_ptr<SparseFullExtractor<Value_, Index_> > dense_row(const Options<Index_>& options) const {
        return populate<true, DimensionSelectionType::FULL, true>(options);
    }

    std::unique_ptr<SparseBlockExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options<Index_>& options) const {
        return populate<true, DimensionSelectionType::BLOCK, true>(options, block_start, block_length);
    }

    std::unique_ptr<SparseIndexExtractor<Value_, Index_> > dense_row(const Index_* index_start, size_t index_length, const Options<Index_>& options) const {
        return populate<true, DimensionSelectionType::INDEX, true>(options, index_start, index_length);
    }

    std::unique_ptr<SparseFullExtractor<Value_, Index_> > dense_column(const Options<Index_>& options) const {
        return populate<false, DimensionSelectionType::FULL, true>(options);
    }

    std::unique_ptr<SparseBlockExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options<Index_>& options) const {
        return populate<false, DimensionSelectionType::BLOCK, true>(options, block_start, block_length);
    }

    std::unique_ptr<SparseIndexExtractor<Value_, Index_> > dense_column(const Index_* index_start, size_t index_length, const Options<Index_>& options) const {
        return populate<false, DimensionSelectionType::INDEX, true>(options, index_start, index_length);
    }
};

}

#endif
