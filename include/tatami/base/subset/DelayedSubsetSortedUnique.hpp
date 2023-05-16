#ifndef TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP
#define TATAMI_DELAYED_SUBSET_SORTED_UNIQUE_HPP

#include "utils.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetSortedUnique.hpp
 *
 * @brief Delayed subsetting with sorted and unique row/column indices.
 */

namespace tatami {

/**
 * @brief Delayed subsetting of a matrix with sorted, unique indices.
 *
 * Implements delayed subsetting (i.e., slicing) on the rows or columns of a matrix, given a vector of sorted and unique indices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam IndexStorage_ Vector containing the subset indices.
 */
template<int margin_, typename Value_, typename Index_, class IndexStorage_>
class DelayedSubsetSortedUnique : public Matrix<Value_, Index_> {
    static constexpr bool storage_has_data = has_data<Index_, IndexStorage_>::value;
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param idx Vector of 0-based indices to use for subsetting on the rows (if `margin_ = 0`) or columns (if `margin_ = 1`).
     * This should be sorted and unique.
     * @param check Whether to check `idx` for sorted and unique values.
     */
    DelayedSubsetSortedUnique(std::shared_ptr<const Matrix<Value_, Index_> > p, IndexStorage_ idx, bool check = true) : mat(std::move(p)) {
        if constexpr(storage_has_data) {
            indices = std::move(idx);
        } else {
            indices = std::vector<Index_>(idx.begin(), idx.end());
        }

        if (check) {
            for (Index_ i = 1, end = indices.size(); i < end; ++i) {
                if (indices[i] <= indices[i-1]) {
                    throw std::runtime_error("indices should be unique and sorted");
                }
            }
        }

        Index_ mapping_dim = margin_ == 0 ? mat->nrow() : mat->ncol();
        mapping_single.resize(mapping_dim);
        for (Index_ i = 0, end = indices.size(); i < end; ++i) {
            mapping_single[indices[i]] = i;
        }
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    typename std::conditional<storage_has_data, IndexStorage_, std::vector<Index_> >::type indices;
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

    /*********************************************
     ************ Parallel extraction ************
     *********************************************/
private:
    template<DimensionSelectionType selection_, bool sparse_>
    struct ParallelWorkspaceBase : public Extractor<selection_, sparse_, Value_, Index_> {
        ParallelWorkspaceBase(const DelayedSubsetSortedUnique* parent, const Options& opt) {
            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->full_length = parent->indices.size();
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::vector<Index_>(parent->indices.begin(), parent->indices.end()), opt); // copy of indices is deliberate.
            }
        }

        ParallelWorkspaceBase(const DelayedSubsetSortedUnique* parent, const Options& opt, Index_ bs, Index_ bl) {
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;

                auto pistart = parent->indices.begin() + bs;
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::vector<Index_>(pistart, pistart + bl), opt);
            }
        }

        ParallelWorkspaceBase(const DelayedSubsetSortedUnique* parent, const Options& opt, std::vector<Index_> idx) {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                Index_ il = idx.size();
                this->index_length = il;
                indices = std::move(idx);

                // Reusing 'indices' to store the inner indices.
                std::vector<Index_> local;
                local.reserve(il);
                for (auto i : indices) {
                    local.push_back(parent->indices[i]);
                }
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::move(local), opt);
            }
        }

    protected:
        std::unique_ptr<Extractor<DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        typename std::conditional<selection_ == DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            internal->set_oracle(std::move(o));
        }
    };

    template<DimensionSelectionType selection_>
    struct DenseParallelWorkspace : public ParallelWorkspaceBase<selection_, false> {
        template<typename ... Args_>
        DenseParallelWorkspace(const DelayedSubsetSortedUnique* parent, const Options& opt, Args_&& ... args) : 
            ParallelWorkspaceBase<selection_, false>(parent, opt, std::forward<Args_>(args)...) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            return this->internal->fetch(i, buffer);
        }
    };

    template<DimensionSelectionType selection_>
    struct SparseParallelWorkspace : public ParallelWorkspaceBase<selection_, true> {
        template<typename ... Args_>
        SparseParallelWorkspace(const DelayedSubsetSortedUnique* p, const Options& opt, Args_&& ... args) : 
            ParallelWorkspaceBase<selection_, true>(p, opt, std::forward<Args_>(args)...), parent(p) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto raw = this->internal->fetch(i, vbuffer, ibuffer);
            if (raw.index) {
                auto icopy = ibuffer;
                for (Index_ i = 0; i < raw.number; ++i, ++icopy) {
                    *icopy = parent->mapping_single[raw.index[i]];
                }
                raw.index = ibuffer;
            } else {
                raw.index = NULL;
            }
            return raw;
        }

    protected:
        const DelayedSubsetSortedUnique* parent;
    };

    /******************************************
     ************ Public overrides ************
     ******************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate(const Options& opt, Args_&& ... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(accrow_ == (margin_ == 0)) {
            // TODO: fiddle with the access limits in 'opt'.
            return subset_utils::populate_perpendicular<accrow_, selection_, sparse_>(mat.get(), indices, opt, std::forward<Args_>(args)...);
        } else {
            if constexpr(sparse_) {
                output.reset(new SparseParallelWorkspace<selection_>(this, opt, std::forward<Args_>(args)...));
            } else {
                output.reset(new DenseParallelWorkspace<selection_>(this, opt, std::forward<Args_>(args)...));
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

}

#endif
