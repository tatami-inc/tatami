#ifndef TATAMI_DELAYED_SUBSET_BLOCK
#define TATAMI_DELAYED_SUBSET_BLOCK

#include "../base/Matrix.hpp"
#include "../base/utils.hpp"
#include "../utils/FixedOracle.hpp"

#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetBlock.hpp
 *
 * @brief Delayed subsetting to a single contiguous block.
 *
 * This is a specialized implementation that is more efficient than the `tatami::DelayedSubset` class.
 */

namespace tatami {

/**
 * @brief Delayed subsetting to a contiguous block.
 *
 * Implements delayed subsetting (i.e., slicing) of a matrix to a single contiguous block of rows or columns.
 * This is a specialized implementation that is more efficient than the `tatami::DelayedSubset` class.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam margin_ Dimension along which the subsetting is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 */
template<int margin_, typename Value_, typename Index_>
class DelayedSubsetBlock : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param s Index of the start of the block. This should be a row index if `margin_ = 0` and a column index otherwise.
     * @param l Length of the block, in terms of the number of rows (if `margin_ = 0`) or columns (otherwise).
     */
    DelayedSubsetBlock(std::shared_ptr<const Matrix<Value_, Index_> > p, Index_ s, Index_ l) : mat(std::move(p)), block_start(s), block_length(l) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    Index_ block_start, block_length;

public:
    Index_ nrow() const {
        if constexpr(margin_==0) {
            return block_length;
        } else {
            return mat->nrow();
        }
    }

    Index_ ncol() const {
        if constexpr(margin_==0) {
            return mat->ncol();
        } else {
            return block_length;
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

    /**********************************
     ***** Along-subset extractor *****
     **********************************/
private:
    template<bool oracle_, bool sparse_>
    struct AlongFullExtractor : public GeneralizedExtractor<oracle_, DimensionSelectionType::FULL, sparse_, Value_, Index_> {
        template<class OraclePointer_>
        AlongFullExtractor(const DelayedSubsetBlock* parent, OraclePointer_ oracle, bool, const Options& opt) {
            this->full_length = (margin_ == 0 ? parent->nrow() : parent->ncol());
            if constexpr(oracle_) {
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::move(oracle), parent->block_start, parent->block_length, opt);
            } else {
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), parent->block_start, parent->block_length, opt);
            }
        }

    protected:
        std::unique_ptr<GeneralizedExtractor<oracle_, DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > internal;
    };

    template<bool oracle_, bool sparse_>
    struct AlongBlockExtractor : public GeneralizedExtractor<oracle_, DimensionSelectionType::BLOCK, sparse_, Value_, Index_> {
        template<class OraclePointer_>
        AlongBlockExtractor(const DelayedSubsetBlock* parent, OraclePointer_ oracle, std::pair<Index_, Index_> block, const Options& opt) {
            this->block_start = block.first;
            this->block_length = block.second;
            if constexpr(oracle_) {
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::move(oracle), block.first + parent->block_start, block.second, opt);
            } else {
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), block.first + parent->block_start, block.second, opt);
            }
        }

    protected:
        std::unique_ptr<GeneralizedExtractor<oracle_, DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > internal;
    };

    template<bool oracle_, bool sparse_>
    struct AlongIndexExtractor : public GeneralizedExtractor<oracle_, DimensionSelectionType::INDEX, sparse_, Value_, Index_> {
        template<class OraclePointer_>
        AlongIndexExtractor(const DelayedSubsetBlock* parent, OraclePointer_ oracle, std::vector<Index_> idx, const Options& opt) : indices(std::move(idx)) {
            this->index_length = indices.size();

            // Shifting the block for the underlying matrix.
            std::vector<Index_> local = indices;
            for (auto& x : local) {
                x += parent->block_start;
            }

            if constexpr(oracle_) {
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::move(oracle), std::move(local), opt);
            } else {
                internal = new_extractor<margin_ != 0, sparse_>(parent->mat.get(), std::move(local), opt);
            }
        }

    protected:
        std::unique_ptr<GeneralizedExtractor<oracle_, DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;
        std::vector<Index_> indices;

    public:
        const Index_* index_start() const {
            return indices.data();
        }
    };

    template<bool oracle_, DimensionSelectionType selection_, bool sparse_>
    using SelectionBasedAlongExtractor = typename std::conditional<
            selection_ == DimensionSelectionType::FULL,
            AlongFullExtractor<oracle_, sparse_>,
            typename std::conditional<
                selection_ == DimensionSelectionType::BLOCK,
                AlongBlockExtractor<oracle_, sparse_>,
                AlongIndexExtractor<oracle_, sparse_>
            >::type
        >::type;

    template<DimensionSelectionType selection_>
    struct SimpleDenseAlongExtractor : public SelectionBasedAlongExtractor<false, selection_, false> {
        template<typename SelectionArg_>
        SimpleDenseAlongExtractor(const DelayedSubsetBlock* p, bool oracle, SelectionArg_&& sarg, const Options& opt) : 
            SelectionBasedAlongExtractor<false, selection_, false>(p, oracle, std::forward<SelectionArg_>(sarg), opt) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            return this->internal->fetch(i, buffer);
        }
    };

    template<DimensionSelectionType selection_>
    struct OracleAwareDenseAlongExtractor : public SelectionBasedAlongExtractor<true, selection_, false> {
        template<typename SelectionArg_>
        OracleAwareDenseAlongExtractor(const DelayedSubsetBlock* p, std::shared_ptr<Oracle<Index_> > oracle, SelectionArg_&& sarg, const Options& opt) : 
            SelectionBasedAlongExtractor<true, selection_, false>(p, std::move(oracle), std::forward<SelectionArg_>(sarg), opt) 
        {
            this->total_predictions = this->internal->total_predictions;
        }

        const Value_* fetch(Index_& i, Value_* buffer) {
            ++(this->used_predictions);
            return this->internal->fetch(i, buffer);
        }

        const Oracle<Index_>* oracle() const {
            return this->internal->oracle();
        }
    };

    template<DimensionSelectionType selection_>
    struct SimpleSparseAlongExtractor : public SelectionBasedAlongExtractor<false, selection_, true> {
        template<typename SelectionArg_>
        SimpleSparseAlongExtractor(const DelayedSubsetBlock* p, bool oracle, SelectionArg_&& sarg, const Options& opt) :
            SelectionBasedAlongExtractor<false, selection_, true>(p, oracle, std::forward<SelectionArg_>(sarg), opt), offset(p->block_start) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto out = this->internal->fetch(i, vbuffer, ibuffer);
            if (out.index && offset) {
                for (Index_ j = 0; j < out.number; ++j) {
                    ibuffer[j] = out.index[j] - offset;
                }
                out.index = ibuffer;
            }
            return out;
        }

    private:
        Index_ offset;
    };

    template<DimensionSelectionType selection_>
    struct OracleAwareSparseAlongExtractor : public SelectionBasedAlongExtractor<true, selection_, true> {
        template<typename SelectionArg_>
        OracleAwareSparseAlongExtractor(const DelayedSubsetBlock* p, std::shared_ptr<Oracle<Index_> > oracle, SelectionArg_&& sarg, const Options& opt) :
            SelectionBasedAlongExtractor<true, selection_, true>(p, std::move(oracle), std::forward<SelectionArg_>(sarg), opt), offset(p->block_start) 
        {
            this->total_predictions = this->internal->total_predictions;
        }

        SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
            auto out = this->internal->fetch(i, vbuffer, ibuffer);
            if (out.index && offset) {
                for (Index_ j = 0; j < out.number; ++j) {
                    ibuffer[j] = out.index[j] - offset;
                }
                out.index = ibuffer;
            }
            ++(this->used_predictions);
            return out;
        }

        const Oracle<Index_>* oracle() const {
            return this->internal->oracle();
        }

    private:
        Index_ offset;
    };

    template<bool oracle_, DimensionSelectionType selection_, bool sparse_>
    using AlongExtractor = typename std::conditional<
            sparse_,
            typename std::conditional<
                    oracle_, 
                    OracleAwareSparseAlongExtractor<selection_>, 
                    SimpleSparseAlongExtractor<selection_> 
                >::type,
            typename std::conditional<
                    oracle_, 
                    OracleAwareDenseAlongExtractor<selection_>, 
                    SimpleDenseAlongExtractor<selection_> 
                >::type
        >::type;

    /***********************************
     ***** Across-subset extractor *****
     ***********************************/
private:
    struct SubsetOracle : public Oracle<Index_> {
        SubsetOracle(std::shared_ptr<Oracle<Index_> > input, Index_ shift) : input(std::move(input)), shift(shift) {}

        size_t total() const {
            return input->total();
        }

        Index_ get(size_t i) const {
            return input->get(i) + shift;
        }

    private:
        std::shared_ptr<Oracle<Index_> > input;
        Index_ shift = 0;
    };

    template<bool oracle_, bool sparse_>
    struct AcrossFullExtractor : public GeneralizedExtractor<oracle_, DimensionSelectionType::FULL, sparse_, Value_, Index_> {
        template<class OraclePointer_>
        AcrossFullExtractor(const DelayedSubsetBlock* parent, OraclePointer_ oracle, bool, const Options& opt) {
            this->full_length = (margin_ == 0 ? parent->ncol() : parent->nrow());
            if constexpr(oracle_) {
                internal = new_extractor<margin_ == 0, sparse_>(parent->mat.get(), std::make_shared<SubsetOracle>(std::move(oracle), parent->block_start), opt);
            } else {
                internal = new_extractor<margin_ == 0, sparse_>(parent->mat.get(), opt);
            }
        }

    protected:
        std::unique_ptr<GeneralizedExtractor<oracle_, DimensionSelectionType::FULL, sparse_, Value_, Index_> > internal;
    };

    template<bool oracle_, bool sparse_>
    struct AcrossBlockExtractor : public GeneralizedExtractor<oracle_, DimensionSelectionType::BLOCK, sparse_, Value_, Index_> {
        template<class OraclePointer_>
        AcrossBlockExtractor(const DelayedSubsetBlock* parent, OraclePointer_ oracle, std::pair<Index_, Index_> block, const Options& opt) {
            this->block_start = block.first;
            this->block_length = block.second;
            if constexpr(oracle_) {
                internal = new_extractor<margin_ == 0, sparse_>(parent->mat.get(), std::make_shared<SubsetOracle>(std::move(oracle), parent->block_start), block.first, block.second, opt);
            } else {
                internal = new_extractor<margin_ == 0, sparse_>(parent->mat.get(), block.first, block.second, opt);
            }
        }

    protected:
        std::unique_ptr<GeneralizedExtractor<oracle_, DimensionSelectionType::BLOCK, sparse_, Value_, Index_> > internal;
    };

    template<bool oracle_, bool sparse_>
    struct AcrossIndexExtractor : public GeneralizedExtractor<oracle_, DimensionSelectionType::INDEX, sparse_, Value_, Index_> {
        template<class OraclePointer_>
        AcrossIndexExtractor(const DelayedSubsetBlock* parent, OraclePointer_ oracle, std::vector<Index_> indices, const Options& opt)  {
            this->index_length = indices.size();
            if constexpr(oracle_) {
                internal = new_extractor<margin_ == 0, sparse_>(parent->mat.get(), std::make_shared<SubsetOracle>(std::move(oracle), parent->block_start), std::move(indices), opt);
            } else {
                internal = new_extractor<margin_ == 0, sparse_>(parent->mat.get(), std::move(indices), opt);
            }
        }

    protected:
        std::unique_ptr<GeneralizedExtractor<oracle_, DimensionSelectionType::INDEX, sparse_, Value_, Index_> > internal;

    public:
        const Index_* index_start() const {
            return internal->index_start();
        }
    };

    template<bool oracle_, DimensionSelectionType selection_, bool sparse_>
    using SelectionBasedAcrossExtractor = typename std::conditional<
            selection_ == DimensionSelectionType::FULL,
            AcrossFullExtractor<oracle_, sparse_>,
            typename std::conditional<
                selection_ == DimensionSelectionType::BLOCK,
                AcrossBlockExtractor<oracle_, sparse_>,
                AcrossIndexExtractor<oracle_, sparse_>
            >::type
        >::type;

    template<DimensionSelectionType selection_>
    struct SimpleDenseAcrossExtractor : public SelectionBasedAcrossExtractor<false, selection_, false> {
        template<typename SelectionArg_>
        SimpleDenseAcrossExtractor(const DelayedSubsetBlock* p, bool oracle, SelectionArg_&& sarg, const Options& opt) : 
            SelectionBasedAcrossExtractor<false, selection_, false>(p, oracle, std::forward<SelectionArg_>(sarg), opt), offset(p->block_start) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            return this->internal->fetch(i + this->offset, buffer);
        }

    private:
        Index_ offset = 0;
    };

    template<DimensionSelectionType selection_>
    struct OracleAwareDenseAcrossExtractor : public SelectionBasedAcrossExtractor<true, selection_, false> {
        template<typename SelectionArg_>
        OracleAwareDenseAcrossExtractor(const DelayedSubsetBlock* p, std::shared_ptr<Oracle<Index_> > input_oracle, SelectionArg_&& sarg, const Options& opt) : 
            SelectionBasedAcrossExtractor<true, selection_, false>(p, input_oracle, std::forward<SelectionArg_>(sarg), opt), my_oracle(std::move(input_oracle)) {}

        const Value_* fetch(Index_& i, Value_* buffer) {
            auto output = this->internal->fetch(i, buffer);
            i = my_oracle->get(this->used_predictions);
            ++(this->used_predictions);
            return output;
        }

        const Oracle<Index_>* oracle() const {
            return my_oracle.get();
        }

    private:
        std::shared_ptr<Oracle<Index_> > my_oracle;
    };

    template<DimensionSelectionType selection_>
    struct SimpleSparseAcrossExtractor : public SelectionBasedAcrossExtractor<false, selection_, true> {
        template<typename SelectionArg_>
        SimpleSparseAcrossExtractor(const DelayedSubsetBlock* p, bool oracle, SelectionArg_&& sarg, const Options& opt) :
            SelectionBasedAcrossExtractor<false, selection_, true>(p, oracle, std::forward<SelectionArg_>(sarg), opt), offset(p->block_start) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            return this->internal->fetch(i + this->offset, vbuffer, ibuffer);
        }

    private:
        Index_ offset = 0;
    };

    template<DimensionSelectionType selection_>
    struct OracleAwareSparseAcrossExtractor : public SelectionBasedAcrossExtractor<true, selection_, true> {
        template<typename SelectionArg_>
        OracleAwareSparseAcrossExtractor(const DelayedSubsetBlock* p, std::shared_ptr<Oracle<Index_> > input_oracle, SelectionArg_&& sarg, const Options& opt) :
            SelectionBasedAcrossExtractor<true, selection_, true>(p, input_oracle, std::forward<SelectionArg_>(sarg), opt), my_oracle(std::move(input_oracle)) {}

        SparseRange<Value_, Index_> fetch(Index_& i, Value_* vbuffer, Index_* ibuffer) {
            auto output = this->internal->fetch(i, vbuffer, ibuffer);
            i = my_oracle->get(this->used_predictions);
            ++(this->used_predictions);
            return output;
        }

        const Oracle<Index_>* oracle() const {
            return my_oracle.get();
        }

    private:
        std::shared_ptr<Oracle<Index_> > my_oracle;
    };

    template<bool oracle_, DimensionSelectionType selection_, bool sparse_>
    using AcrossExtractor = typename std::conditional<
            sparse_,
            typename std::conditional<
                    oracle_, 
                    OracleAwareSparseAcrossExtractor<selection_>, 
                    SimpleSparseAcrossExtractor<selection_> 
                >::type,
            typename std::conditional<
                    oracle_, 
                    OracleAwareDenseAcrossExtractor<selection_>, 
                    SimpleDenseAcrossExtractor<selection_> 
                >::type
        >::type;

    /********************************
     ***** Populating extractor *****
     ********************************/
private:
    template<bool accrow_, bool oracle_, DimensionSelectionType selection_, bool sparse_, typename OraclePointer_, typename SelectionArg_>
    std::unique_ptr<GeneralizedExtractor<oracle_, selection_, sparse_, Value_, Index_> > populate(OraclePointer_&& oracle, SelectionArg_&& sarg, const Options& opt) const {
        std::unique_ptr<GeneralizedExtractor<oracle_, selection_, sparse_, Value_, Index_> > output;
        if constexpr(accrow_ != (margin_ == 0)) {
            output.reset(new AlongExtractor<oracle_, selection_, sparse_>(this, std::forward<OraclePointer_>(oracle), std::forward<SelectionArg_>(sarg), opt));
        } else {
            output.reset(new AcrossExtractor<oracle_, selection_, sparse_>(this, std::forward<OraclePointer_>(oracle), std::forward<SelectionArg_>(sarg), opt));
        }
        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return populate<true, false, DimensionSelectionType::FULL, false>(false, false, opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, false, DimensionSelectionType::BLOCK, false>(false, std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, false, DimensionSelectionType::INDEX, false>(false, std::move(indices), opt);
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return populate<false, false, DimensionSelectionType::FULL, false>(false, false, opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, false, DimensionSelectionType::BLOCK, false>(false, std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, false, DimensionSelectionType::INDEX, false>(false, std::move(indices), opt);
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return populate<true, false, DimensionSelectionType::FULL, true>(false, false, opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, false, DimensionSelectionType::BLOCK, true>(false, std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return populate<true, false, DimensionSelectionType::INDEX, true>(false, std::move(indices), opt);
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return populate<false, false, DimensionSelectionType::FULL, true>(false, false, opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, false, DimensionSelectionType::BLOCK, true>(false, std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return populate<false, false, DimensionSelectionType::INDEX, true>(false, std::move(indices), opt);
    }

public:
    std::unique_ptr<FullDenseOracleAwareExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate<true, true, DimensionSelectionType::FULL, false>(std::move(oracle), false, opt);
    }

    std::unique_ptr<BlockDenseOracleAwareExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, true, DimensionSelectionType::BLOCK, false>(std::move(oracle), std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexDenseOracleAwareExtractor<Value_, Index_> > dense_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate<true, true, DimensionSelectionType::INDEX, false>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<FullDenseOracleAwareExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate<false, true, DimensionSelectionType::FULL, false>(std::move(oracle), true, opt);
    }

    std::unique_ptr<BlockDenseOracleAwareExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, true, DimensionSelectionType::BLOCK, false>(std::move(oracle), std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexDenseOracleAwareExtractor<Value_, Index_> > dense_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate<false, true, DimensionSelectionType::INDEX, false>(std::move(oracle), std::move(indices), opt);
    }

public:
    std::unique_ptr<FullSparseOracleAwareExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate<true, true, DimensionSelectionType::FULL, true>(std::move(oracle), true, opt);
    }

    std::unique_ptr<BlockSparseOracleAwareExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<true, true, DimensionSelectionType::BLOCK, true>(std::move(oracle), std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexSparseOracleAwareExtractor<Value_, Index_> > sparse_row(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate<true, true, DimensionSelectionType::INDEX, true>(std::move(oracle), std::move(indices), opt);
    }

    std::unique_ptr<FullSparseOracleAwareExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, const Options& opt) const {
        return populate<false, true, DimensionSelectionType::FULL, true>(std::move(oracle), true, opt);
    }

    std::unique_ptr<BlockSparseOracleAwareExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, Index_ block_start, Index_ block_length, const Options& opt) const {
        return populate<false, true, DimensionSelectionType::BLOCK, true>(std::move(oracle), std::make_pair(block_start, block_length), opt);
    }

    std::unique_ptr<IndexSparseOracleAwareExtractor<Value_, Index_> > sparse_column(std::shared_ptr<Oracle<Index_> > oracle, std::vector<Index_> indices, const Options& opt) const {
        return populate<false, true, DimensionSelectionType::INDEX, true>(std::move(oracle), std::move(indices), opt);
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam margin_ Dimension along which the addition is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 *
 * @param p Pointer to the underlying (pre-subset) `Matrix`.
 * @param f Index of the start of the block. This should be a row index if `margin_ = 0` and a column index otherwise.
 * @param l Index of the one-past-the-end of the block.
 *
 * @return A pointer to a `DelayedSubsetBlock` instance.
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<const Matrix<Value_, Index_> > p, Index_ f, Index_ l) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedSubsetBlock<margin_, Value_, Index_>(std::move(p), f, l));
}

/**
 * @cond
 */
template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedSubsetBlock(std::shared_ptr<Matrix<Value_, Index_> > p, Index_ f, Index_ l) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedSubsetBlock<margin_, Value_, Index_>(std::move(p), f, l));
}
/**
 * @endcond
 */

}

#endif
