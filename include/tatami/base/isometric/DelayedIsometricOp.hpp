#ifndef TATAMI_DELAYED_ISOMETRIC_OP_H
#define TATAMI_DELAYED_ISOMETRIC_OP_H

#include <memory>
#include "../Matrix.hpp"
#include "../utils.hpp"

/**
 * @file DelayedIsometricOp.hpp
 *
 * @brief Delayed isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on a matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Functor class implementing the operation.
 * This should accept the row index, column index and value, and return the modified value after applying the operation. 
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) : mat(p), operation(std::move(op)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    Operation_ operation;
    static_assert(std::is_same<Value_, decltype(operation(0, 0, 0))>::value);

public:
    Index_ nrow() const {
        return mat->nrow();
    }
    
    Index_ ncol() const {
        return mat->ncol();
    }

    /**
     * @return `true` if both the underlying (pre-operation) matrix is sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(Operation_::sparse_) {
            return mat->sparse();
        } else {
            return false;
        }
    }

    bool prefer_rows() const { 
        return mat->prefer_rows();
    }

    bool uses_oracle(bool row) const {
        return mat->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<DimensionSelectionType selection_, bool sparse_, bool inner_sparse_>
    struct IsometricExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        IsometricExtractorBase(const DelayedIsometricOp* p, std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > i) : parent(p), internal(std::move(i)) {
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

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            internal->set_oracle(std::move(o));
        }

    protected:
        const DelayedIsometricOp* parent;

        std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > internal;

        friend class DelayedIsometricOp;
    };

    /**************************************
     ********** Dense extraction **********
     **************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor : public IsometricExtractorBase<selection_, false, false> {
        DenseIsometricExtractor(const DelayedIsometricOp* p, std::unique_ptr<Extractor<selection_, false, Value_, Index_> > i) : 
            IsometricExtractorBase<selection_, false, false>(p, std::move(i)) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ptr = this->internal->fetch(i, buffer);
            this->parent->template mutate_dense<accrow_, selection_>(ptr, buffer, i, this);
            return buffer;
        }
    };

    template<bool accrow_, DimensionSelectionType selection_, class IsometricExtractor_>
    void mutate_dense(const Value_* ptr, Value_* buffer, Index_ i, const IsometricExtractor_* exptr) const {
        if constexpr(selection_ == DimensionSelectionType::FULL) {
            auto end = exptr->full_length;
            for (Index_ j = 0; j < end; ++j, ++ptr) {
                if constexpr(accrow_) {
                    buffer[j] = operation(i, j, *ptr);
                } else {
                    buffer[j] = operation(j, i, *ptr);
                }
            }

        } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
            auto end = exptr->block_length;
            auto shift = exptr->block_start;
            for (Index_ j = 0; j < end; ++j, ++ptr) {
                if constexpr(accrow_) {
                    buffer[j] = operation(i, j + shift, *ptr);
                } else {
                    buffer[j] = operation(j + shift, i, *ptr);
                }
            }

        } else {
            auto xptr = exptr->internal->index_start();
            auto end = exptr->index_length;
            for (Index_ j = 0; j < end; ++j, ++ptr, ++xptr) {
                if constexpr(accrow_) {
                    buffer[j] = operation(i, *xptr, *ptr);
                } else {
                    buffer[j] = operation(*xptr, i, *ptr);
                }
            }
        }
    }

    /***************************************
     ********** Sparse extraction **********
     ***************************************/
private:
    // This extractor doesn't actually require the indices on the extraction
    // dimension e.g., because the operator doesn't rely on position. So we
    // don't have to guarantee the extraction of the indices if the user didn't
    // provide space in the 'ibuffer'.
    template<bool accrow_, DimensionSelectionType selection_> 
    struct SimpleSparseIsometricExtractor : public IsometricExtractorBase<selection_, true, true> {
        SimpleSparseIsometricExtractor(const DelayedIsometricOp* p, std::unique_ptr<Extractor<selection_, true, Value_, Index_> > i) : 
            IsometricExtractorBase<selection_, true, true>(p, std::move(i)) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto raw = this->internal->fetch(i, vbuffer, ibuffer);

            if (raw.value) {
                const auto& OP = this->parent->operation;
                for (Index_ j = 0; j < raw.number; ++j) {
                    if constexpr(accrow_) {
                        vbuffer[j] = OP(i, 0, raw.value[j]); // no-op value, we don't need the column indices.
                    } else {
                        vbuffer[j] = OP(0, i, raw.value[j]); // no-op value, we don't need the row indices.
                    }
                }
                raw.value = vbuffer;
            }

            return raw;
        }
    };

    template<bool accrow_, DimensionSelectionType selection_> 
    struct RegularSparseIsometricExtractor : public IsometricExtractorBase<selection_, true, true> {
        RegularSparseIsometricExtractor(const DelayedIsometricOp* p, std::unique_ptr<Extractor<selection_, true, Value_, Index_> > i, bool ri, bool report_value) : 
            IsometricExtractorBase<selection_, true, true>(p, std::move(i)), report_index(ri)
        {
            if (!report_index && report_value) {
                // We only need an internal ibuffer if the user wants the
                // values but didn't provide enough space to store the indices
                // (which we need to pass to the operation's functor).
                internal_ibuffer.resize(extracted_length<selection_, Index_>(*this));
            }
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            // If the internal ibuffer is empty, we're either extracting the indices
            // directly into the user's ibuffer, or we don't need the indices at all.
            // Either way, it doesn't hurt to use the user's ibuffer.
            Index_* iin = (internal_ibuffer.empty() ? ibuffer : internal_ibuffer.data());

            auto raw = this->internal->fetch(i, vbuffer, iin);

            if (raw.value) {
                const auto& OP = this->parent->operation;
                for (Index_ j = 0; j < raw.number; ++j) {
                    if constexpr(accrow_) {
                        vbuffer[j] = OP(i, raw.index[j], raw.value[j]);
                    } else {
                        vbuffer[j] = OP(raw.index[j], i, raw.value[j]);
                    }
                }
                raw.value = vbuffer;
            }

            if (!report_index) {
                raw.index = NULL;
            }

            return raw;
        }

    protected:
        std::vector<Index_> internal_ibuffer;
        bool report_index = false;
    };

    /*******************************************
     ********** Un-sparsed extraction **********
     *******************************************/
private:
    // Technically, we could avoid constructing the internal extractor if
    // we don't want the values, but that's a pretty niche optimization,
    // so we wont' bother doing that.
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DensifiedSparseIsometricExtractor : public IsometricExtractorBase<selection_, true, false> {
        DensifiedSparseIsometricExtractor(const DelayedIsometricOp* p, std::unique_ptr<Extractor<selection_, false, Value_, Index_> > i, bool ri, bool rv) : 
            IsometricExtractorBase<selection_, true, false>(p, std::move(i)), report_index(ri), report_value(rv)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            SparseRange<Value_, Index_> raw(extracted_length<selection_, Index_>(*this), NULL, NULL);

            if (report_value) {
                auto ptr = this->internal->fetch(i, vbuffer);
                this->parent->template mutate_dense<accrow_, selection_>(ptr, vbuffer, i, this);
                raw.value = vbuffer;
            }

            if (report_index) {
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    std::iota(ibuffer, ibuffer + this->full_length, 0);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    std::iota(ibuffer, ibuffer + this->block_length, this->block_start);
                } else {
                    auto xptr = this->internal->index_start();
                    std::copy(xptr, xptr + this->index_length, ibuffer);
                }
                raw.index = ibuffer;
            }

            return raw;
        }

    protected:
        bool report_index = false;
        bool report_value = false;
    };

    /**********************************************
     ********** Public extractor methods **********
     **********************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > propagate(const Options& opt, Args_&& ... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(!sparse_) {
            auto inner = new_extractor<accrow_, false>(mat.get(), std::forward<Args_>(args)..., opt);
            output.reset(new DenseIsometricExtractor<accrow_, selection_>(this, std::move(inner)));

        } else if constexpr(!Operation_::sparse_) {
            bool report_value = opt.sparse_extract_value;
            bool report_index = opt.sparse_extract_index;
            auto inner = new_extractor<accrow_, false>(mat.get(), std::forward<Args_>(args)..., opt);
            output.reset(new DensifiedSparseIsometricExtractor<accrow_, selection_>(this, std::move(inner), report_index, report_value));

        } else if constexpr((accrow_ && !Operation_::needs_column_) || (!accrow_ && !Operation_::needs_row_)) { // i.e., extraction indices aren't required for the operation.
            auto inner = new_extractor<accrow_, true>(mat.get(), std::forward<Args_>(args)..., opt);
            output.reset(new SimpleSparseIsometricExtractor<accrow_, selection_>(this, std::move(inner)));

        } else {
            bool report_value = opt.sparse_extract_value;
            bool report_index = opt.sparse_extract_index;

            std::unique_ptr<Extractor<selection_, true, Value_, Index_> > inner;
            if (report_value && !report_index) {
                auto optcopy = opt;
                optcopy.sparse_extract_index = true; // if we get to this clause, we need the indices; otherwise we'd have constructed a SimpleSparseIsometricExtractor.
                inner = new_extractor<accrow_, true>(mat.get(), std::forward<Args_>(args)..., optcopy);
            } else {
                inner = new_extractor<accrow_, true>(mat.get(), std::forward<Args_>(args)..., opt);
            }

            output.reset(new RegularSparseIsometricExtractor<accrow_, selection_>(this, std::move(inner), report_index, report_value));
        }

        return output;
    }

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options& opt) const {
        return propagate<true, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<true, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options& opt) const {
        return propagate<true, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options& opt) const {
        return propagate<false, DimensionSelectionType::FULL, false>(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<false, DimensionSelectionType::BLOCK, false>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options& opt) const {
        return propagate<false, DimensionSelectionType::INDEX, false>(opt, std::move(indices));
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options& opt) const {
        return propagate<true, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<true, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options& opt) const {
        return propagate<true, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options& opt) const {
        return propagate<false, DimensionSelectionType::FULL, true>(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options& opt) const {
        return propagate<false, DimensionSelectionType::BLOCK, true>(opt, block_start, block_length);
    }
    
    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options& opt) const {
        return propagate<false, DimensionSelectionType::INDEX, true>(opt, std::move(indices));
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam Matrix_ A realized `Matrix` class, possibly one that is `const`.
 * @tparam Operation_ Helper class defining the operation.
 *
 * @param p Pointer to a `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedIsometricOp` clas.
 */
template<class Matrix_, class Operation_>
std::shared_ptr<Matrix<typename Matrix_::value_type, typename Matrix_::index_type> > make_DelayedIsometricOp(std::shared_ptr<Matrix_> p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix_>(new DelayedIsometricOp<typename Matrix_::value_type, typename Matrix_::index_type, Op_>(p, std::move(op)));
}

}

#include "arith_scalar_helpers.hpp"

#include "arith_vector_helpers.hpp"

#include "math_helpers.hpp"

#endif
