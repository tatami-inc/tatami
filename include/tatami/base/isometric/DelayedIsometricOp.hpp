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

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<DimensionSelectionType selection_, bool sparse_, bool inner_sparse_>
    struct IsometricExtractorBase : public Extractor<sparse_, Value_, Index_> {
        IsometricExtractorBase(std::unique_ptr<Extractor<inner_sparse_, Value_, Index_> > i, const DelayedIsometricOp* p) : internal(std::move(i)), parent(p) {
            this->extracted_selection = internal->extracted_selection;
            this->extracted_length = internal->extracted_length;
            if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->extracted_block = internal->extracted_block;
            }
        }

        const Index_* extracted_index() const {
            return internal->extracted_index();
        }

    protected:
        std::unique_ptr<Extractor<inner_sparse_, Value_, Index_> > internal;

        const DelayedIsometricOp* parent;

        friend class DelayedIsometricOp;
    };

    /**************************************
     ********** Dense extraction **********
     **************************************/
public:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor : public IsometricExtractorBase<selection_, false, false> {
        DenseIsometricExtractor(std::unique_ptr<Extractor<false, Value_, Index_> > i, const DelayedIsometricOp* p) : 
            IsometricExtractorBase<selection_, false, false>(std::move(i), p) 
        {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            auto ptr = this->internal->fetch(i, buffer);
            this->parent->template mutate_dense<accrow_, selection_>(ptr, buffer, i, this);
            return buffer;
        }
    };

    template<bool accrow_, DimensionSelectionType selection_, class IsometricExtractor_>
    void mutate_dense(const Value_* ptr, Value_* buffer, Index_ i, const IsometricExtractor_* exptr) const {
        auto end = exptr->extracted_length;

        if constexpr(selection_ == DimensionSelectionType::FULL) {
            for (Index_ j = 0; j < end; ++j, ++ptr) {
                if constexpr(accrow_) {
                    buffer[j] = operation(i, j, *ptr);
                } else {
                    buffer[j] = operation(j, i, *ptr);
                }
            }

        } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
            auto shift = exptr->extracted_block;
            for (Index_ j = 0; j < end; ++j, ++ptr) {
                if constexpr(accrow_) {
                    buffer[j] = operation(i, j + shift, *ptr);
                } else {
                    buffer[j] = operation(j + shift, i, *ptr);
                }
            }

        } else {
            auto xptr = exptr->internal->extracted_index();
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
        SimpleSparseIsometricExtractor(std::unique_ptr<Extractor<true, Value_, Index_> > i, const DelayedIsometricOp* p) : 
            IsometricExtractorBase<selection_, true, true>(std::move(i), p)
        {}

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
        RegularSparseIsometricExtractor(std::unique_ptr<Extractor<true, Value_, Index_> > i, const DelayedIsometricOp* p, bool ri, bool report_value) : 
            IsometricExtractorBase<selection_, true, true>(std::move(i), p), report_index(ri)
        {
            if (!report_index && report_value) {
                // We only need an internal ibuffer if the user wants the
                // values but didn't provide enough space to store the indices
                // (which we need to pass to the operation's functor).
                internal_ibuffer.resize(this->extracted_length);
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
        bool report_index = false;
        std::vector<Index_> internal_ibuffer;
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
        DensifiedSparseIsometricExtractor(std::unique_ptr<Extractor<false, Value_, Index_> > i, const DelayedIsometricOp* p, bool ri, bool rv) : 
            IsometricExtractorBase<selection_, true, false>(std::move(i), p), report_index(ri), report_value(rv)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            SparseRange<Value_, Index_> raw(this->extracted_length, NULL, NULL);

            if (report_value) {
                auto ptr = this->internal->fetch(i, vbuffer);
                this->parent->template mutate_dense<accrow_, selection_>(ptr, vbuffer, i, this);
                raw.value = vbuffer;
            }

            if (report_index) {
                if constexpr(selection_ == DimensionSelectionType::FULL) {
                    std::iota(ibuffer, ibuffer + this->extracted_length, 0);
                } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                    std::iota(ibuffer, ibuffer + this->extracted_length, this->extracted_block);
                } else {
                    auto xptr = this->internal->extracted_index();
                    std::copy(xptr, xptr + this->extracted_length, ibuffer);
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
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_>
    std::unique_ptr<Extractor<sparse_, Value_, Index_> > propagate_by_operation(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        std::unique_ptr<Extractor<sparse_, Value_, Index_> > output;

        if constexpr(!sparse_) {
            auto inner = new_extractor<accrow_, false>(mat.get(), std::move(iopt), std::move(eopt));
            output.reset(new DenseIsometricExtractor<accrow_, selection_>(std::move(inner), this));

        } else if constexpr(!Operation_::sparse_) {
            // Make copies to avoid invalidation on move(eopt).
            bool report_value = eopt.sparse_extract_value;
            bool report_index = eopt.sparse_extract_index;

            auto inner = new_extractor<accrow_, false>(mat.get(), std::move(iopt), std::move(eopt));
            output.reset(new DensifiedSparseIsometricExtractor<accrow_, selection_>(std::move(inner), this, report_index, report_value));

        } else if constexpr((accrow_ && !Operation_::needs_column_) || (!accrow_ && !Operation_::needs_row_)) { // i.e., extraction indices aren't required for the operation.
            auto inner = new_extractor<accrow_, true>(mat.get(), std::move(iopt), std::move(eopt));
            output.reset(new SimpleSparseIsometricExtractor<accrow_, selection_>(std::move(inner), this));

        } else {
            bool report_value = eopt.sparse_extract_value;
            bool report_index = eopt.sparse_extract_index;
            if (report_value && !report_index) {
                eopt.sparse_extract_index = true; // if we get to this clause, we need the indices; otherwise we'd have constructed a SimpleSparseIsometricExtractor.
            }

            auto inner = new_extractor<accrow_, true>(mat.get(), std::move(iopt), std::move(eopt));
            output.reset(new RegularSparseIsometricExtractor<accrow_, selection_>(std::move(inner), this, report_index, report_value));
        }

        return output;
    }

    template<bool accrow_, bool sparse_>
    std::unique_ptr<Extractor<sparse_, Value_, Index_> > propagate(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        switch (eopt.selection.type) {
            case DimensionSelectionType::FULL:
                return propagate_by_operation<accrow_, DimensionSelectionType::FULL, sparse_>(std::move(iopt), std::move(eopt));
            case DimensionSelectionType::BLOCK:
                return propagate_by_operation<accrow_, DimensionSelectionType::BLOCK, sparse_>(std::move(iopt), std::move(eopt));
            case DimensionSelectionType::INDEX:
                return propagate_by_operation<accrow_, DimensionSelectionType::INDEX, sparse_>(std::move(iopt), std::move(eopt));
        }
    }

public:
    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_row(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        return propagate<true, false>(std::move(iopt), std::move(eopt));
    }

    std::unique_ptr<DenseExtractor<Value_, Index_> > dense_column(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        return propagate<false, false>(std::move(iopt), std::move(eopt));
    }

    std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_row(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        return propagate<true, true>(std::move(iopt), std::move(eopt));
    }

    std::unique_ptr<SparseExtractor<Value_, Index_> > sparse_column(IterationOptions<Index_> iopt, ExtractionOptions<Index_> eopt) const {
        return propagate<false, true>(std::move(iopt), std::move(eopt));
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
