#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H

#include <memory>
#include "../../base/Matrix.hpp"
#include "../../base/utils.hpp"

/**
 * @file DelayedUnaryIsometricOp.hpp
 *
 * @brief Delayed unary isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on a single matrix.
 *
 * Implements any operation that preserves the shape of the matrix and operates on each matrix value independently.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `DenseExtractor::fetch()` or friends.
 * We only consider "unary" operations that involve a single `Matrix` - see `DelayedBinaryIsometricOp` for operations between two `Matrix` instances.
 *
 * The `Operation_` class is expected to provide the following static `constexpr` member variables:
 *
 * - `needs_row`: whether the operation needs the row index, e.g., if the operation varies by row.
 * - `needs_column`: whether the operation needs the column index, e.g., if the operation varies by column.
 * - `always_dense`: whether the operation always yields a dense result when applied on sparse data.
 * - `always_sparse`: whether the operation always yields a sparse result when applied on sparse data.
 *   This should not be `true` if `always_dense` is `true`.
 * 
 * The class should implement the following method:
 *
 * - `void dense<row_>(Index_ i, Index_ start, Index_ length, Value_* buffer) const`: 
 *   This method should apply the operation to all values in `buffer`, which contains a contiguous block of elements from row `i` (when `row_ = true`).
 *   This block spans columns from `[start, start + length)`.
 *   If `row_ = false`, `i` is instead a column and the block contains rows.
 * - `void dense<row_>(Index_ i, const Index_* indices, Index_ length, Value_* buffer) const`: 
 *   This method should apply the operation to all values in `buffer`, which contains to a subset of elements from row `i` (when `row_ = true`).
 *   Column indices of the elements in the subset are defined by `indices`, which is a sorted and unique array of length `length`w.
 *   If `row_ = false`, `i` is instead a column and `indices` contains rows.
 * 
 * If neither of `always_dense` or `always_sparse` is `true`, the class should also implement:
 *
 * - `bool actual_sparse() const`: whether this particular instance of the operation yields a sparse result when applied on sparse data.
 *   For example, an addition operation remains sparse if the added value is zero.
 * 
 * If `always_sparse = false`, the class should implement:
 *
 * - `void expanded<row_>(Index_ i, Index_ start, Index_ length, Value_* buffer) const`: 
 *   same as the corresponding `dense()` method, but with opportunities for optimization where it is known that many input values in `buffer` are zero.
 * - `void expanded<row_>(Index_ i, const Index_* indices, Index_ length, Value_* buffer) const`: 
 *   same as the corresponding `dense()` method, but with opportunities for optimization where it is known that many input values in `buffer` are zero.
 *
 * If `always_dense = false`, the class should implement:
 *
 * - `void sparse<row_>(Index_ i, Index_ number, Value_* buffer, const Index_* indices) const`:
 *   This method should apply the operation to all values in `buffer`, which contains `number` elements from row `i` (when `row_ = true`).
 *   The column indices of all elements is contained in `indices` - unless `needs_column = false`, in which case `indices` may be `NULL`.
 *   If `row_ = false`, `i` is instead a column and `indices` contains row indices (or is `NULL`, if `needs_row = false`).
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedUnaryIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying matrix.
     * @param op Instance of the functor class.
     */
    DelayedUnaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) : mat(std::move(p)), operation(std::move(op)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;
    Operation_ operation;
    static_assert(!Operation_::always_dense || !Operation_::always_sparse);

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
        if constexpr(Operation_::always_sparse) {
            return mat->sparse();
        } else if constexpr(!Operation_::always_dense) {
            if (operation.actual_sparse()) {
                return mat->sparse();
            } 
        }
        return false;
    }

    double sparse_proportion() const {
        if constexpr(Operation_::always_sparse) {
            return mat->sparse_proportion();
        } else if constexpr(!Operation_::always_dense) {
            if (operation.actual_sparse()) {
                return mat->sparse_proportion();
            }
        }
        return 0;
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

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<DimensionSelectionType selection_, bool sparse_, bool inner_sparse_ = sparse_>
    struct IsometricExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        IsometricExtractorBase(const DelayedUnaryIsometricOp* p, std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > i) : parent(p), internal(std::move(i)) {
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
        const DelayedUnaryIsometricOp* parent;

        std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > internal;
    };

    /**************************************
     ********** Dense extraction **********
     **************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor : public IsometricExtractorBase<selection_, false> {
        DenseIsometricExtractor(const DelayedUnaryIsometricOp* p, std::unique_ptr<Extractor<selection_, false, Value_, Index_> > i) : 
            IsometricExtractorBase<selection_, false, false>(p, std::move(i)) {}

        const Value_* fetch(Index_ i, Value_* buffer) {
            this->internal->fetch_copy(i, buffer);

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->parent->operation.template dense<accrow_>(i, 0, this->full_length, buffer);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->parent->operation.template dense<accrow_>(i, this->block_start, this->block_length, buffer);
            } else {
                this->parent->operation.template dense<accrow_>(i, this->internal->index_start(), this->index_length, buffer);
            }

            return buffer;
        }
    };

    /***************************************
     ********** Sparse extraction **********
     ***************************************/
private:
    // This extractor doesn't actually require the indices on the extraction
    // dimension e.g., because the operator doesn't rely on position. So we
    // don't have to guarantee the extraction of the indices if the user didn't
    // provide space in the 'ibuffer'.
    template<bool accrow_, DimensionSelectionType selection_> 
    struct SimpleSparseIsometricExtractor : public IsometricExtractorBase<selection_, true> {
        SimpleSparseIsometricExtractor(const DelayedUnaryIsometricOp* p, std::unique_ptr<Extractor<selection_, true, Value_, Index_> > i) : 
            IsometricExtractorBase<selection_, true, true>(p, std::move(i)) {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto raw = this->internal->fetch(i, vbuffer, ibuffer);

            if (raw.value) {
                if (raw.value != vbuffer) {
                    // Don't use fetch_copy as we don't want to copy the indices.
                    std::copy(raw.value, raw.value + raw.number, vbuffer);
                }

                if constexpr(!Operation_::always_dense) { // avoid the need for operations to define sparse if they're always dense.
                    this->parent->operation.template sparse<accrow_>(i, raw.number, vbuffer, raw.index);
                }
                raw.value = vbuffer;
            }

            return raw;
        }
    };

    template<bool accrow_, DimensionSelectionType selection_> 
    struct RegularSparseIsometricExtractor : public IsometricExtractorBase<selection_, true> {
        RegularSparseIsometricExtractor(const DelayedUnaryIsometricOp* p, std::unique_ptr<Extractor<selection_, true, Value_, Index_> > i, bool ri, bool report_value) : 
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
                if (raw.value != vbuffer) {
                    // Don't use fetch_copy as we don't want to copy the indices.
                    std::copy(raw.value, raw.value + raw.number, vbuffer);
                }

                if constexpr(!Operation_::always_dense) { // avoid the need for operations to define sparse if they're always dense.
                    this->parent->operation.template sparse<accrow_>(i, raw.number, vbuffer, raw.index);
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
    // so we won't bother doing that.
    template<bool accrow_, DimensionSelectionType selection_>
    struct DensifiedSparseIsometricExtractor : public IsometricExtractorBase<selection_, true, false> {
        DensifiedSparseIsometricExtractor(const DelayedUnaryIsometricOp* p, std::unique_ptr<Extractor<selection_, false, Value_, Index_> > i, bool ri, bool rv) :
            IsometricExtractorBase<selection_, true, false>(p, std::move(i)), report_index(ri), report_value(rv)
        {}

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            SparseRange<Value_, Index_> output(extracted_length<selection_, Index_>(*this), NULL, NULL);

            if (report_value) {
                this->internal->fetch_copy(i, vbuffer);

                // 'expanded' is a special case of 'dense' where we assume that most elements are zero.
                // This provides some opportunities to precompute the zero transform, if such a value exists.
                // However, if always_sparse is true, operations can omit the definition of expanded. 
                if constexpr(!Operation_::always_sparse) {
                    if constexpr(selection_ == DimensionSelectionType::FULL) {
                        this->parent->operation.template expanded<accrow_>(i, 0, this->full_length, vbuffer);
                    } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                        this->parent->operation.template expanded<accrow_>(i, this->block_start, this->block_length, vbuffer);
                    } else {
                        this->parent->operation.template expanded<accrow_>(i, this->internal->index_start(), this->index_length, vbuffer);
                    }
                }

                output.value = vbuffer;
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

                output.index = ibuffer;
            }

            return output;
        }

    protected:
        bool report_index = false;
        bool report_value = false;
    };

    /**********************************************
     ********** Public extractor methods **********
     **********************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, typename ... Args_>
    void fill_densified(std::unique_ptr<Extractor<selection_, true, Value_, Index_> >& output, const Options& opt, Args_&& ... args) const {
        bool report_value = opt.sparse_extract_value;
        bool report_index = opt.sparse_extract_index;
        auto inner = new_extractor<accrow_, false>(mat.get(), std::forward<Args_>(args)..., opt);
        output.reset(new DensifiedSparseIsometricExtractor<accrow_, selection_>(this, std::move(inner), report_index, report_value));
    }

    template<bool accrow_, DimensionSelectionType selection_, typename ... Args_>
    void fill_sparse(std::unique_ptr<Extractor<selection_, true, Value_, Index_> >& output, const Options& opt, Args_&& ... args) const {
        if constexpr((accrow_ && !Operation_::needs_column) || (!accrow_ && !Operation_::needs_row)) { // i.e., extraction indices aren't required for the operation.
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
    }

    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > propagate(const Options& opt, Args_&& ... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(!sparse_) {
            auto inner = new_extractor<accrow_, false>(mat.get(), std::forward<Args_>(args)..., opt);
            output.reset(new DenseIsometricExtractor<accrow_, selection_>(this, std::move(inner)));
        } else if constexpr(Operation_::always_sparse) {
            fill_sparse<accrow_, selection_>(output, opt, std::forward<Args_>(args)...);
        } else if constexpr(Operation_::always_dense) {
            fill_densified<accrow_, selection_>(output, opt, std::forward<Args_>(args)...);
        } else if (operation.actual_sparse()) {
            fill_sparse<accrow_, selection_>(output, opt, std::forward<Args_>(args)...);
        } else {
            fill_densified<accrow_, selection_>(output, opt, std::forward<Args_>(args)...);
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
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Helper class defining the operation.
 *
 * @param p Pointer to a (possibly `const`) `Matrix`.
 * @param op Instance of the operation helper class.
 *
 * @return Instance of a `DelayedUnaryIsometricOp` clas.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedUnaryIsometricOp(std::shared_ptr<Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedUnaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

// For back-compatibility.
template<typename ... Args_>
auto make_DelayedUnaryIsometricOp(Args_&&... args) {
    return make_DelayedUnaryIsometricOp(std::forward<Args_>(args)...);
}

template<typename Value_, typename Index_, class Operation_>
using DelayedIsometricOp = DelayedUnaryIsometricOp<Value_, Index_, Operation_>;
/**
 * @endcond
 */

}

#include "arith_helpers.hpp"

#include "math_helpers.hpp"

#include "compare_helpers.hpp"

#include "boolean_helpers.hpp"

#endif
