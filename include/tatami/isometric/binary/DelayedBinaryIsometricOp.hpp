#ifndef TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H
#define TATAMI_DELAYED_UNARY_ISOMETRIC_OP_H

#include <memory>
#include "../../base/Matrix.hpp"
#include "../../base/utils.hpp"

/**
 * @file DelayedBinaryIsometricOp.hpp
 *
 * @brief Delayed binary isometric operations.
 *
 * This is equivalent to the class of the same name in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed isometric operations on two matrices
 *
 * Implements any operation that takes two matrices of the same shape and returns another matrix of that shape.
 * Each entry of the output matrix is a function of the corresponding values in the two input matrices.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `DenseExtractor::fetch()` or friends.
 *
 * The `Operation_` class is expected to provide the following static `constexpr` member variables:
 *
 * - `sparse`: whether the operation can be optimized if both input matrices are sparse.
 * 
 * The class should implement the following method:
 *
 * - `void dense<row_>(Index_ i, Index_ start, Index_ length, Value_* left_buffer, const Value_* right_buffer) const`: 
 *   This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
 *   each of which contain a contiguous block of elements from row `i` of the left and right matrices, respectively (when `row_ = true`).
 *   The result of the operation should be stored in `left_buffer`.
 *   The block starts at column `start` and is of length `length`.
 *   If `row_ = false`, `i` is instead a column and the block starts at row `start`.
 * - `void dense<row_>(Index_ i, const Index_* indices, Index_ length, Value_* buffer1, const Value_* buffer2) const`: 
 *   This method should apply the operation to corresponding values of `left_buffer` and `right_buffer`,
 *   each of which contain a subset of elements from row `i` of the left and right matrices, respectively (when `row_ = true`).
 *   The result of the operation should be stored in `left_buffer`.
 *   The subset is defined by column indices in the `indices` array of length `length`.
 *   If `row_ = false`, `i` is instead a column and `indices` contains rows.
 * 
 * If `sparse = true`, the class should implement:
 *
 * - `Index_ sparse<row_, needs_value, needs_index>(Index_ i, const SparseRange<Value_, Index_>& left, const SparseRange<Value_, Index_>& right, Value_* value_buffer, Index_* index_buffer) const`:
 *   This method should apply the operation to the sparse values in `left` and `right`, 
 *   consisting of the contents of row `i` from the left and right matrices, respectively (when `row_ = true`).
 *   All non-zero values resulting from the operation should be stored in `value_buffer` if `needs_value = true`, otherwise `value_buffer = NULL` and should be ignored.
 *   The corresponding indices of those values should be stored in `index_buffer` if `needs_index = true`, otherwise `index_buffer = NULL` and should be ignored.
 *   The return value should be the number of structural non-zero elements in the output buffers.
 *   If `row_ = false`, the contents of `left` and `right` are taken from column `i` instead.
 *   Note that all values in `left` and `right` are already sorted by increasing index.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 * @tparam Operation_ Class implementing the operation.
 */
template<typename Value_, typename Index_, class Operation_>
class DelayedBinaryIsometricOp : public Matrix<Value_, Index_> {
public:
    /**
     * @param l Pointer to the left matrix.
     * @param r Pointer to the right matrix.
     * @param op Instance of the functor class.
     */
    DelayedBinaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > l, std::shared_ptr<const Matrix<Value_, Index_> > r, Operation_ op) : 
        left(std::move(l)), right(std::move(r)), operation(std::move(op)) 
    {
        if (left->nrow() != right->nrow() || left->ncol() != right->ncol()) {
            throw std::runtime_error("shape of the left and right matrices should be the same");
        }

        prefer_rows_proportion_internal = (left->prefer_rows_proportion() + right->prefer_rows_proportion()) / 2;
    }

private:
    std::shared_ptr<const Matrix<Value_, Index_> > left, right;
    Operation_ operation;
    double prefer_rows_proportion_internal;

public:
    Index_ nrow() const {
        return mat->nrow();
    }

    Index_ ncol() const {
        return mat->ncol();
    }

    /**
     * @return `true` if both underlying (pre-operation) matrices are sparse and the operation preserves sparsity.
     * Otherwise returns `false`.
     */
    bool sparse() const {
        if constexpr(Operation_::sparse) {
            return left->sparse() && right->sparse();
        }
        return false;
    }

    double sparse_proportion() const {
        if constexpr(Operation_::sparse) {
            // Well, better than nothing.
            return (left->sparse_proportion() + right->sparse_proportion())/2;
        }
        return 0;
    }

    bool prefer_rows() const { 
        return prefer_rows_proportion_internal > 0.5;
    }

    double prefer_rows_proportion() const { 
        return prefer_rows_proportion_internal;
    }

    bool uses_oracle(bool row) const {
        // TODO: insert oracle.
        return false;
        //return left->uses_oracle(row) || right->uses_oracle(row);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

private:
    template<DimensionSelectionType selection_, bool sparse_, bool inner_sparse_ = sparse_>
    struct IsometricExtractorBase : public Extractor<selection_, sparse_, Value_, Index_> {
        IsometricExtractorBase(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > l,
            std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > r
        ) : 
            parent(p), 
            left_internal(std::move(l)),
            right_internal(std::move(r))
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

        void set_oracle(std::unique_ptr<Oracle<Index_> > o) {
            // Need to propagate the oracle!
        }

    protected:
        const DelayedBinaryIsometricOp* parent;

        std::unique_ptr<Extractor<selection_, inner_sparse_, Value_, Index_> > left_internal, right_internal;
    };

    /**************************************
     ********** Dense extraction **********
     **************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct DenseIsometricExtractor : public IsometricExtractorBase<selection_, false> {
        DenseIsometricExtractor(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > l, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > r 
        ) : 
            IsometricExtractorBase<selection_, false, false>(p, std::move(l), std::move(r))
        {
            holding_buffer.resize(extracted_length<selection_, Index_>(this));
        }

        const Value_* fetch(Index_ i, Value_* buffer) {
            this->left_internal->fetch_copy(i, buffer);
            auto rptr = this->right_internal->fetch(i, holding_buffer.data());

            if constexpr(selection_ == DimensionSelectionType::FULL) {
                this->parent->operation.template dense<accrow_>(i, 0, this->full_length, buffer, rptr);
            } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                this->parent->operation.template dense<accrow_>(i, this->block_start, this->block_length, buffer, rptr);
            } else {
                this->parent->operation.template dense<accrow_>(i, this->internal->index_start(), this->index_length, buffer, rptr);
            }

            return buffer;
        }

    private:
        std::vector<Value_> holding_buffer;
    };

    /***************************************
     ********** Sparse extraction **********
     ***************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_> 
    struct RegularSparseIsometricExtractor : public IsometricExtractorBase<selection_, true> {
        RegularSparseIsometricExtractor(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, true, Value_, Index_> > l, 
            std::unique_ptr<Extractor<selection_, true, Value_, Index_> > r, 
            bool rv,
            bool ri
        ) : 
            IsometricExtractorBase<selection_, true, true>(p, std::move(l), std::move(r)), 
            report_value(rv),
            report_index(ri)
        {
            auto n = extracted_length<selection_, Index_>(*this):
            left_internal_ibuffer.resize(n);
            right_internal_ibuffer.resize(n);

            if (report_value) {
                left_internal_vbuffer.resize(n);
                right_internal_vbuffer.resize(n);
            }
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto left_ranges = this->left_internal->fetch(i, left_internal_vbuffer.data(), left_internal_ibuffer.data());
            auto right_ranges = this->right_internal->fetch(i, right_internal_vbuffer.data(), right_internal_ibuffer.data());

            SparseRange<Value_, Index_> output(0, NULL, NULL);

            if (report_value && report_index) {
                output.number = operation.sparse<accrow_, true, true>(i, left_ranges, right_ranges, vbuffer, ibuffer);
                output.value = vbuffer;
                output.index = ibuffer;
            } else if (report_value) {
                output.number = operation.sparse<accrow_, true, false>(i, left_ranges, right_ranges, vbuffer, NULL);
                output.value = vbuffer;
            } else if (report_index) {
                output.number = operation.sparse<accrow_, false, true>(i, left_ranges, right_ranges, NULL, ibuffer);
                output.index = ibuffer;
            } else {
                output.number = operation.sparse<accrow_, false, false>(i, left_ranges, right_ranges, NULL, NULL);
            }

            return raw;
        }

    protected:
        std::vector<Value_> left_internal_vbuffer, right_internal_vbuffer;
        std::vector<Index_> left_internal_ibuffer, right_internal_ibuffer;
        bool report_value = false;
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
        DensifiedSparseIsometricExtractor(
            const DelayedBinaryIsometricOp* p, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > l, 
            std::unique_ptr<Extractor<selection_, false, Value_, Index_> > r 
            bool rv,
            bool ri
        ) :
            IsometricExtractorBase<selection_, true, false>(p, std::move(l), std::move(r)), 
            report_value(rv),
            report_index(ri) 
        {
            holding_buffer.resize(extracted_length<selection_, Index_>(this));
        }

        SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            SparseRange<Value_, Index_> output(extracted_length<selection_, Index_>(*this), NULL, NULL);

            if (report_value) {
                this->left_internal->fetch_copy(i, vbuffer);
                auto rptr = this->right_internal->fetch(i, holding_buffer.data());

                if constexpr(!Operation_::always_sparse) {
                    if constexpr(selection_ == DimensionSelectionType::FULL) {
                        this->parent->operation.template expanded<accrow_>(i, 0, this->full_length, vbuffer, rptr);
                    } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
                        this->parent->operation.template expanded<accrow_>(i, this->block_start, this->block_length, vbuffer, rptr);
                    } else {
                        this->parent->operation.template expanded<accrow_>(i, this->internal->index_start(), this->index_length, vbuffer, rptr);
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
        std::vector<Value_> holding_buffer;
        bool report_value = false;
        bool report_index = false;
    };

    /**********************************************
     ********** Public extractor methods **********
     **********************************************/
private:
    template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename ... Args_>
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > propagate(const Options& opt, Args_&& ... args) const {
        std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

        if constexpr(!sparse_) {
            auto left_inner = new_extractor<accrow_, false>(left.get(), std::forward<Args_>(args)..., opt);
            auto right_inner = new_extractor<accrow_, false>(right.get(), std::forward<Args_>(args)..., opt);
            output.reset(new DenseIsometricExtractor<accrow_, selection_>(this, std::move(left_inner), std::move(right_inner)));

        } else if constexpr(Operation_::sparse) {
            bool report_value = opt.sparse_extract_value;
            bool report_index = opt.sparse_extract_index;

            auto optcopy = opt;
            optcopy.sparse_extract_index = true; // We need the indices to combine things properly.
            optcopy.sparse_ordered_index = true; // Make life easier for operation implementers.

            auto left_inner = new_extractor<accrow_, true>(left.get(), std::forward<Args_>(args)..., optcopy);
            auto right_inner = new_extractor<accrow_, true>(right.get(), std::forward<Args_>(args)..., optcopy);
            output.reset(new RegularSparseIsometricExtractor<accrow_, selection_>(this, std::move(left_inner), std::move(right_inner), report_value, report_index));

        } else {
            bool report_value = opt.sparse_extract_value;
            bool report_index = opt.sparse_extract_index;
            auto left_inner = new_extractor<accrow_, false>(left.get(), std::forward<Args_>(args)..., opt);
            auto right_inner = new_extractor<accrow_, false>(right.get(), std::forward<Args_>(args)..., opt);
            output.reset(new DensifiedSparseIsometricExtractor<accrow_, selection_>(this, std::move(left_inner), std::move(right_inner), report_value, report_index));
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
 * @return Instance of a `DelayedBinaryIsometricOp` clas.
 */
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBinaryIsometricOp(std::shared_ptr<const Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBinaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

/**
 * @cond
 */
// For automatic template deduction with non-const pointers.
template<typename Value_, typename Index_, class Operation_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBinaryIsometricOp(std::shared_ptr<Matrix<Value_, Index_> > p, Operation_ op) {
    typedef typename std::remove_reference<Operation_>::type Op_;
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBinaryIsometricOp<Value_, Index_, Op_>(std::move(p), std::move(op)));
}

// For back-compatibility.
template<typename ... Args_>
auto make_DelayedBinaryIsometricOp(Args_&&... args) {
    return make_DelayedBinaryIsometricOp(std::forward<Args_>(args)...);
}

template<typename Value_, typename Index_, class Operation_>
using DelayedIsometricOp = DelayedBinaryIsometricOp<Value_, Index_, Operation_>;
/**
 * @endcond
 */

}

#include "arith_helpers.hpp"

#endif
