#ifndef TATAMI_DELAYED_SUBSET_BLOCK
#define TATAMI_DELAYED_SUBSET_BLOCK

#include "typed_matrix.hpp"
#include <algorithm>
#include <memory>

/**
 * @file DelayedSubsetBlock.hpp
 *
 * Delayed subsetting to a single contiguous block.
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
 * @tparam MARGIN Dimension along which the addition is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam T Type of matrix value.
 * @tparam IDX Type of index value.
 */
template<int MARGIN, typename T, typename IDX>
class DelayedSubsetBlock : public typed_matrix<T, IDX> {
public:
    /**
     * @param p Pointer to the underlying (pre-subset) matrix.
     * @param f Index of the start of the block. This should be a row index if `MARGIN = 0` and a column index otherwise.
     * @param l Index of the one-past-the-end of the block.
     */
    DelayedSubsetBlock(std::shared_ptr<const typed_matrix<T, IDX> > p, size_t f, size_t l) : mat(p), first(f), last(l) {}

    /**
     * @copydoc DelayedSubsetBlock
     */
    DelayedSubsetBlock(std::shared_ptr<typed_matrix<T, IDX> > p, size_t f, size_t l) : mat(p), first(f), last(l) {}

    ~DelayedSubsetBlock() {}

public:
    const T* row(size_t r, T* buffer, size_t start, size_t end, workspace* work=nullptr) const {
        if constexpr(MARGIN == 0) {
            return mat->row(first + r, buffer, start, end, work);
        } else {
            return mat->row(r, buffer, first + start, first + end, work);
        }
    }

    const T* column(size_t c, T* buffer, size_t start, size_t end, workspace* work=nullptr) const {
        if constexpr(MARGIN == 0) {
            return mat->column(c, buffer, first + start, first + end, work);
        } else {
            return mat->column(first + c, buffer, start, end, work);
        }
    }

    using typed_matrix<T, IDX>::column;

    using typed_matrix<T, IDX>::row;

public:
    sparse_range<T, IDX> sparse_row(size_t r, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work=nullptr, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return mat->sparse_row(first + r, out_values, out_indices, start, end, work, sorted);
        } else {
            return subset_sparse<true>(r, out_values, out_indices, start, end, work, sorted);
        }
    }

    sparse_range<T, IDX> sparse_column(size_t c, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work=nullptr, bool sorted=true) const {
        if constexpr(MARGIN==0) {
            return subset_sparse<false>(c, out_values, out_indices, start, end, work, sorted);
        } else {
            return mat->sparse_column(c + first, out_values, out_indices, start, end, work, sorted);
        }
    }

    using typed_matrix<T, IDX>::sparse_column;

    using typed_matrix<T, IDX>::sparse_row;

public:
    /**
     * @return Number of rows after any subsetting is applied.
     */
    size_t nrow() const {
        if constexpr(MARGIN==0) {
            return last - first;
        } else {
            return mat->nrow();
        }
    }
    
    /**
     * @return Number of columns after any subsetting is applied.
     */
    size_t ncol() const {
        if constexpr(MARGIN==0) {
            return mat->ncol();
        } else {
            return last - first;
        }
    }

    /**
     * @return A null pointer or a shared pointer to a `workspace` object, depending on the underlying (pre-subsetted) matrix.
     */
    std::shared_ptr<workspace> new_workspace() const {
        return mat->new_workspace();
    }

    /**
     * @return The sparsity status of the underlying (pre-subsetted) matrix.
     */
    bool sparse() const {
        return mat->sparse();
    }

    /**
     * @return Whether the underlying (pre-subsetted) matrix prefers row access.
     */
    bool prefer_rows() const {
        return mat->prefer_rows();
    }

private:
    std::shared_ptr<const typed_matrix<T, IDX> > mat;
    size_t first, last;

    template<bool ROW>
    sparse_range<T, IDX> subset_sparse(size_t i, T* out_values, IDX* out_indices, size_t start, size_t end, workspace* work, bool sorted) const {
        sparse_range<T, IDX> output;

        if constexpr(ROW) {
            output = mat->sparse_row(i, out_values, out_indices, start + first, end + first, work, sorted);
        } else {
            output = mat->sparse_column(i, out_values, out_indices, start + first, end + first, work, sorted);
        }

        if (first) {
            if (out_indices != output.index) {
                std::copy(output.index, output.index + output.number, out_indices);
                output.index = out_indices;
            }
            for (size_t i = 0; i < output.number; ++i) {
                out_indices[i] -= first;
            }
        }

        return output;
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam MARGIN Dimension along which the addition is to occur.
 * If 0, the subset is applied to the rows; if 1, the subset is applied to the columns.
 * @tparam MAT A specialized `typed_matrix`, to be automatically deducted.
 * @tparam V Vector containing the subset indices, to be automatically deducted.
 *
 * @param p Pointer to the underlying (pre-subset) `typed_matrix`.
 * @param f Index of the start of the block. This should be a row index if `MARGIN = 0` and a column index otherwise.
 * @param l Index of the one-past-the-end of the block.
 *
 * @return A pointer to a `DelayedSubsetBlock` instance.
 */
template<int MARGIN, class MAT>
std::shared_ptr<MAT> make_DelayedSubsetBlock(std::shared_ptr<MAT> p, size_t f, size_t l) {
    return std::shared_ptr<MAT>(
        new DelayedSubsetBlock<MARGIN, typename MAT::value, typename MAT::index>(
            p,
            f,
            l
        )
    );
}

}

#endif
