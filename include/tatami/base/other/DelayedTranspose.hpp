#ifndef TATAMI_DELAYED_TRANSPOSE
#define TATAMI_DELAYED_TRANSPOSE

#include "../Matrix.hpp"
#include <memory>

/**
 * @file DelayedTranspose.hpp
 *
 * @brief Delayed transposition.
 *
 * This is equivalent to the `DelayedAperm` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @brief Delayed transposition of a matrix.
 *
 * Implements delayed transposition of a matrix.
 * This operation is "delayed" in that it is only evaluated on request, e.g., with `row()` or friends.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 */
template<typename Value_, typename Index_>
class DelayedTranspose : public Matrix<Value_, Index_> {
public:
    /**
     * @param p Pointer to the underlying (pre-transpose) matrix.
     */
    DelayedTranspose(std::shared_ptr<const Matrix<Value_, Index_> > p) : mat(std::move(p)) {}

private:
    std::shared_ptr<const Matrix<Value_, Index_> > mat;

public:
    Index_ nrow() const {
        return mat->ncol();
    }

    Index_ ncol() const {
        return mat->nrow();
    }

    bool sparse() const {
        return mat->sparse();
    }

    bool prefer_rows() const {
        return !mat->prefer_rows();
    }

    std::pair<double, double> dimension_preference () const {
        auto temp = mat->dimension_preference();
        return std::pair<double, double>(temp.second, temp.first);
    }

    using Matrix<Value_, Index_>::dense_row;

    using Matrix<Value_, Index_>::dense_column;

    using Matrix<Value_, Index_>::sparse_row;

    using Matrix<Value_, Index_>::sparse_column;

public:
    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_row(const Options<Index_>& opt) const {
        return mat->dense_column(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return mat->dense_column(block_start, block_length, opt);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const Options<Index_>& opt) const {
        return mat->dense_column(std::move(indices), opt);
    }

    std::unique_ptr<FullDenseExtractor<Value_, Index_> > dense_column(const Options<Index_>& opt) const {
        return mat->dense_row(opt);
    }

    std::unique_ptr<BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return mat->dense_row(block_start, block_length, opt);
    }

    std::unique_ptr<IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const Options<Index_>& opt) const {
        return mat->dense_row(std::move(indices), opt);
    }

public:
    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_row(const Options<Index_>& opt) const {
        return mat->sparse_column(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return mat->sparse_column(block_start, block_length, opt);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const Options<Index_>& opt) const {
        return mat->sparse_column(std::move(indices), opt);
    }

    std::unique_ptr<FullSparseExtractor<Value_, Index_> > sparse_column(const Options<Index_>& opt) const {
        return mat->sparse_row(opt);
    }

    std::unique_ptr<BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const Options<Index_>& opt) const {
        return mat->sparse_row(block_start, block_length, opt);
    }

    std::unique_ptr<IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const Options<Index_>& opt) const {
        return mat->sparse_row(std::move(indices), opt);
    }
};

/**
 * A `make_*` helper function to enable partial template deduction of supplied types.
 *
 * @tparam Matrix_ A realized `Matrix` class, possibly one that is `const`.
 *
 * @param p Pointer to a `Matrix` instance.
 *
 * @return A pointer to a `DelayedTranspose` instance.
 */
template<class Matrix_>
std::shared_ptr<Matrix_> make_DelayedTranspose(std::shared_ptr<Matrix_> p) {
    typedef typename Matrix_::value_type Value_;
    typedef typename Matrix_::index_type Index_;
    return std::shared_ptr<Matrix_>(new DelayedTranspose<Value_, Index_>(std::move(p)));
}

}

#endif
