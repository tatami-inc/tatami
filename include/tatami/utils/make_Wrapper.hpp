#ifndef TATAMI_BUFFERED_EXTRACTOR_HPP
#define TATAMI_BUFFERED_EXTRACTOR_HPP

#include <memory>
#include <vector>

#include "../base/Extractor.hpp"
#include "../base/Matrix.hpp"

/**
 * @file make_Wrapper.hpp
 * @brief User-friendly wrappers around `Extractor` instances.
 */

namespace tatami {

/********************
 *** Myopic dense ***
 ********************/

/**
 * @brief User-friendly wrapper around a `MyopicDenseExtractor`.
 *
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * Each instance of this class automatically allocates the buffer to be used in `MyopicDenseExtractor::fetch()`.
 */ 
template<typename Value_, typename Index_>
class MyopicDenseWrapper {
public:
    /**
     * @cond
     */
    MyopicDenseWrapper(std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > ext, Index_ length) :
        my_ext(std::move(ext)),
        my_buffer(tatami::cast_Index_to_container_size<decltype(my_buffer)>(length))
    {}
    /**
     * @endcond
     */

    /**
     * @return Length of the array extracted by `fetch()`.
     * This is constant for a given instance of this class. 
     */
    Index_ length() const {
        return my_buffer.size();
    }

    /**
     * @param i Index of the row/column of interest.
     * @return Pointer to an array of length `length()`, containing the values of the `i`-th row/column.
     * This is only guaranteed to be valid until the next call to this method.
     */
    const Value_* fetch(Index_ i) {
        return my_ext->fetch(i, my_buffer.data());
    }
private:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > my_ext;
    std::vector<Value_> my_buffer;
};

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract the full extent of each element of the target dimension.
 */
template<typename Value_, typename Index_>
MyopicDenseWrapper<Value_, Index_> make_DenseWrapper( 
    const Matrix<Value_, Index_>& matrix,
    bool row,
    const Options& opt
) {
    return MyopicDenseWrapper<Value_, Index_>(matrix.dense(row, opt), (row ? matrix.ncol() : matrix.nrow()));
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param block_start Start of the block on the non-target dimension.
 * @param block_length Length of the block on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract a contiguous block of each element of the target dimension.
 */
template<typename Value_, typename Index_>
MyopicDenseWrapper<Value_, Index_> make_DenseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    Index_ block_start,
    Index_ block_length,
    const Options& opt
) {
    return MyopicDenseWrapper<Value_, Index_>(matrix.dense(row, block_start, block_length, opt), block_length);
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param indices_ptr Pointer to an indexed subset on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract an indexed subset of each element of the target dimension.
 */
template<typename Value_, typename Index_>
MyopicDenseWrapper<Value_, Index_> make_DenseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    VectorPtr<Index_> indices_ptr,
    const Options& opt
) {
    const auto ilen = indices_ptr->size();
    return MyopicDenseWrapper<Value_, Index_>(matrix.dense(row, std::move(indices_ptr), opt), ilen);
}

/**********************
 *** Oracular dense ***
 **********************/

/**
 * @brief User-friendly wrapper around a `OracularDenseExtractor`.
 *
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * Each instance of this class automatically allocates the buffer to be used in `OracularDenseExtractor::fetch()`.
 */ 
template<typename Value_, typename Index_>
class OracularDenseWrapper {
public:
    /**
     * @cond
     */
    OracularDenseWrapper(std::unique_ptr<OracularDenseExtractor<Value_, Index_> > ext, Index_ length) :
        my_ext(std::move(ext)),
        my_buffer(tatami::cast_Index_to_container_size<decltype(my_buffer)>(length))
    {}
    /**
     * @endcond
     */

    /**
     * @return Length of the array extracted by `next()`.
     * This is constant for a given instance of this class. 
     */
    Index_ length() const {
        return my_buffer.size();
    }

    /**
     * @return Pointer to an array of length `length()`, containing the values of the next predicted row/column.
     * This is only guaranteed to be valid until the next call to this method.
     */
    const Value_* next() {
        return my_ext->fetch(my_buffer.data());
    }

private:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > my_ext;
    std::vector<Value_> my_buffer;
};

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param oracle Pointer to an `Oracle` that predicts the elements of the target dimensions to be extracted.
 * @param opt Further options.
 *
 * @return Wrapper to extract the full extent of each predicted element of the target dimension.
 */
template<typename Value_, typename Index_>
OracularDenseWrapper<Value_, Index_> make_DenseWrapper( 
    const Matrix<Value_, Index_>& matrix,
    bool row,
    std::shared_ptr<const Oracle<Index_> > oracle,
    const Options& opt
) {
    return OracularDenseWrapper<Value_, Index_>(matrix.dense(row, std::move(oracle), opt), (row ? matrix.ncol() : matrix.nrow()));
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param oracle Pointer to an `Oracle` that predicts the elements of the target dimensions to be extracted.
 * @param block_start Start of the block on the non-target dimension.
 * @param block_length Length of the block on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract a contiguous block of each predicted element of the target dimension.
 */
template<typename Value_, typename Index_>
OracularDenseWrapper<Value_, Index_> make_DenseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    std::shared_ptr<const Oracle<Index_> > oracle,
    Index_ block_start,
    Index_ block_length,
    const Options& opt
) {
    return OracularDenseWrapper<Value_, Index_>(matrix.dense(row, std::move(oracle), block_start, block_length, opt), block_length);
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param oracle Pointer to an `Oracle` that predicts the elements of the target dimensions to be extracted.
 * @param indices_ptr Pointer to an indexed subset on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract an indexed subset of each predicted element of the target dimension.
 */
template<typename Value_, typename Index_>
OracularDenseWrapper<Value_, Index_> make_DenseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    std::shared_ptr<const Oracle<Index_> > oracle,
    VectorPtr<Index_> indices_ptr,
    const Options& opt
) {
    const auto ilen = indices_ptr->size();
    return OracularDenseWrapper<Value_, Index_>(matrix.dense(row, std::move(oracle), std::move(indices_ptr), opt), ilen);
}

/*********************
 *** Myopic sparse ***
 *********************/

/**
 * @brief User-friendly wrapper around a `MyopicSparseExtractor`.
 *
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * Each instance of this class automatically allocates the buffer(s) to be used in `MyopicSparseExtractor::fetch()`.
 */ 
template<typename Value_, typename Index_>
class MyopicSparseWrapper {
public:
    /**
     * @cond
     */
    MyopicSparseWrapper(std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > ext, Index_ length, bool needs_value, bool needs_index) :
        my_ext(std::move(ext)),
        my_length(length)
    {
        if (needs_value) {
            tatami::resize_container_to_Index_size(my_vbuffer, length);
        }
        if (needs_index) {
            tatami::resize_container_to_Index_size(my_ibuffer, length);
        }
    }
    /**
     * @endcond
     */

    /**
     * @return Length of the arrays in the `SparseRange` returned by `fetch()`.
     * This is constant for a given instance of this class. 
     */
    Index_ length() const {
        return my_length;
    }

    /**
     * @param i Index of the row/column of interest.
     * @return `SparseRange` containing the number, values and/or indices of the structural non-zeros of the `i`-th row/column.
     * Pointers within this object are only guaranteed to be valid until the next call to this method.
     */
    SparseRange<Value_, Index_> fetch(Index_ i) {
        return my_ext->fetch(i, my_vbuffer.data(), my_ibuffer.data());
    }

private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > my_ext;
    Index_ my_length;
    std::vector<Value_> my_vbuffer;
    std::vector<Index_> my_ibuffer;
};

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract the full extent of each element of the target dimension.
 */
template<typename Value_, typename Index_>
MyopicSparseWrapper<Value_, Index_> make_SparseWrapper( 
    const Matrix<Value_, Index_>& matrix,
    bool row,
    const Options& opt
) {
    return MyopicSparseWrapper<Value_, Index_>(
        matrix.sparse(row, opt),
        (row ? matrix.ncol() : matrix.nrow()),
        opt.sparse_extract_value,
        opt.sparse_extract_index
    );
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param block_start Start of the block on the non-target dimension.
 * @param block_length Length of the block on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract a contiguous block of each element of the target dimension.
 */
template<typename Value_, typename Index_>
MyopicSparseWrapper<Value_, Index_> make_SparseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    Index_ block_start,
    Index_ block_length,
    const Options& opt
) {
    return MyopicSparseWrapper<Value_, Index_>(
        matrix.sparse(row, block_start, block_length, opt),
        block_length,
        opt.sparse_extract_value,
        opt.sparse_extract_index
    );
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param indices_ptr Pointer to an indexed subset on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract an indexed subset of each element of the target dimension.
 */
template<typename Value_, typename Index_>
MyopicSparseWrapper<Value_, Index_> make_SparseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    VectorPtr<Index_> indices_ptr,
    const Options& opt
) {
    const auto ilen = indices_ptr->size();
    return MyopicSparseWrapper<Value_, Index_>(
        matrix.sparse(row, std::move(indices_ptr), opt),
        ilen,
        opt.sparse_extract_value,
        opt.sparse_extract_index
    );
}

/***********************
 *** Oracular sparse ***
 ***********************/

/**
 * @brief User-friendly wrapper around a `OracularSparseExtractor`.
 *
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * Each instance of this class automatically allocates the buffer(s) to be used in `OracularSparseExtractor::fetch()`.
 */ 
template<typename Value_, typename Index_>
class OracularSparseWrapper {
public:
    /**
     * @cond
     */
    OracularSparseWrapper(std::unique_ptr<OracularSparseExtractor<Value_, Index_> > ext, Index_ length, bool needs_value, bool needs_index) :
        my_ext(std::move(ext)),
        my_length(length)
    {
        if (needs_value) {
            tatami::resize_container_to_Index_size(my_vbuffer, length);
        }
        if (needs_index) {
            tatami::resize_container_to_Index_size(my_ibuffer, length);
        }
    }
    /**
     * @endcond
     */

    /**
     * @return Length of the arrays in the `SparseRange` returned by `next()`.
     * This is constant for a given instance of this class. 
     */
    Index_ length() const {
        return my_length;
    }

    /**
     * @return `SparseRange` containing the number, values and/or indices of the structural non-zeros of the next predicted row/column.
     * Pointers within this object are only guaranteed to be valid until the next call to this method.
     */
    SparseRange<Value_, Index_> next() {
        return my_ext->fetch(my_vbuffer.data(), my_ibuffer.data());
    }

private:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > my_ext;
    Index_ my_length;
    std::vector<Value_> my_vbuffer;
    std::vector<Index_> my_ibuffer;
};

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param oracle Pointer to an `Oracle` that predicts the elements of the target dimensions to be extracted.
 * @param opt Further options.
 *
 * @return Wrapper to extract the full extent of each predicted element of the target dimension.
 */
template<typename Value_, typename Index_>
OracularSparseWrapper<Value_, Index_> make_SparseWrapper( 
    const Matrix<Value_, Index_>& matrix,
    bool row,
    std::shared_ptr<const Oracle<Index_> > oracle,
    const Options& opt
) {
    return OracularSparseWrapper<Value_, Index_>(
        matrix.sparse(row, std::move(oracle), opt),
        (row ? matrix.ncol() : matrix.nrow()),
        opt.sparse_extract_value,
        opt.sparse_extract_index
    );
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param oracle Pointer to an `Oracle` that predicts the elements of the target dimensions to be extracted.
 * @param block_start Start of the block on the non-target dimension.
 * @param block_length Length of the block on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract a contiguous block of each predicted element of the target dimension.
 */
template<typename Value_, typename Index_>
OracularSparseWrapper<Value_, Index_> make_SparseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    std::shared_ptr<const Oracle<Index_> > oracle,
    Index_ block_start,
    Index_ block_length,
    const Options& opt
) {
    return OracularSparseWrapper<Value_, Index_>(
        matrix.sparse(row, std::move(oracle), block_start, block_length, opt),
        block_length,
        opt.sparse_extract_value,
        opt.sparse_extract_index
    );
}

/**
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param matrix Matrix of interest.
 * @param row Whether the rows are the target dimension.
 * If `false`, the columns are the target dimension.
 * @param oracle Pointer to an `Oracle` that predicts the elements of the target dimensions to be extracted.
 * @param indices_ptr Pointer to an indexed subset on the non-target dimension.
 * @param opt Further options.
 *
 * @return Wrapper to extract an indexed subset of each predicted element of the target dimension.
 */
template<typename Value_, typename Index_>
OracularSparseWrapper<Value_, Index_> make_SparseWrapper(
    const Matrix<Value_, Index_>& matrix,
    bool row,
    std::shared_ptr<const Oracle<Index_> > oracle,
    VectorPtr<Index_> indices_ptr,
    const Options& opt
) {
    const auto ilen = indices_ptr->size();
    return OracularSparseWrapper<Value_, Index_>(
        matrix.sparse(row, std::move(oracle), std::move(indices_ptr), opt),
        ilen,
        opt.sparse_extract_value,
        opt.sparse_extract_index
    );
}

}

#endif
