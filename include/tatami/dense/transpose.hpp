#ifndef TATAMI_TRANSPOSE_HPP
#define TATAMI_TRANSPOSE_HPP

#include <algorithm>

/**
 * @file transpose.hpp
 * @brief Transpose a dense array.
 */

namespace tatami {

/**
 * @tparam Input_ Input type.
 * @tparam Output_ Output type.
 * @param[in] input Pointer to an array containing a row-major matrix with `nrow` rows and `ncol` columns.
 * Elements within each row should be contiguous but consecutive rows can be separated by a constant stride, see `input_stride`.
 * The array should have at least `(nrow - 1) * input_stride + ncol` addressable elements.
 * @param nrow Number of rows in the matrix stored at `input`.
 * @param ncol Number of columns in the matrix stored at `input`.
 * @param input_stride Distance between corresponding entries on consecutive rows of the `input` matrix.
 * This should be greater than or equal to `ncol`.
 * @param[out] output Pointer to an array in which to store the transpose of the matrix in `input`.
 * On output, this stores a row-major matrix with `ncol` rows and `nrow` columns.
 * Elements within each row should be contiguous but consecutive rows can be separated by a constant stride, see `output_stride`.
 * The array should have at least `(ncol - 1) * output_stride + nrow` addressable elements.
 * It is assumed that `output` does not alias `input`.
 * @param output_stride Distance between corresponding entries on consecutive rows of the `output` matrix.
 * This should be greater than or equal to `nrow`.
 *
 * This function is intended for developers of `Matrix` subclasses who need to do some transposition, e.g., for dense chunks during caching.
 * The `*_stride` arguments allow `input` and `output` to refer to submatrices of larger arrays.
 *
 * The argument descriptions refer to row-major matrices only for the sake of convenience.
 * This function is equally applicable to column-major matrices, just replace all instances of "row" with "column" and vice versa. 
 */
template<typename Input_, typename Output_>
void transpose(const Input_* input, std::size_t nrow, std::size_t ncol, std::size_t input_stride, Output_* output, std::size_t output_stride) {
    if ((nrow == 1 && output_stride == 1) || (ncol == 1 && input_stride == 1)) {
        std::copy_n(input, nrow * ncol, output);
        return;
    }

    // Using a blockwise strategy to perform the transposition,
    // in order to be more input-friendly.
    constexpr std::size_t block = 16;
    std::size_t col_start = 0;
    while (col_start < ncol) {
        std::size_t col_end = col_start + std::min(block, ncol - col_start);

        std::size_t row_start = 0;
        while (row_start < nrow) {
            std::size_t row_end = row_start + std::min(block, nrow - row_start);
            for (std::size_t c = col_start; c < col_end; ++c) {
                for (std::size_t r = row_start; r < row_end; ++r) {
                    output[c * output_stride + r] = input[r * input_stride + c];
                }
            }

            row_start = row_end;
        }
        col_start = col_end;
    }
}

/**
 * @tparam Input_ Input type.
 * @tparam Output_ Output type.
 * @param[in] input Pointer to an array containing a row-major matrix with `nrow` rows and `ncol` columns.
 * The array should have at least `nrow * ncol` addressable elements, and all elements should be stored contiguously in the array.
 * @param nrow Number of rows in the matrix stored at `input`.
 * @param ncol Number of columns in the matrix stored at `input`.
 * @param[out] output Pointer to an array of length `nrow * ncol`.
 * On output, this will hold the transpose of the matrix represented by `input`,
 * i.e., a row-major matrix with `ncol` rows and `nrow` columns.
 * It is assumed that `output` does not alias `input`.
 *
 * This function is intended for developers of `Matrix` subclasses who need to do some transposition, e.g., for dense chunks during caching.
 * Users should instead construct a `DelayedTranspose` object to perform a memory-efficient delayed transposition,
 * or use `convert_to_dense()` to convert their dense data into the desired storage layout.
 *
 * The argument descriptions refer to row-major matrices only for the sake of convenience.
 * This function is equally applicable to column-major matrices, just replace all instances of "row" with "column" and vice versa. 
 */
template<typename Input_, typename Output_>
void transpose(const Input_* input, std::size_t nrow, std::size_t ncol, Output_* output) {
    transpose(input, nrow, ncol, ncol, output, nrow);
}

// COMMENT:
// I tried really hard to make an in-place version, but it's too frigging complicated for non-square matrices.
// It can be done, but I can't see a way to do it efficiently as you end up hopping all over the matrix (a la in-place reordering).
// There doesn't seem to be any opportunity to do it in blocks for cache-friendliness;
// the displaced values from the original matrix don't form a corrresponding block in the transposed matrix.
// Perhaps this is a skill issue but at least Eigen agrees with my assessment, as they just make a new copy when the matrix is not square.
// (See https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/Core/Transpose.h#L293 for the relevant code.)
// I'm not going to optimize the square case because tatami rarely, if ever, deals in situations with square matrices.

}

#endif
