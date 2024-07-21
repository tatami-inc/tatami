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
 * @param nrow Number of rows in the matrix stored at `input`.
 * @param ncol Number of columns in the matrix stored at `input`.
 * @param input_stride Distance between corresponding entries on consecutive rows of the `input` matrix.
 * @param[out] output Pointer to an array of length `nrow * ncol`.
 * On output, this will hold the transpose of the matrix represented by `input`.
 * @param output_stride Distance between corresponding entries on consecutive columns of the `output` matrix.
 *
 * This function is intended for developers of `Matrix` subclasses who need to do some transposition, e.g., for dense chunks during caching.
 * The `*_stride` arguments allow `input` and `output` to refer to submatrices of larger arrays.
 */
template<typename Input_, typename Output_>
void transpose(const Input_* input, size_t nrow, size_t ncol, size_t input_stride, Output_* output, size_t output_stride) {
    if ((nrow == 1 && output_stride == 1) || (ncol == 1 && input_stride == 1)) {
        std::copy_n(input, nrow * ncol, output);
        return;
    }

    // Using a blockwise strategy to perform the transposition,
    // in order to be more input-friendly.
    constexpr size_t block = 16;
    size_t col_start = 0;
    while (col_start < ncol) {
        size_t col_end = col_start + std::min(block, ncol - col_start);

        size_t row_start = 0;
        while (row_start < nrow) {
            size_t row_end = row_start + std::min(block, nrow - row_start);

            // We use offsets instead of directly performing pointer
            // arithmetic, to avoid creating an invalid pointer address (which
            // is UB, even if unused) after the last inner loop iteration. 
            size_t input_offset = col_start + row_start * input_stride;
            size_t output_offset = col_start * output_stride + row_start;

            for (size_t c = col_start; c < col_end; ++c, ++input_offset, output_offset += output_stride) {
                auto input_offset_copy = input_offset;
                auto output_offset_copy = output_offset;

                for (size_t r = row_start; r < row_end; ++r, input_offset_copy += input_stride, ++output_offset_copy) {
                    output[output_offset_copy] = input[input_offset_copy];
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
 * @param nrow Number of rows in the matrix stored at `input`.
 * @param ncol Number of columns in the matrix stored at `input`.
 * @param[out] output Pointer to an array of length `nrow * ncol`.
 * On output, this will hold the transpose of the matrix represented by `input`.
 *
 * This function is intended for developers of `Matrix` subclasses who need to do some transposition, e.g., for dense chunks during caching.
 * Users should instead construct a `DelayedTranspose` object to perform a memory-efficient delayed transposition,
 * or use `convert_to_dense()` to convert their dense data into the desired storage layout.
 */
template<typename Input_, typename Output_>
void transpose(const Input_* input, size_t nrow, size_t ncol, Output_* output) {
    transpose(input, nrow, ncol, ncol, output, nrow);
}

}

#endif
