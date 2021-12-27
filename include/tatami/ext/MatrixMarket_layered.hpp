#ifndef TATAMI_MATRIX_MARKET_LAYERED_HPP
#define TATAMI_MATRIX_MARKET_LAYERED_HPP

#include <limits>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cctype>
#include <string>

#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

#include "buffin/parse_text_file.hpp"

#ifdef TATAMI_USE_ZLIB
#include "buffin/parse_gzip_file.hpp"
#include "buffin/parse_zlib_buffer.hpp"
#endif

/**
 * @file MatrixMarket_layered.hpp
 *
 * @brief Create a layered sparse matrix from the Matrix Market coordinate format.
 */

namespace tatami {

namespace MatrixMarket {

/**
 * @cond
 */
struct LineAssignments {
    LineAssignments() : rows_per_category(3), lines_per_category(3) {}

public:
    std::vector<uint8_t> category;
    std::vector<size_t> index;
    std::vector<size_t> rows_per_category;
    std::vector<size_t> lines_per_category;
    std::vector<size_t> lines_per_row;
    std::vector<size_t> permutation;

    size_t nrows, ncols;

public:
    void setdim(size_t currow, size_t curcol, size_t) {
        nrows = currow;
        ncols = curcol;

        category.resize(nrows);
        index.resize(nrows);
        permutation.resize(nrows);
        lines_per_row.resize(nrows);
        return;
    }

    void addline(size_t row, size_t col, size_t data, size_t) {
        // Assigning each line to a block based on its integer size.
        constexpr int max8 = std::numeric_limits<uint8_t>::max();
        constexpr int max16 = std::numeric_limits<uint16_t>::max();

        if (data > max16) {
            category[row] = std::max(category[row], static_cast<uint8_t>(2));
        } else if (data > max8) {
            category[row] = std::max(category[row], static_cast<uint8_t>(1));
        }

        ++lines_per_row[row];
        return;
    }

public:
    BaseMMParser base;

    template<typename B>
    void add(const B* buffer, size_t n) {
        base.add(buffer, n, *this);
        return;
    }

    void finish() {
        base.finish(*this);
        
        // Computing the number of features and lines in each block.
        auto iIt = index.begin();
        auto cIt = lines_per_row.begin();
        for (auto f : category) {
            auto& current = rows_per_category[f];

            (*iIt) = current;
            ++current;
            ++iIt;

            lines_per_category[f] += *cIt;
            ++cIt;
        }

        // Computing the permutation.
        auto cumsum = rows_per_category;
        size_t last = 0;
        for (auto& x : cumsum) {
            std::swap(x, last);
            last += x;
        }
        for (size_t i = 0; i < nrows; ++i) {
            permutation[i] = cumsum[category[i]] + index[i];
        }

        return;
    }
};
/**
 * @endcond
 */

/**
 * @cond
 */
template<typename ROW>
struct LayeredBuilder {
    LayeredBuilder(const LineAssignments* ass) : 
        assign(ass),
        row8(assign->lines_per_category[0]),
        col8(assign->lines_per_category[0]),
        dat8(assign->lines_per_category[0]),
        row16(assign->lines_per_category[1]),
        col16(assign->lines_per_category[1]),
        dat16(assign->lines_per_category[1]),
        row32(assign->lines_per_category[2]),
        col32(assign->lines_per_category[2]),
        dat32(assign->lines_per_category[2]) {}

public:
    const LineAssignments* assign;

    std::vector<ROW>      row8;
    std::vector<uint32_t> col8;
    std::vector<uint8_t>  dat8;
    size_t counter8 = 0;

    std::vector<ROW>      row16;
    std::vector<uint32_t> col16;
    std::vector<uint16_t> dat16;
    size_t counter16 = 0;

    std::vector<ROW>      row32;
    std::vector<uint32_t> col32;
    std::vector<uint32_t> dat32;
    size_t counter32 = 0;

public:
    void setdim(size_t, size_t, size_t) {}

    void addline(size_t row, size_t col, size_t data, size_t line) {
        auto idx = assign->index[row];
        switch (assign->category[row]) {
        case 0:
            row8[counter8] = idx;
            col8[counter8] = col;
            dat8[counter8] = data;
            ++counter8;
            break;
        case 1:
            row16[counter16] = idx;
            col16[counter16] = col;
            dat16[counter16] = data;
            ++counter16;
            break;
        case 2:
            row32[counter32] = idx;
            col32[counter32] = col;
            dat32[counter32] = data;
            ++counter32;
            break;
        }
    }

private:
    template<typename T, typename IDX, typename U, typename V, typename W>
    static std::shared_ptr<Matrix<T, IDX> > create_sparse_matrix(size_t nr, size_t nc, U& values, V& rows, W& cols) {
        auto indptrs = compress_sparse_triplets<false>(nr, nc, values, rows, cols);
        return std::shared_ptr<Matrix<T, IDX> >(new CompressedSparseColumnMatrix<T, IDX, U, V, decltype(indptrs)>(nr, nc, std::move(values), std::move(rows), std::move(indptrs)));
    }

public:
    BaseMMParser base;

    template<typename B>
    void add(const B* buffer, size_t n) {
        base.add(buffer, n, *this);
        return;
    }

    template<typename T, typename IDX>
    std::shared_ptr<Matrix<T, IDX> > finish() {
        base.finish(*this);

        std::vector<std::shared_ptr<Matrix<T, IDX> > > collated;
        const auto& ass = *assign;

        if (ass.rows_per_category[0]) {
            collated.push_back(create_sparse_matrix<T, IDX>(ass.rows_per_category[0], ass.ncols, dat8, row8, col8));
        }
        if (ass.rows_per_category[1]) {
            collated.push_back(create_sparse_matrix<T, IDX>(ass.rows_per_category[1], ass.ncols, dat16, row16, col16));
        }
        if (ass.rows_per_category[2]) {
            collated.push_back(create_sparse_matrix<T, IDX>(ass.rows_per_category[2], ass.ncols, dat32, row32, col32));
        }

        if (collated.size() == 0) {
            return create_sparse_matrix<T, IDX>(0, ass.ncols, dat8, row8, col8);
        } else if (collated.size() == 1) { 
            return collated[0];
        } else {
            return make_DelayedBind<0>(std::move(collated));
        }
    }
};
/**
 * @endcond
 */

/**
 * @brief Pointer and permutations for a layered sparse matrix.
 *
 * This holds a pointer to a "layered sparse matrix", see `load_layered_sparse_matrix()` for details.
 * It also holds the permutation vector that was generated as part of the layering process.
 *
 * @tparam T Type of the matrix values.
 * @tparam IDX Integer type for the index.
 */
template<typename T = double, typename IDX = int>
struct LayeredMatrixData {
    /**
     * Pointer to the matrix data.
     * Note that rows will not follow their original order, see `permutation`.
     */
    std::shared_ptr<Matrix<T, IDX> > matrix;

    /**
     * Permutation vector indicating the position of each original row in `matrix`.
     * Specifically, for an original row `r`, its new row index in `matrix` is defined as `permutation[i]`.
     */
    std::vector<size_t> permutation;
};

/**
 * @cond
 */
template<typename T = double, typename IDX = int, class Function>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_internal(Function process) {
    LineAssignments ass;
    process(ass);
    ass.finish();

    LayeredMatrixData<T, IDX> output;
    output.permutation = ass.permutation;

    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (ass.nrows > max16) {
        LayeredBuilder<uint16_t> builder(&ass);
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    } else {
        LayeredBuilder<IDX> builder(&ass);
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    }

    return output;
}
/**
 * @endcond
 */

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param bufsize Size of the buffer (in bytes) to use when reading from file.
 * 
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This loads a sparse integer matrix from a Matrix Market coordinate file, with a twist:
 * values are stored in the smallest integer type that will fit all entries on the same row.
 * This is achieved by reordering the rows into separate `tatami::CompressedSparseMatrix` submatrices,
 * where each submatrix contains all rows that can fit into a certain integer type.
 * Submatrices are then combined together using an instance of the `tatami::DelayedBind` class.
 * The idea is to reduce memory usage in common situations involving small integers.
 *
 * Note that the rows of the output matrix will be in a different order from the data in the file.
 * This is a consequence of the reordering to ensure that rows of the same type can be put into the same submatrix.
 * The returned `LayeredMatriData` object will contain a permutation vector indicating the new position of each original row.
 * This can be passed to `tatami::DelayedSubset` to restore the original order, if so desired.
 *
 * Currently, only unsigned integers are supported.
 */
template<typename T = double, typename IDX = int>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix(const char * filepath, size_t bufsize = 65536) {
    return load_layered_sparse_matrix_internal([&](auto& obj) -> void {
        buffin::parse_text_file(filepath, obj, bufsize);        
    });
}

/**
 * @param buffer Array containing the contents of a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param n Length of the array.
 * 
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 * @tparam B Type of the buffer, usually `char` or `unsigned char`.
 *
 * This is equivalent to `load_layered_sparse_matrix()` but assumes that the entire file has been read into `buffer`.
 */
template<typename T = double, typename IDX = int, typename B>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_from_buffer(const B* buffer, size_t n) {
    return load_layered_sparse_matrix_internal([&](auto& obj) -> void {
        obj.add(buffer, n);
    });
}

#ifdef TATAMI_USE_ZLIB

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored with Gzip compression.
 * @param bufsize Size of the buffer to use for decompression, in bytes.
 *
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This is a version of `load_layered_sparse_matrix()` for loading in Gzip-compressed Matrix Market files.
 * To make this function available, make sure to define `TATAMI_USE_ZLIB` and compile with **zlib** support.
 */
template<typename T = double, typename IDX = int>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_gzip(const char * filepath, int bufsize = 65536) {
    return load_layered_sparse_matrix_internal([&](auto& obj) -> void {
        buffin::parse_gzip_file(filepath, obj, bufsize);
    });
}

/**
 * @param buffer Array containing the contents of a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param n Length of the array.
 * @param bufsize Size of the buffer to use for decompression, in bytes.
 * 
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This is equivalent to `load_layered_sparse_matrix_gzip()` but assumes that the entire file has been read into `buffer`.
 */
template<typename T = double, typename IDX = int>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_from_buffer_gzip(const unsigned char* buffer, size_t n, size_t bufsize = 65536) {
    auto rebuffer = const_cast<unsigned char*>(buffer);
    return load_layered_sparse_matrix_internal([&](auto& obj) -> void {
        buffin::parse_zlib_buffer(rebuffer, n, obj, 3, bufsize);
    });
}

#endif

}

}

#endif
