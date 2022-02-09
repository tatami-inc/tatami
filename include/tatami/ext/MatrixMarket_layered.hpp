#ifndef TATAMI_MATRIX_MARKET_LAYERED_HPP
#define TATAMI_MATRIX_MARKET_LAYERED_HPP

#include <limits>
#include <cstdint>
#include <algorithm>
#include <vector>

#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

#include "MatrixMarket.hpp"
#include "LayeredMatrixData.hpp"
#include "layered_utils.hpp"

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
    bool sorted_by_column = true, sorted_by_row = true;
    size_t lastrow = 0, lastcol = 0;

    // For use with row-major extraction.
    std::vector<size_t> by_row_ptr8, by_row_ptr16, by_row_ptr32;
    std::vector<size_t> by_row_buffered_index;

    // For general use, computed on the run.
    std::vector<int> category;
    std::vector<size_t> lines_per_row;

    // For later use, computed in finish().
    std::vector<size_t> lines_per_category;
    std::vector<size_t> index;
    std::vector<size_t> rows_per_category;
    std::vector<size_t> permutation;

    size_t nrows, ncols;

private:
    void update_row_data() {
        if (by_row_buffered_index.empty()) {
            return;
        }

        int cat = category[lastrow];
        for (auto idx : by_row_buffered_index) {
            ++idx; // shift up 1 for the CSC pointers.
            switch (cat) {
                case 0:
                    ++(by_row_ptr8[idx]);
                    break;
                case 1:
                    ++(by_row_ptr16[idx]);
                    break;
                case 2:
                    ++(by_row_ptr32[idx]);
                    break;
                default:
                    break;
            }
        }
        by_row_buffered_index.clear();
    }

public:
    void setdim(size_t currow, size_t curcol, size_t) {
        nrows = currow;
        ncols = curcol;

        category.resize(nrows);
        lines_per_row.resize(nrows);

        by_row_ptr8.resize(ncols + 1);
        by_row_ptr16.resize(ncols + 1);
        by_row_ptr32.resize(ncols + 1);
        by_row_buffered_index.reserve(ncols);

        return;
    }

    void addline(size_t row, size_t col, size_t data, size_t) {
        // No need to check for negative values in 'data' when 'data' is an
        // unsigned type. To be honest, not much point checking for zeros
        // either, but we might as well try to prune them out.
        if (data == 0) {
            return;
        }

        auto chosen = layered_utils::categorize_row(data);
        if (chosen > category[row]) {
            category[row] = chosen;
        }

        if (sorted_by_row && (row < lastrow || (row == lastrow && col <= lastcol))) {
            sorted_by_row = false;
        } else {
            if (row != lastrow) {
                update_row_data();
            }
            by_row_buffered_index.push_back(col);
        }

        if (sorted_by_column && (col < lastcol || (col == lastcol && row <= lastrow))) {
            sorted_by_column = false;
        }

        lastrow = row;
        lastcol = col;

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

        if (sorted_by_row) {
            update_row_data();
            layered_utils::counts_to_offsets(by_row_ptr8);
            layered_utils::counts_to_offsets(by_row_ptr16);
            layered_utils::counts_to_offsets(by_row_ptr32);
        }
        
        // Computing the number of features and lines in each block.
        for (size_t i = 0; i < category.size(); ++i) {
            lines_per_category[category[i]] += lines_per_row[i];
        }

        // Computing the permutation.
        auto remap = layered_utils::compute_new_indices<size_t>(category);
        rows_per_category.swap(remap.per_category);
        index.swap(remap.new_indices);
        permutation.swap(remap.permutation);

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
struct LayeredBuilderByRow {
    LayeredBuilderByRow(const LineAssignments* ass) :
        assign(ass),

        off8(ass->by_row_ptr8),
        row8(off8.back()),
        dat8(off8.back()),

        off16(ass->by_row_ptr16),
        row16(off16.back()),
        dat16(off16.back()),

        off32(ass->by_row_ptr32),
        row32(off32.back()),
        dat32(off32.back())
    {}

public:
    const LineAssignments* assign;

    std::vector<size_t>   off8;
    std::vector<ROW>      row8;
    std::vector<uint8_t>  dat8;

    std::vector<size_t>   off16;
    std::vector<ROW>      row16;
    std::vector<uint16_t> dat16;

    std::vector<size_t>   off32;
    std::vector<ROW>      row32;
    std::vector<uint32_t> dat32;

public:
    void setdim(size_t, size_t, size_t) {}

    void addline(size_t row, size_t col, size_t data, size_t) {
        if (data == 0) { // match the LayeredAssignments behavior.
            return;
        }

        auto idx = assign->index[row];
        switch (assign->category[row]) {
            case 0:
                {
                    auto& counter8 = off8[col];
                    row8[counter8] = idx;
                    dat8[counter8] = data;
                    ++counter8;
                }
                break;
            case 1:
                {
                    auto& counter16 = off16[col];
                    row16[counter16] = idx;
                    dat16[counter16] = data;
                    ++counter16;
                }
                break;
            case 2:
                {
                    auto& counter32 = off32[col];
                    row32[counter32] = idx;
                    dat32[counter32] = data;
                    ++counter32;
                }
                break;
            default:
                break;
        }
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
        
        return layered_utils::consolidate_submatrices<T, IDX>(
            assign->rows_per_category,
            std::move(row8),
            std::move(dat8),
            std::move(assign->by_row_ptr8),
            std::move(row16),
            std::move(dat16),
            std::move(assign->by_row_ptr16),
            std::move(row32),
            std::move(dat32),
            std::move(assign->by_row_ptr32));
    }
};
/**
 * @endcond
 */

/**
 * @cond
 */
template<typename ROW>
struct LayeredBuilderByColumn {
    LayeredBuilderByColumn(const LineAssignments* ass) : 
        assign(ass),
        ptr8(ass->ncols + 1),
        ptr16(ass->ncols + 1),
        ptr32(ass->ncols + 1)
    {
        row8.reserve(assign->lines_per_category[0]);
        dat8.reserve(assign->lines_per_category[0]);

        row16.reserve(assign->lines_per_category[1]);
        dat16.reserve(assign->lines_per_category[1]);

        row32.reserve(assign->lines_per_category[2]);
        dat32.reserve(assign->lines_per_category[2]);

        return;
    }

public:
    const LineAssignments* assign;

    std::vector<ROW>      row8;
    std::vector<uint8_t>  dat8;
    std::vector<size_t>   ptr8;

    std::vector<ROW>      row16;
    std::vector<uint16_t> dat16;
    std::vector<size_t>   ptr16;

    std::vector<ROW>      row32;
    std::vector<uint32_t> dat32;
    std::vector<size_t>   ptr32;

public:
    void setdim(size_t, size_t, size_t) {}

    void addline(size_t row, size_t col, size_t data, size_t) {
        if (data == 0) { // match the LayeredAssignments behavior.
            return;
        }

        auto idx = assign->index[row];
        ++col;

        switch (assign->category[row]) {
            case 0:
                row8.push_back(idx);
                dat8.push_back(data);
                ++ptr8[col];
                break;
            case 1:
                row16.push_back(idx);
                dat16.push_back(data);
                ++ptr16[col];
                break;
            case 2:
                row32.push_back(idx);
                dat32.push_back(data);
                ++ptr32[col];
                break;
            default:
                break;
        }
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
        
        layered_utils::counts_to_offsets(ptr8);
        layered_utils::counts_to_offsets(ptr16);
        layered_utils::counts_to_offsets(ptr32);

        return layered_utils::consolidate_submatrices<T, IDX>(
            assign->rows_per_category,
            std::move(row8),
            std::move(dat8),
            std::move(ptr8),
            std::move(row16),
            std::move(dat16),
            std::move(ptr16),
            std::move(row32),
            std::move(dat32),
            std::move(ptr32));
    }
};
/**
 * @endcond
 */

/**
 * @cond
 */
template<typename ROW>
struct LayeredBuilderUnsorted {
    LayeredBuilderUnsorted(const LineAssignments* ass) : assign(ass) {
        row8.reserve(assign->lines_per_category[0]);
        col8.reserve(assign->lines_per_category[0]);
        dat8.reserve(assign->lines_per_category[0]);

        row16.reserve(assign->lines_per_category[1]);
        col16.reserve(assign->lines_per_category[1]);
        dat16.reserve(assign->lines_per_category[1]);

        row32.reserve(assign->lines_per_category[2]);
        col32.reserve(assign->lines_per_category[2]);
        dat32.reserve(assign->lines_per_category[2]);
        return;
    }

public:
    const LineAssignments* assign;

    std::vector<ROW>      row8;
    std::vector<uint32_t> col8;
    std::vector<uint8_t>  dat8;

    std::vector<ROW>      row16;
    std::vector<uint32_t> col16;
    std::vector<uint16_t> dat16;

    std::vector<ROW>      row32;
    std::vector<uint32_t> col32;
    std::vector<uint32_t> dat32;

public:
    void setdim(size_t, size_t, size_t) {}

    void addline(size_t row, size_t col, size_t data, size_t line) {
        if (data == 0) { // match the LayeredAssignments behavior.
            return;
        }

        auto idx = assign->index[row];
        switch (assign->category[row]) {
        case 0:
            row8.push_back(idx);
            col8.push_back(col);
            dat8.push_back(data);
            break;
        case 1:
            row16.push_back(idx);
            col16.push_back(col);
            dat16.push_back(data);
            break;
        case 2:
            row32.push_back(idx);
            col32.push_back(col);
            dat32.push_back(data);
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

        return layered_utils::consolidate_submatrices(std::move(collated), ass.ncols);
    }
};
/**
 * @endcond
 */

/**
 * @cond
 */
template<typename T = double, typename IDX = int, typename RowIndex, class Function>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_internal(LineAssignments& ass, Function& process) {
    LayeredMatrixData<T, IDX> output;
    output.permutation = ass.permutation;

    if (ass.sorted_by_column) {
        LayeredBuilderByColumn<RowIndex> builder(&ass);
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    } else if (ass.sorted_by_row) {
        LayeredBuilderByRow<RowIndex> builder(&ass);
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    } else {
        LayeredBuilderUnsorted<RowIndex> builder(&ass);
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    }

    return output;
}

template<typename T = double, typename IDX = int, class Function>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_internal(Function process) {
    LineAssignments ass;
    process(ass);
    ass.finish();

    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (ass.nrows <= max16) {
        return load_layered_sparse_matrix_internal<T, IDX, uint16_t>(ass, process);
    } else {
        return load_layered_sparse_matrix_internal<T, IDX, IDX>(ass, process);
    }
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
 * This function loads a layered sparse integer matrix from a Matrix Market coordinate file.
 * The aim is to reduce memory usage by storing each gene's counts in the smallest unsigned integer type that can hold them.
 * See the documentation for `LayeredMatrixData` for more details.
 *
 * On a related note, this function will automatically use `uint16_t` values to store the internal row indices if the number of rows is less than 65536.
 * This aims to further reduce memory usage in most cases, e.g., gene count matrices usually have fewer than 50000 rows.
 * Note that the internal storage is orthogonal to the choice of `IDX` in the `tatami::Matrix` interface.
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
