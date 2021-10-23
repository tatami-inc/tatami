#ifndef SCRAN_LOAD_MATRIX_MARKET_HPP
#define SCRAN_LOAD_MATRIX_MARKET_HPP

#include <limits>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cctype>

#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

#ifdef TATAMI_USE_ZLIB
#include "zlib.h"
#include <array>
#endif

/**
 * @file MatrixMarket.hpp
 *
 * @brief Read a sparse matrix from the Matrix Market coordinate format.
 */

namespace tatami {

namespace MatrixMarket {

/**
 * @cond
 */
template<typename T>
size_t process_triplet_line(const char* buffer, T& arg1, T& arg2, T& arg3, size_t n) {
    auto read = [=](size_t& i, T& arg) -> bool {
        // These had better be positive integers, 
        // otherwise this will throw a bunch of errors.
        while (i < n && std::isspace(buffer[i])) { 
            ++i;
        }
        if (i == n) {
            return false;
        }

        auto prev = i;
        arg = 0;
        while (i < n && std::isdigit(buffer[i])) {
            arg *= 10;
            arg += (buffer[i] - '0');
            ++i;
        }
        if (i == n || i == prev) {
            return false;
        }

        if (i < n && !std::isspace(buffer[i]) && buffer[i] != '\0') {
            throw std::runtime_error("values should be non-negative integers");
        }
        return true;
    };

    size_t i = 0;
    if (!read(i, arg1)) {
        return 0;
    }
    if (!read(i, arg2)) {
        return 0;
    }
    if (!read(i, arg3)) {
        return 0;
    }

    return i + 1;
}

inline size_t read_to_eol(const char * buffer, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (buffer[i] == '\n' || buffer[i] == '\0') {
            return i + 1;
        }
    }
    return 0; // out of range
}

struct TextReader {
    TextReader(const char* path) : filepath(path) {}

    template<class Object>
    void operator()(Object& obj) {
        std::ifstream in(filepath);
        if (!in) {
            throw std::runtime_error("failed to open file");
        }
        std::string line;
        int counter = 0;
        while (std::getline(in, line)) {
            auto status = obj.add(line.c_str(), line.size() + 1); // for the null terminator
            if (!status) {
                throw std::runtime_error(std::string("failed to process line ") + std::to_string(counter + 1));
            }
            ++counter;
        }
        return;
    }
private:
    const char* filepath;
};

#ifdef TATAMI_USE_ZLIB

struct GzipReader {
    GzipReader (const char* path, const int size = 16348) : filepath(path), bufsize(size), buffer(bufsize) {}
    const int bufsize;
    const char* filepath;
    std::vector<unsigned char> buffer;

    struct GZFile {
        GZFile(const char* path) : handle(gzopen(path, "rb")) {
            if (!handle) {
                throw std::runtime_error("failed to open file");
            }
        }

        ~GZFile() {
            gzclose(handle);
        }

        // Delete the remaining constructors.
        GZFile(const GZFile&) = delete;
        GZFile(GZFile&&) = delete;
        GZFile& operator=(const GZFile&) = delete;
        GZFile& operator=(GZFile&&) = delete;

        gzFile handle;
    };

    template<class Object>
    void operator()(Object& obj) {
        GZFile gz(filepath);

        size_t read = 1; // any positive value, to get the ball rolling.
        size_t leftovers = 0;

        while (read) {
            read = gzread(gz.handle, buffer.data() + leftovers, bufsize - leftovers);
            size_t current_stored = read + leftovers;

            if (read == 0) {
                // Making sure we have a terminating newline on the last line.
                if (gzeof(gz.handle)) { 
                    if (current_stored && buffer[current_stored-1]!='\n') {
                        buffer[current_stored] = '\n';
                        ++current_stored;
                    }
                } else {
                    int dummy;
                    throw std::runtime_error(gzerror(gz.handle, &dummy));
                }
            }

            // Adding whole lines.
            size_t last_processed = 0, total_processed = 0;
            do {
                last_processed = obj.add((char*)buffer.data() + total_processed, current_stored - total_processed);
                total_processed += last_processed;
            } while (last_processed);

            // Rotating what's left to the front for the next cycle.
            leftovers = current_stored - total_processed;
            for (size_t i = 0; i < leftovers; ++i) {
                buffer[i] = buffer[total_processed + i];
            }
        }

        return;
    }
};

#endif

/**
 * @endcond
 */

/*******************************************************************************
 **************************** Simple loading functions *************************
 *******************************************************************************/

/**
 * @cond
 */
template<typename T, typename IDX> 
struct SimpleBuilder {
    std::vector<uint16_t> short_rows;
    std::vector<IDX> long_rows;
    bool use_short_rows = false;

    std::vector<uint16_t> short_cols;
    std::vector<IDX> long_cols;
    bool use_short_cols = false;

    std::vector<T> values;
    int current_line = 0;

    size_t nrows, ncols, nlines;
    bool passed_preamble = false;

    size_t add(const char* buffer, size_t n) {
        if (buffer[0] == '%') {
            // TODO: should probably check for 'coordinate integer'.
            return read_to_eol(buffer, n);

        } else if (!passed_preamble) {
            passed_preamble = true;
            auto read = process_triplet_line(buffer, nrows, ncols, nlines, n);
            if (read == 0) {
                return 0;
            }

            constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
            if (nrows <= max16) {
                short_rows.resize(nlines);
                use_short_rows = true;
            } else {
                long_rows.resize(nlines);
            }

            if (ncols <= max16) {
                short_cols.resize(nlines);
                use_short_cols = true;
            } else {
                long_cols.resize(nlines);
            }

            values.resize(nlines);
            return read;
        }

        int row, col, data;
        auto read = process_triplet_line(buffer, row, col, data, n);
        if (read == 0){
            return 0;
        }

        --row;
        if (use_short_rows) {
            short_rows[current_line] = row;
        } else {
            long_rows[current_line] = row;
        }

        --col;
        if (use_short_cols) {
            short_cols[current_line] = col;
        } else {
            long_cols[current_line] = col;
        }

        values[current_line] = data;
        ++current_line;

        return read;
    }

    std::shared_ptr<tatami::Matrix<T, IDX> > finish() {
        auto create_matrix = [&](auto& rows, auto& cols) -> auto {
            auto idptrs = compress_sparse_triplets<false>(nrows, ncols, values, rows, cols);

            // Jesus, so much typedefing.
            typedef typename std::remove_reference<decltype(rows)>::type RowType;
            typedef typename std::remove_reference<decltype(idptrs)>::type IndType;
            typedef CompressedSparseColumnMatrix<T, IDX, decltype(values), RowType, IndType> SparseMat;

            return std::shared_ptr<tatami::Matrix<T, IDX> >(new SparseMat(nrows, ncols, std::move(values), std::move(rows), std::move(idptrs), false));
        };

        if (use_short_rows) {
            if (use_short_cols) {
                return create_matrix(short_rows, short_cols);
            } else {
                return create_matrix(short_rows, long_cols);
            }
        } else {
            if (use_short_cols) {
                return create_matrix(long_rows, short_cols);
            } else {
                return create_matrix(long_rows, long_cols);
            }
        }
    }
};
/**
 * @endcond
 */

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * 
 * @return A pointer to a `tatami::Matrix` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This loads a sparse integer matrix from a Matrix Market coordinate file.
 * It will store the data in memory as a compressed sparse column matrix,
 * and is smart enough to use a smaller integer type if the number of rows is less than the maximum value of `uint16_t`.
 * Currently, only unsigned integers are supported.
 */
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix(const char * filepath) {
    TextReader txt(filepath);
    SimpleBuilder<T, IDX> build;
    txt(build);
    return build.finish();
}

#ifdef TATAMI_USE_ZLIB

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored with Gzip compression.
 * @param buffer Size of the buffer to use for decompression, in bytes.
 *
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This is a version of `load_sparse_matrix()` for loading in Gzip-compressed Matrix Market files.
 * To make this function available, make sure to define `TATAMI_USE_ZLIB` and compile with **zlib** support.
 */
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_gzip(const char * filepath, int buffer = 16384) {
    GzipReader unz(filepath, buffer);
    SimpleBuilder<T, IDX> build;
    unz(build);
    return build.finish();
}

#endif

/*******************************************************************************
 *********************** Layered matrix loading functions **********************
 *******************************************************************************/

/**
 * @cond
 */
struct LineAssignments {
    LineAssignments() : rows_per_category(3), lines_per_category(3) {}

    std::vector<uint8_t> category;
    std::vector<size_t> index;
    std::vector<size_t> rows_per_category;
    std::vector<size_t> lines_per_category;
    std::vector<size_t> lines_per_row;
    std::vector<size_t> permutation;

    size_t nrows, ncols, nlines;
    bool passed_preamble = false;

    size_t add(const char* buffer, size_t n) {
        if (buffer[0] == '%') {
            // TODO: should probably check for 'coordinate integer'.
            return read_to_eol(buffer, n);

        } else if (!passed_preamble) {
            passed_preamble = true;
            auto read = process_triplet_line(buffer, nrows, ncols, nlines, n);
            if (read == 0) {
                return 0;
            }

            category.resize(nrows);
            index.resize(nrows);
            permutation.resize(nrows);
            lines_per_row.resize(nrows);
            return read;
        }

        // Assigning each line to a block based on its integer size.
        constexpr int max8 = std::numeric_limits<uint8_t>::max();
        constexpr int max16 = std::numeric_limits<uint16_t>::max();

        int row, col, data;
        auto read = process_triplet_line(buffer, row, col, data, n);
        if (read == 0){
            return 0;
        }
       
        --row; // 1-based.
        if (data > max16) {
            category[row] = std::max(category[row], static_cast<uint8_t>(2));
        } else if (data > max8) {
            category[row] = std::max(category[row], static_cast<uint8_t>(1));
        }
        ++lines_per_row[row];

        return read;
    }

    void finish() {
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

template<typename ROW>
struct LayeredBuilder {
private:
    LineAssignments assign;

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

    bool passed_preamble = false;
public:
    LayeredBuilder(LineAssignments ass) : 
        assign(std::move(ass)),
        row8(assign.lines_per_category[0]),
        col8(assign.lines_per_category[0]),
        dat8(assign.lines_per_category[0]),
        row16(assign.lines_per_category[1]),
        col16(assign.lines_per_category[1]),
        dat16(assign.lines_per_category[1]),
        row32(assign.lines_per_category[2]),
        col32(assign.lines_per_category[2]),
        dat32(assign.lines_per_category[2]) {}

public:
    size_t add(const char* buffer, size_t n) {
        if (buffer[0] == '%') {
            return read_to_eol(buffer, n);
        } else if (!passed_preamble) {
            passed_preamble = true;
            return read_to_eol(buffer, n);
        }

        int row, col, data;
        auto read = process_triplet_line(buffer, row, col, data, n);
        if (read == 0){
            return 0;
        }

        --row; // 1-based.
        --col; 

        auto idx = assign.index[row];
        switch (assign.category[row]) {
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
        
        return read;
    }

private:
    template<typename T, typename IDX, typename U, typename V, typename W>
    std::shared_ptr<Matrix<T, IDX> > create_sparse_matrix(size_t nr, size_t nc, U& values, V& rows, W& cols) {
        auto indptrs = compress_sparse_triplets<false>(nr, nc, values, rows, cols);
        return std::shared_ptr<Matrix<T, IDX> >(new CompressedSparseColumnMatrix<T, IDX, U, V, decltype(indptrs)>(nr, nc, std::move(values), std::move(rows), std::move(indptrs)));
    }

public:
    template<typename T, typename IDX>
    std::shared_ptr<Matrix<T, IDX> > finish() {
        std::vector<std::shared_ptr<Matrix<T, IDX> > > collated;
        if (assign.rows_per_category[0]) {
            collated.push_back(create_sparse_matrix<T, IDX>(assign.rows_per_category[0], assign.ncols, dat8, row8, col8));
        }
        if (assign.rows_per_category[1]) {
            collated.push_back(create_sparse_matrix<T, IDX>(assign.rows_per_category[1], assign.ncols, dat16, row16, col16));
        }
        if (assign.rows_per_category[2]) {
            collated.push_back(create_sparse_matrix<T, IDX>(assign.rows_per_category[2], assign.ncols, dat32, row32, col32));
        }

        if (collated.size() == 0) {
            return create_sparse_matrix<T, IDX>(0, assign.ncols, dat8, row8, col8);
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
template<typename T = double, typename IDX = int, class Functor>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_internal(Functor& process) {
#ifdef TATAMI_PROGRESS_PRINTER
    TATAMI_PROGRESS_PRINTER("tatami::load_layered_sparse_matrix", 1, 3, "Assigning lines to submatrices")
#endif

    LineAssignments ass;
    process(ass);
    ass.finish();

    LayeredMatrixData<T, IDX> output;
    output.permutation = ass.permutation;

#ifdef TATAMI_PROGRESS_PRINTER
    TATAMI_PROGRESS_PRINTER("tatami::load_layered_sparse_matrix", 2, 3, "Building the submatrices")
#endif

    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (ass.nrows > max16) {
        LayeredBuilder<uint16_t> builder(std::move(ass));
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    } else {
        LayeredBuilder<IDX> builder(std::move(ass));
        process(builder);
        output.matrix = builder.template finish<T, IDX>();
    }

#ifdef TATAMI_PROGRESS_PRINTER
    TATAMI_PROGRESS_PRINTER("tatami::load_layered_sparse_matrix", 3, 3, "Done")
#endif

    return output;
}
/**
 * @endcond
 */

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
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
LayeredMatrixData<T, IDX> load_layered_sparse_matrix(const char * filepath) {
    TextReader txt(filepath);
    return load_layered_sparse_matrix_internal(txt);
}

#ifdef TATAMI_USE_ZLIB

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored with Gzip compression.
 * @param buffer Size of the buffer to use for decompression, in bytes.
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
LayeredMatrixData<T, IDX> load_layered_sparse_matrix_gzip(const char * filepath, int buffer = 16384) {
    GzipReader unz(filepath, buffer);
    return load_layered_sparse_matrix_internal(unz);
}

#endif

}

}

#endif
