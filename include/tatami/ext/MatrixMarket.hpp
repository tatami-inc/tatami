#ifndef SCRAN_LOAD_MATRIX_MARKET_HPP
#define SCRAN_LOAD_MATRIX_MARKET_HPP

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

#ifdef TATAMI_USE_ZLIB
#include "zlib.h"
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
template<typename B, typename T>
size_t process_triplet_line(const B* buffer, size_t n, T& arg1, T& arg2, T& arg3, size_t line) {
    auto read = [=](size_t& i, T& arg) -> bool {
        // Chomping up any preceding whitespace.
        while (i < n && std::isspace(buffer[i]) && buffer[i] != '\n') { 
            ++i;
        }
        if (i == n) {
            return false;
        }
        if (buffer[i] == '\n') {
            throw std::runtime_error("premature termination of triplet on line " + std::to_string(line));
        }

        // These had better be non-negative integers, 
        // otherwise this will throw a bunch of errors.
        bool terminated = false;
        arg = 0;
        while (i < n) {
            if (!std::isdigit(buffer[i])) {
                if (std::isspace(buffer[i])) {
                    terminated = true;
                    break;
                } else {
                    throw std::runtime_error("values should be non-negative integers on line " + std::to_string(line));
                }
            }

            arg *= 10;
            arg += (buffer[i] - '0');
            ++i;
        }

        // If it didn't end with some whitespace, there might be more stuff
        // coming, so we return false to restart once the buffer is topped up.
        return terminated;
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

    // Check that we terminate on a newline. If we can't find the newline,
    // we bail out without an error; but if we find something other than
    // a new line, then we need to error.
    while (i < n && std::isspace(buffer[i]) && buffer[i] != '\n') { 
        ++i;
    }
    if (i == n) {
        return 0;
    } else if (buffer[i] != '\n') {
        throw std::runtime_error("triplet should terminate with a newline on " + std::to_string(line));
    }

    return i + 1;
}

template<typename B>
size_t read_to_eol(const B* buffer, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (buffer[i] == '\n' || buffer[i] == '\0') {
            return i + 1;
        }
    }
    return 0; // out of range
}

template<typename B, class Object>
size_t add_from_buffer(B* buffer, size_t length, Object& obj) {
    // Adding whole lines.
    size_t last_processed = 0, total_processed = 0;
    do {
        last_processed = obj.add(buffer + total_processed, length - total_processed);
        total_processed += last_processed;
    } while (last_processed);

    // If nothing was actually stored, this is an error because there's not
    // going to be space to add anything to the buffer. The only exception
    // is the edge case where length == 0 at the EOF (i.e., the previous read
    // caught the last newline and the current read caught the EOF).
    if (total_processed == 0 && length > 0) {
        throw std::runtime_error("failed to read line " + std::to_string(obj.line_number() + 1));
    }

    // Rotating what's left to the front for the next cycle.
    size_t leftovers = length - total_processed;
    for (size_t i = 0; i < leftovers; ++i) {
        buffer[i] = buffer[total_processed + i];
    }

    return leftovers;
}

template<typename B>
void set_terminating_newline(B* buffer, size_t& available) {
    if (buffer[available - 1] != '\n') {
        buffer[available] = '\n';
        ++available;
    }
}

/*******************************************************************************
 ******************************** Reader classes *******************************
 *******************************************************************************/

template<typename B>
struct BufferReader {
    BufferReader(const B* buf, size_t n) : buffer(buf), length(n) {}

    template<class Object>
    void operator()(Object& obj) const {
        size_t last_processed = 0, total_processed = 0;
        do {
            last_processed = obj.add(buffer + total_processed, length - total_processed);
            total_processed += last_processed;
        } while (last_processed);

        size_t leftovers = length - total_processed;
        if (leftovers) {
            std::vector<B> temp(leftovers + 1); // adding a safety newline, if required.
            std::copy(buffer + total_processed, buffer + length, temp.data());
            set_terminating_newline(temp.data(), leftovers);
            obj.add(temp.data(), leftovers);
        }

        return;
    }

private:
    const B* buffer;
    size_t length;
};

struct TextReader {
    TextReader(const char* path, size_t bufsize = 65536) : filepath(path), buffer_size(bufsize) {}

    template<class Object>
    void operator()(Object& obj) const {
        std::ifstream in(filepath);
        if (!in) {
            throw std::runtime_error("failed to open file at '" + std::string(filepath) + "'");
        }

        std::vector<char> buffer(buffer_size + 1); // adding an extra space for the terminating newline.
        size_t leftovers = 0;
        
        while (1) {
            in.read(buffer.data() + leftovers, buffer_size - leftovers);
            size_t available = in.gcount() + leftovers;

            if (!in) {
                set_terminating_newline(buffer.data(), available); // adding a terminating newline, just in case.
                add_from_buffer(buffer.data(), available, obj);
                break;
            } else {
                leftovers = add_from_buffer(buffer.data(), available, obj);
            }
        }

        return;
    }
private:
    const char* filepath;
    size_t buffer_size;
};

#ifdef TATAMI_USE_ZLIB

struct GzipReader {
    GzipReader (const char* path, size_t bufsize = 65536) : filepath(path), buffer_size(bufsize) {}
    const char* filepath;
    const size_t buffer_size;

    struct GZFile {
        GZFile(const char* path) : handle(gzopen(path, "rb")) {
            if (!handle) {
                throw std::runtime_error("failed to open file at '" + std::string(path) + "'");
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
    void operator()(Object& obj) const {
        GZFile gz(filepath);
        std::vector<unsigned char> buffer(buffer_size + 1); // extra space for the terminating newline.

        size_t leftovers = 0;
        while (1) {
            size_t read = gzread(gz.handle, buffer.data() + leftovers, buffer_size - leftovers);
            size_t current_stored = read + leftovers;

            if (read == 0) {
                if (!gzeof(gz.handle)) { 
                    int dummy;
                    throw std::runtime_error(gzerror(gz.handle, &dummy));
                }

                if (current_stored) {
                    // Making sure we have a terminating newline on the last line.
                    set_terminating_newline(buffer.data(), current_stored);
                    add_from_buffer(buffer.data(), current_stored, obj);
                }
                break;
            } else {
                leftovers = add_from_buffer(buffer.data(), current_stored, obj);
            }
        }

        return;
    }
};

// Stolen from 'inf()' at http://www.zlib.net/zpipe.c,
// with some shuffling of code to make it a bit more C++-like.
struct GzipBufferReader {
    GzipBufferReader(const unsigned char* buf, size_t n, size_t bufsize = 65536) : buffer(buf), len(n), buffer_size(bufsize) {}
    const unsigned char* buffer;
    size_t len, buffer_size;

    struct ZStream {
        ZStream() {
            /* allocate inflate state */
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            strm.avail_in = 0;
            strm.next_in = Z_NULL;

            // https://stackoverflow.com/questions/1838699/how-can-i-decompress-a-gzip-stream-with-zlib
            int ret = inflateInit2(&strm, 16+MAX_WBITS); 
            if (ret != Z_OK) {
                throw 1;
            }
        }

        ~ZStream() {
            (void)inflateEnd(&strm);
            return;
        }

        // Delete the remaining constructors.
        ZStream(const ZStream&) = delete;
        ZStream(ZStream&&) = delete;
        ZStream& operator=(const ZStream&) = delete;
        ZStream& operator=(ZStream&&) = delete;

        z_stream strm;
    };

    template<class OBJECT>
    void operator()(OBJECT& obj) {
        std::vector<unsigned char> output(buffer_size + 1); // for a safety newline at EOF, see below.

        ZStream zstr;
        zstr.strm.avail_in = len;
        zstr.strm.next_in = const_cast<unsigned char*>(buffer); // because C interfaces don't have const'ness.

        size_t leftovers = 0;

        /* run inflate() on input until output buffer not full */
        while (1) {
            zstr.strm.avail_out = buffer_size - leftovers;
            zstr.strm.next_out = output.data() + leftovers;
            int ret = inflate(&(zstr.strm), Z_NO_FLUSH);

            switch (ret) {
                case Z_STREAM_ERROR:
                case Z_NEED_DICT:
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    throw std::runtime_error("zlib error");
            }

            size_t current_stored = buffer_size - zstr.strm.avail_out;
            
            // Making sure we have a terminating newline.
            if (ret == Z_STREAM_END) {
                if (current_stored) {
                    set_terminating_newline(output.data(), current_stored);
                    add_from_buffer(output.data(), current_stored, obj);
                }
                break;
            } else {
                leftovers = add_from_buffer(output.data(), current_stored, obj);
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
    size_t current_line = 0;

    size_t nrows, ncols, nlines;
    bool passed_preamble = false;

    template<typename B>
    size_t add(const B* buffer, size_t n) {
        if (buffer[0] == '%') {
            // TODO: should probably check for 'coordinate integer'.
            return read_to_eol(buffer, n);

        } else if (!passed_preamble) {
            auto read = process_triplet_line(buffer, n, nrows, ncols, nlines, current_line);
            if (read == 0) {
                return 0;
            }
            passed_preamble = true;

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
        auto read = process_triplet_line(buffer, n, row, col, data, current_line);
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

    size_t line_number() const {
        return current_line;
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
 * @param bufsize Size of the buffer (in bytes) to use when reading from file.
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
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix(const char * filepath, size_t bufsize = 65536) {
    TextReader txt(filepath, bufsize);
    SimpleBuilder<T, IDX> build;
    txt(build);
    return build.finish();
}

/**
 * @param buffer Array containing the contents of a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param n Length of the array.
 * 
 * @return A pointer to a `tatami::Matrix` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 * @tparam B Type of the buffer, usually `char` or `unsigned char`.
 *
 * This is equivalent to `load_sparse_matrix()` but assumes that the entire file has been read into `buffer`.
 */
template<typename T = double, typename IDX = int, typename B>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_from_buffer(const B* buffer, size_t n) {
    BufferReader raw(buffer, n);
    SimpleBuilder<T, IDX> build;
    raw(build);
    return build.finish();
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
 * This is a version of `load_sparse_matrix()` for loading in Gzip-compressed Matrix Market files.
 * To make this function available, make sure to define `TATAMI_USE_ZLIB` and compile with **zlib** support.
 */
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_gzip(const char * filepath, size_t bufsize = 65536) {
    GzipReader unz(filepath, bufsize);
    SimpleBuilder<T, IDX> build;
    unz(build);
    return build.finish();
}

/**
 * @param buffer Array containing the contents of a Matrix Market file.
 * The file should contain integer data in the coordinate format, stored in text without any compression.
 * @param n Size of the `buffer` array.
 * @param bufsize Size of the buffer to use for decompression, in bytes.
 *
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This is equivalent to `load_sparse_matrix_gzip()` but assumes that the entire file has been read into `buffer`.
 */
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_from_buffer_gzip(const unsigned char * buffer, size_t n, size_t bufsize = 65536) {
    GzipBufferReader unz(buffer, n, bufsize);
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
    size_t current_line = 0;

    template<typename B>
    size_t add(const B* buffer, size_t n) {
        if (buffer[0] == '%') {
            // TODO: should probably check for 'coordinate integer'.
            return read_to_eol(buffer, n);

        } else if (!passed_preamble) {
            auto read = process_triplet_line(buffer, n, nrows, ncols, nlines, current_line);
            if (read == 0) {
                return 0;
            }
            passed_preamble = true;

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
        auto read = process_triplet_line(buffer, n, row, col, data, current_line);
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
        ++current_line;

        return read;
    }

    size_t line_number() const {
        return current_line;
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
    size_t current_line = 0;
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
    template<typename B>
    size_t add(const B* buffer, size_t n) {
        if (buffer[0] == '%') {
            return read_to_eol(buffer, n);
        } else if (!passed_preamble) {
            size_t read = read_to_eol(buffer, n);
            if (read) {
                passed_preamble = true;
            }
            return read;
        }

        int row, col, data;
        auto read = process_triplet_line(buffer, n, row, col, data, current_line);
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

        ++current_line;
        return read;
    }

    size_t line_number() const {
        return current_line;
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
    TextReader txt(filepath, bufsize);
    return load_layered_sparse_matrix_internal(txt);
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
    BufferReader raw(buffer, n);
    return load_layered_sparse_matrix_internal(raw);
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
    GzipReader unz(filepath, bufsize);
    return load_layered_sparse_matrix_internal(unz);
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
    GzipBufferReader raw(buffer, n, bufsize);
    return load_layered_sparse_matrix_internal(raw);
}

#endif

}

}

#endif
