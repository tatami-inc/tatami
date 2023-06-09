#ifndef TATAMI_MATRIX_MARKET_SIMPLE_HPP
#define TATAMI_MATRIX_MARKET_SIMPLE_HPP

#include <limits>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <cctype>
#include <string>

#include "../../sparse/CompressedSparseMatrix.hpp"
#include "../../utils/compress_sparse_triplets.hpp"

#include "byteme/RawFileReader.hpp"
#include "byteme/RawBufferReader.hpp"

#ifdef TATAMI_USE_ZLIB
#include "byteme/GzipFileReader.hpp"
#include "byteme/ZlibBufferReader.hpp"
#include "byteme/SomeFileReader.hpp"
#include "byteme/SomeBufferReader.hpp"
#endif

/**
 * @file simple.hpp
 *
 * @brief Read a sparse non-negative integer matrix from the Matrix Market coordinate format.
 */

namespace tatami {

/**
 * @namespace tatami::MatrixMarket
 * @brief Classes and utilites to handle Matrix Market files.
 */
namespace MatrixMarket {

/**
 * @cond
 */
struct BaseMMParser {
private:
    size_t current_line = 0;
    size_t current_data_line = 0;
    bool passed_preamble = false;
    bool in_comment = false;

    int onto = 0;
    bool non_empty = false;

    size_t currow = 0, curcol = 0, curval = 0;
    size_t nrows, ncols, nlines;

public:
    // Some SFINAE nonsense for fun and profit.
    template<class X, typename = int>
    struct preamble_only {
        constexpr static bool value = false;
    };

    template<class X>
    struct preamble_only<X, decltype((void) X::preamble_only, 0)> {
        constexpr static bool value = X::preamble_only;
    };

private:
    template<class Store>
    void new_line(Store& store) {
        if (in_comment) {
            in_comment = false;
        } else {
            if ((onto == 3 && !non_empty) || (onto == 2 && non_empty)) {
                // i.e., there are three fields.
                ;
            } else {
                throw std::runtime_error("line " + std::to_string(current_line + 1) + " should contain three values");
            }

            if (!passed_preamble) {
                nrows = currow;
                ncols = curcol;
                nlines = curval;

                store.setdim(currow, curcol, curval);
                passed_preamble = true;
            } else {
                if (!currow) {
                    throw std::runtime_error("row index must be positive on line " + std::to_string(current_line + 1));
                }
                if (currow > nrows) {
                    throw std::runtime_error("row index out of range on line " + std::to_string(current_line + 1));
                }

                if (!curcol) {
                    throw std::runtime_error("column index must be positive on line " + std::to_string(current_line + 1));
                }
                if (curcol > ncols) {
                    throw std::runtime_error("column index out of range on line " + std::to_string(current_line + 1));
                }

                if (current_data_line >= nlines) {
                    throw std::runtime_error("more lines present than specified in the header (" + std::to_string(nlines) + ")");
                }

                if constexpr(!preamble_only<Store>::value) {
                    store.addline(currow - 1, curcol - 1, curval, current_data_line);
                }
                ++current_data_line;
            }

            onto = 0;
            currow = 0;
            curcol = 0;
            curval = 0;
            non_empty = false;
        }
        ++current_line;
    }

public:
    template<class Reader, class Store>
    void operator()(Reader& reader, Store& store) {
        bool remaining;

        do {
            remaining = reader();
            auto buffer = reinterpret_cast<const char*>(reader.buffer());
            auto n = reader.available();

            size_t i = 0;
            while (i < n) {
                if (buffer[i] == '%') {
                    in_comment = true;

                    // Try to get quickly to the next line.
                    do {
                        ++i;
                    } while (i < n && buffer[i] != '\n');
                } else if (buffer[i] == '\n') {
                    new_line(store);
                    ++i;

                    // Breaking out early if we only care about the preamble.
                    if constexpr(preamble_only<Store>::value) {
                        if (passed_preamble) {
                            return;
                        }
                    }
                } else if (in_comment) {
                    // Try to get to the next line, again.
                    do {
                        ++i;
                    } while (i < n && buffer[i] != '\n');
                } else {
                    while (i < n) {
                        if (!std::isdigit(buffer[i])) {
                            if (!std::isspace(buffer[i])) {
                                throw std::runtime_error("values should be non-negative integers on line " + std::to_string(current_line));
                            }

                            if (buffer[i] == '\n') {
                                break;
                            }

                            if (non_empty) {
                                ++onto;
                                non_empty = false;
                            }
                            ++i;
                            continue;
                        }

                        non_empty = true; // check that something was _actually_ present in this field.
                        int delta = buffer[i] - '0';
                        switch(onto) {
                            case 0:
                                currow *= 10;
                                currow += delta;
                                break;
                            case 1:
                                curcol *= 10;
                                curcol += delta;
                                break;
                            case 2:
                                curval *= 10;
                                curval += delta;
                                break;
                        }
                        ++i;
                    }
                }
            }
        } while (remaining);

        // If onto = 0 and non_empty = false, we ended on a newline, so 
        // there's no extra entry to add. Otherwise, we try to add the 
        // last line that was _not_ terminated by a newline.
        if (onto != 0 || non_empty) { 
            new_line(store); 
        }

        if (!passed_preamble) {
            throw std::runtime_error("no header line specifying the dimensions");
        }

        if (current_data_line != nlines) {
            throw std::runtime_error("detected " + std::to_string(current_data_line) + " lines but " + std::to_string(nlines) + " lines specified in the header");
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
template<typename T, typename IDX> 
struct SimpleBuilder {
    struct Core {
        std::vector<uint16_t> short_rows;
        std::vector<IDX> long_rows;
        bool use_short_rows = false;

        std::vector<uint16_t> short_cols;
        std::vector<IDX> long_cols;
        bool use_short_cols = false;

        size_t nrows, ncols;
        std::vector<T> values;

        void setdim(size_t currow, size_t curcol, size_t nlines) {
            nrows = currow;
            ncols = curcol;

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
        }

        void addline(size_t row, size_t col, size_t val, size_t line) {
            if (use_short_rows) {
                short_rows[line] = row;
            } else {
                long_rows[line] = row;
            }

            if (use_short_cols) {
                short_cols[line] = col;
            } else {
                long_cols[line] = col;
            }

            values[line] = val;
        }
    };

    template<class Reader>
    static std::shared_ptr<tatami::Matrix<T, IDX> > build(Reader& reader) {
        BaseMMParser parser;
        Core store;
        parser(reader, store);

        auto create_matrix = [&](auto& rows, auto& cols) -> auto {
            auto idptrs = compress_sparse_triplets<false>(store.nrows, store.ncols, store.values, rows, cols);

            // Jesus, so much typedefing.
            typedef typename std::remove_reference<decltype(rows)>::type RowType;
            typedef typename std::remove_reference<decltype(idptrs)>::type IndType;
            typedef CompressedSparseColumnMatrix<T, IDX, decltype(store.values), RowType, IndType> SparseMat;

            return std::shared_ptr<tatami::Matrix<T, IDX> >(
                new SparseMat(
                    store.nrows, 
                    store.ncols, 
                    std::move(store.values), 
                    std::move(rows), 
                    std::move(idptrs), 
                    false
                )
            );
        };

        if (store.use_short_rows) {
            if (store.use_short_cols) {
                return create_matrix(store.short_rows, store.short_cols);
            } else {
                return create_matrix(store.short_rows, store.long_cols);
            }
        } else {
            if (store.use_short_cols) {
                return create_matrix(store.long_rows, store.short_cols);
            } else {
                return create_matrix(store.long_rows, store.long_cols);
            }
        }
    }
};
/**
 * @endcond
 */

/**
 * @cond
 */
template<class Builder>
auto operate_on_file(const char * filepath, int compression, size_t bufsize) {
    if (compression != 0) {
#ifndef TATAMI_USE_ZLIB
        throw std::runtime_error("tatami not compiled with support for non-zero 'compression'");
#else
        if (compression == -1) {
            byteme::SomeFileReader reader(filepath, bufsize);
            return Builder::build(reader);
       } else if (compression == 1) {
            byteme::GzipFileReader reader(filepath, bufsize);
            return Builder::build(reader);
        }
#endif
    }
    byteme::RawFileReader reader(filepath, bufsize);
    return Builder::build(reader);
}

template<class Builder>
auto operate_on_buffer(const unsigned char * buffer, size_t n, int compression, size_t bufsize) {
    if (compression != 0) {
#ifndef TATAMI_USE_ZLIB
        throw std::runtime_error("tatami not compiled with support for non-zero 'compression'");
#else
        if (compression == -1) {
            byteme::SomeBufferReader reader(buffer, n, bufsize);
            return Builder::build(reader);
        } else if (compression == 1) {
            byteme::ZlibBufferReader reader(buffer, n, 3, bufsize);
            return Builder::build(reader);
        }
#endif
    }

    byteme::RawBufferReader reader(buffer, n);
    return Builder::build(reader);
}
/**
 * @endcond
 */

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain non-negative integer data in the coordinate format.
 * @param compression Compression method for the file - no compression (0) or Gzip compression (1).
 * If set to -1, the function will automatically guess the compression based on magic numbers.
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
 * 
 * To support Gzip decompression, make sure to define `TATAMI_USE_ZLIB` and compile with **zlib** support.
 */
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_from_file(const char * filepath, int compression = 0, size_t bufsize = 65536) {
    return operate_on_file<SimpleBuilder<T, IDX> >(filepath, compression, bufsize);
}

// For back-compatibility.
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix(const char * filepath, size_t bufsize = 65536) {
    return load_sparse_matrix_from_file(filepath, 0, bufsize);
}

/**
 * @param buffer Array containing the contents of a Matrix Market file.
 * The file should contain non-negative integer data in the coordinate format.
 * @param n Length of the array.
 * @param compression Compression method for the file contents - no compression (0) or Gzip/Zlib compression (1).
 * If set to -1, the function will automatically guess the compression based on magic numbers.
 * @param bufsize Size of the buffer (in bytes) to use when decompressing the file contents.
 * 
 * @return A pointer to a `tatami::Matrix` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This is equivalent to `load_sparse_matrix_from_file()` but assumes that the entire file has been read into `buffer`.
 * To support Gzip decompression, make sure to define `TATAMI_USE_ZLIB` and compile with **zlib** support.
 */
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_from_buffer(const unsigned char* buffer, size_t n, int compression = 0, size_t bufsize = 65536) {
    return operate_on_buffer<SimpleBuilder<T, IDX> >(buffer, n, compression, bufsize);
}

#ifdef TATAMI_USE_ZLIB

// For back-compatibility.
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_gzip(const char * filepath, size_t bufsize = 65536) {
    return load_sparse_matrix_from_file(filepath, 1, bufsize);
}

// For back-compatibility.
template<typename T = double, typename IDX = int>
std::shared_ptr<tatami::Matrix<T, IDX> > load_sparse_matrix_from_buffer_gzip(const unsigned char * buffer, size_t n, size_t bufsize = 65536) {
    return load_sparse_matrix_from_buffer(buffer, n, 1, bufsize);
}

#endif

/**
 * @brief Details extracted from a Matrix Market header.
 */
struct HeaderDetails {
    /**
     * Number of rows.
     */
    size_t nrow;

    /**
     * Number of columns.
     */
    size_t ncol;

    /**
     * Number of lines.
     */
   size_t nlines;
};

/**
 * @cond
 */
struct Inspector {
    struct Core {
        HeaderDetails header;
        constexpr static bool preamble_only = true;
        void setdim(size_t nr, size_t nc, size_t nl) {
            header.nrow = nr;
            header.ncol = nc;
            header.nlines = nl;
        }
    };

    template<class Reader>
    static HeaderDetails build(Reader& reader) {
        BaseMMParser parser;
        Core store;
        parser(reader, store);
        return store.header;
    }
};
/**
 * @endcond
 */

/**
 * @param filepath Path to a Matrix Market file.
 * The file should contain non-negative integer data in the coordinate format.
 * @param compression Compression method for the file - no compression (0) or Gzip compression (1).
 * If set to -1, the function will automatically guess the compression based on magic numbers.
 * @param bufsize Size of the buffer (in bytes) to use when reading from file. 
 * 
 * @return A `HeaderDetails` object containing details about the file contents.
 */
inline HeaderDetails extract_header_from_file(const char* filepath, int compression = 0, size_t bufsize = 65536) {
    return operate_on_file<Inspector>(filepath, compression, bufsize);
}

/**
 * @param buffer Array containing the contents of a Matrix Market file.
 * The file should contain non-negative integer data in the coordinate format.
 * @param n Length of the array.
 * @param compression Compression method for the file contents - no compression (0) or Gzip/Zlib compression (1).
 * If set to -1, the function will automatically guess the compression based on magic numbers.
 * @param bufsize Size of the buffer (in bytes) to use when decompressing the file contents.
 * 
 * @return A `HeaderDetails` object containing details about the file contents.
 */
inline HeaderDetails extract_header_from_buffer(const unsigned char* buffer, size_t n, int compression = 0, size_t bufsize = 65536) {
    return operate_on_buffer<Inspector>(buffer, n, compression, bufsize);
}

}

}

#endif
