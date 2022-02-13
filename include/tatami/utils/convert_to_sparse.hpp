#ifndef TATAMI_CONVERT_TO_SPARSE_H
#define TATAMI_CONVERT_TO_SPARSE_H

#include "../base/CompressedSparseMatrix.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_sparse.hpp
 *
 * Convert a matrix into a compressed sparse format.
 */

namespace tatami {

/**
 * @tparam row_ Whether to return a compressed sparse row matrix.
 * @tparam MatrixIn Input matrix class, most typically a `tatami::Matrix`.
 * @tparam Data Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`, possibly containing delayed operations.
 * @param reserve The expected density of non-zero values in `incoming`.
 * A slight overestimate will avoid reallocation of the temporary vectors.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row = true`, the matrix is compressed sparse row, otherwise it is compressed sparse column.
 */
template <bool row_, class MatrixIn, typename DataOut = typename MatrixIn::data_type, typename IndexOut = typename MatrixIn::index_type>
inline std::shared_ptr<Matrix<DataOut, IndexOut> > convert_to_sparse(const MatrixIn* incoming, double reserve = 0.1) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    size_t primary = (row_ ? NR : NC);
    size_t secondary = (row_ ? NC : NR);

    std::vector<DataOut> output_v;
    std::vector<IndexOut> output_i;
    std::vector<size_t> indptrs(primary + 1);

    typedef typename MatrixIn::data_type DataIn; 
    typedef typename MatrixIn::index_type IndexIn; 

    if (row_ == incoming->prefer_rows()) {
        size_t reservation = static_cast<double>(NR * NC) * reserve;
        output_v.reserve(reservation);
        output_i.reserve(reservation);

        auto wrk = incoming->new_workspace(row_);
        std::vector<DataIn> buffer_v(secondary);

        if (incoming->sparse()) {
            std::vector<IndexIn> buffer_i(secondary);

            for (size_t p = 0; p < primary; ++p) {
                SparseRange<DataIn, IndexIn> range;
                if constexpr(row_) {
                    range = incoming->sparse_row(p, buffer_v.data(), buffer_i.data(), wrk.get());
                } else {
                    range = incoming->sparse_column(p, buffer_v.data(), buffer_i.data(), wrk.get());
                }

                for (size_t i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                    if (*range.value) {
                        output_v.push_back(*range.value);
                        output_i.push_back(*range.index);
                    }
                }
                indptrs[p + 1] = output_v.size();
            }

        } else {
            // Special conversion from dense to save ourselves from having to make
            // indices that we aren't really interested in.
            for (size_t p = 0; p < primary; ++p) {
                const DataIn* ptr;
                if constexpr(row_) {
                    ptr = incoming->row(p, buffer_v.data(), wrk.get());
                } else {
                    ptr = incoming->column(p, buffer_v.data(), wrk.get());
                }

                for (size_t s = 0; s < secondary; ++s, ++ptr) {
                    if (*ptr) {
                        output_v.push_back(*ptr);
                        output_i.push_back(s);
                    }
                }
                indptrs[p + 1] = output_v.size();
            }
        }

    } else {
        // We iterate on the incoming matrix's preferred dimension, under the
        // assumption that it may be arbitrarily costly to extract in the
        // non-preferred dim; it is thus cheaper to do cache-unfriendly inserts
        // into the output buffer. In this case, there's not much choice but to
        // make extensible vectors for each primary dimension.
        std::vector<std::vector<DataOut> > store_v(primary);
        std::vector<std::vector<DataOut> > store_i(primary);
        size_t reservation = secondary * reserve;
        for (size_t p = 0; p < primary; ++p) {
            store_v[p].reserve(reservation);
            store_i[p].reserve(reservation);
        }

        auto wrk = incoming->new_workspace(!row_);
        std::vector<DataIn> buffer_v(primary);

        if (incoming->sparse()) {
            std::vector<IndexIn> buffer_i(primary);
            for (size_t s = 0; s < secondary; ++s) {
                SparseRange<DataIn, IndexIn> range;
                if constexpr(row_) {
                    range = incoming->sparse_column(s, buffer_v.data(), buffer_i.data(), wrk.get());
                } else {
                    range = incoming->sparse_row(s, buffer_v.data(), buffer_i.data(), wrk.get());
                }

                for (size_t i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                    if (*range.value) {
                        store_v[*range.index].push_back(*range.value);
                        store_i[*range.index].push_back(s);
                    }
                }
            }
        } else {
            for (size_t s = 0; s < secondary; ++s) {
                const DataIn* ptr;
                if constexpr(row_) {
                    ptr = incoming->column(s, buffer_v.data(), wrk.get());
                } else {
                    ptr = incoming->row(s, buffer_v.data(), wrk.get());
                }

                for (size_t p = 0; p < primary; ++p, ++ptr) {
                    if (*ptr) {
                        store_v[p].push_back(*ptr);
                        store_i[p].push_back(s);
                    }
                }
            }
        }

        // Concatenating everything together.
        size_t total_size = 0;
        for (size_t p = 0; p < primary; ++p) {
            total_size += store_v[p].size();
            indptrs[p + 1] = total_size;
        }

        output_v.reserve(total_size);
        output_i.reserve(total_size);
        for (size_t p = 0; p < primary; ++p) {
            output_v.insert(output_v.end(), store_v[p].begin(), store_v[p].end());
            output_i.insert(output_i.end(), store_i[p].begin(), store_i[p].end());
        }
    }

    return std::shared_ptr<Matrix<DataOut, IndexOut> >(new CompressedSparseMatrix<row_, DataOut, IndexOut>(NR, NC, std::move(output_v), std::move(output_i), std::move(indptrs)));
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam MatrixIn Input matrix class, most typically a `tatami::Matrix`.
 * @tparam Data Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 * @tparam T Type of the values in the matrix.
 * @tparam IDX Type of index values.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output matrix - compressed sparse row (0) or column (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <class MatrixIn, typename DataOut = typename MatrixIn::data_type, typename IndexOut = typename MatrixIn::index_type>
std::shared_ptr<Matrix<DataOut, IndexOut> > convert_to_sparse(const MatrixIn* incoming, int order) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows());
    }
    if (order == 0) {
        return convert_to_sparse<true, MatrixIn, DataOut, IndexOut>(incoming);
    } else {
        return convert_to_sparse<false, MatrixIn, DataOut, IndexOut>(incoming);
    }
}

}

#endif
