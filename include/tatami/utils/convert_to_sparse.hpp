#ifndef TATAMI_CONVERT_TO_SPARSE_H
#define TATAMI_CONVERT_TO_SPARSE_H

#include "../base/sparse/CompressedSparseMatrix.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_sparse.hpp
 *
 * @brief Convert a matrix into a compressed sparse format.
 */

namespace tatami {

/**
 * @tparam row_ Whether to return a compressed sparse row matrix.
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputMatrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`, possibly containing delayed operations.
 * @param reserve The expected density of non-zero values in `incoming`.
 * A slight overestimate will avoid reallocation of the temporary vectors.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row_ = true`, the matrix is compressed sparse row, otherwise it is compressed sparse column.
 */
template <bool row_, 
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    class InputMatrix_
>
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const InputMatrix_* incoming, double reserve = 0.1) {
    Index_ NR = incoming->nrow();
    Index_ NC = incoming->ncol();
    Index_ primary = (row_ ? NR : NC);
    Index_ secondary = (row_ ? NC : NR);

    std::vector<StoredValue_> output_v;
    std::vector<StoredIndex_> output_i;
    std::vector<size_t> indptrs(primary + 1);

    typedef typename InputMatrix_::value_type InputData_; 
    typedef typename InputMatrix_::index_type InputIndex_; 

    if (row_ == incoming->prefer_rows()) {
        size_t reservation = static_cast<double>(NR) * static_cast<double>(NC) * reserve;
        output_v.reserve(reservation);
        output_i.reserve(reservation);
        std::vector<InputData_> buffer_v(secondary);

        if (incoming->sparse()) {
            std::vector<InputIndex_> buffer_i(secondary);
            auto wrk = new_extractor<row_, true>(incoming);
            for (size_t p = 0; p < primary; ++p) {
                auto range = wrk->fetch(p, buffer_v.data(), buffer_i.data());
                for (size_t i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                    if (*range.value) {
                        output_v.push_back(*range.value);
                        output_i.push_back(*range.index);
                    }
                }
                indptrs[p + 1] = output_v.size();
            }

        } else {
            auto wrk = new_extractor<row_, false>(incoming);

            // Special conversion from dense to save ourselves from having to make
            // indices that we aren't really interested in.
            for (size_t p = 0; p < primary; ++p) {
                auto ptr = wrk->fetch(p, buffer_v.data());
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
        std::vector<std::vector<StoredValue_> > store_v(primary);
        std::vector<std::vector<StoredValue_> > store_i(primary);
        size_t reservation = secondary * reserve;
        for (size_t p = 0; p < primary; ++p) {
            store_v[p].reserve(reservation);
            store_i[p].reserve(reservation);
        }
        std::vector<InputData_> buffer_v(primary);

        if (incoming->sparse()) {
            auto wrk = new_extractor<!row_, true>(incoming);
            std::vector<InputIndex_> buffer_i(primary);
            for (size_t s = 0; s < secondary; ++s) {
                auto range = wrk->fetch(s, buffer_v.data(), buffer_i.data());
                for (size_t i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                    if (*range.value) {
                        store_v[*range.index].push_back(*range.value);
                        store_i[*range.index].push_back(s);
                    }
                }
            }

        } else {
            auto wrk = new_extractor<!row_, false>(incoming);
            for (size_t s = 0; s < secondary; ++s) {
                auto ptr = wrk->fetch(s, buffer_v.data());
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

    return std::shared_ptr<Matrix<Value_, Index_> >(
        new CompressedSparseMatrix<
            row_, 
            Value_, 
            Index_,
            decltype(output_v),
            decltype(output_i),
            decltype(indptrs)
        >(
            NR, 
            NC, 
            std::move(output_v), 
            std::move(output_i), 
            std::move(indptrs)
        )
    );
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputMatrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output matrix - compressed sparse row (0) or column (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    class InputMatrix_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const InputMatrix_* incoming, int order) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows());
    }
    if (order == 0) {
        return convert_to_sparse<true, Value_, Index_, StoredValue_, StoredIndex_, InputMatrix_>(incoming);
    } else {
        return convert_to_sparse<false, Value_, Index_, StoredValue_, StoredIndex_, InputMatrix_>(incoming);
    }
}

}

#endif
