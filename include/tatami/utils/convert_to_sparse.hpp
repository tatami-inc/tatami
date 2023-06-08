#ifndef TATAMI_CONVERT_TO_SPARSE_H
#define TATAMI_CONVERT_TO_SPARSE_H

#include "../sparse/CompressedSparseMatrix.hpp"
#include "../stats/utils.hpp"

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
 * @param threads Number of threads to use.
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
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const InputMatrix_* incoming, int threads = 1) {
    Index_ NR = incoming->nrow();
    Index_ NC = incoming->ncol();
    Index_ primary = (row_ ? NR : NC);
    Index_ secondary = (row_ ? NC : NR);

    std::vector<std::vector<StoredValue_> > store_v(primary);
    std::vector<std::vector<StoredValue_> > store_i(primary);

    if (threads == 0) {
        threads = 1; // for back-compatibility with calls using the old 'reserve' argument.
    }

    typedef typename InputMatrix_::value_type InputData_; 
    typedef typename InputMatrix_::index_type InputIndex_; 

    if (row_ == incoming->prefer_rows()) {
        if (incoming->sparse()) {
            // Using InputIndex_ here, otherwise template type deduction gets confused for consecutive_extractor.
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<row_, true, InputData_, InputIndex_>(incoming, 0, primary);

                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    auto range = wrk->fetch(p, buffer_v.data(), buffer_i.data());
                    auto& sv = store_v[p];
                    auto& si = store_i[p];
                    sv.reserve(range.number);
                    si.reserve(range.number);

                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        if (*range.value) {
                            sv.push_back(*range.value);
                            si.push_back(*range.index);
                        }
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false, InputData_, InputIndex_>(incoming, 0, primary);

                // Special conversion from dense to save ourselves from having to make
                // indices that we aren't really interested in.
                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    auto& sv = store_v[p];
                    auto& si = store_i[p];

                    for (InputIndex_ s = 0; s < secondary; ++s, ++ptr) {
                        if (*ptr) {
                            sv.push_back(*ptr);
                            si.push_back(s);
                        }
                    }
                }
            }, primary, threads);
        }

    } else {
        // We iterate on the incoming matrix's preferred dimension, under the
        // assumption that it may be arbitrarily costly to extract in the
        // non-preferred dim; it is thus cheaper to do cache-unfriendly inserts
        // into the output buffers. 

        if (incoming->sparse()) {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(primary);
                std::vector<InputIndex_> buffer_i(primary);
                auto wrk = consecutive_extractor<!row_, true, InputData_, InputIndex_>(incoming, 0, secondary, start, length);

                for (InputIndex_ s = 0; s < secondary; ++s) {
                    auto range = wrk->fetch(s, buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        if (*range.value) {
                            store_v[*range.index].push_back(*range.value);
                            store_i[*range.index].push_back(s);
                        }
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                auto wrk = consecutive_extractor<!row_, false, InputData_, InputIndex_>(incoming, 0, secondary, start, length);
                auto len = wrk->block_length;
                std::vector<InputData_> buffer_v(len);

                for (InputIndex_ s = 0; s < secondary; ++s) {
                    auto ptr = wrk->fetch(s, buffer_v.data());
                    for (InputIndex_ p = 0; p < len; ++p, ++ptr) {
                        if (*ptr) {
                            store_v[p + start].push_back(*ptr);
                            store_i[p + start].push_back(s);
                        }
                    }
                }
            }, primary, threads);
        }
    }

    // Concatenating everything together.
    size_t total_size = 0;
    std::vector<size_t> indptrs(primary + 1);
    for (size_t p = 0; p < primary; ++p) {
        total_size += store_v[p].size();
        indptrs[p + 1] = total_size;
    }

    std::vector<StoredValue_> output_v;
    output_v.reserve(total_size);
    std::vector<StoredIndex_> output_i;
    output_i.reserve(total_size);

    for (size_t p = 0; p < primary; ++p) {
        output_v.insert(output_v.end(), store_v[p].begin(), store_v[p].end());
        output_i.insert(output_i.end(), store_i[p].begin(), store_i[p].end());
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
 * @param threads Number of threads to use.
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
std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const InputMatrix_* incoming, int order, int threads = 1) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows());
    }
    if (order == 0) {
        return convert_to_sparse<true, Value_, Index_, StoredValue_, StoredIndex_, InputMatrix_>(incoming, threads);
    } else {
        return convert_to_sparse<false, Value_, Index_, StoredValue_, StoredIndex_, InputMatrix_>(incoming, threads);
    }
}

}

#endif
