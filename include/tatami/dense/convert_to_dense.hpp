#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "./DenseMatrix.hpp"
#include "./transpose.hpp"

#include "../utils/consecutive_extractor.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/copy.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_dense.hpp
 *
 * @brief Convert a matrix into a dense format.
 */

namespace tatami {

/**
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam InputValue_ Type of data values in the input.
 * @tparam InputIndex_ Integer type for the indices in the input.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @tparam row Whether to return a row-major matrix.
 * @param[out] store Pointer to an array of length equal to the product of the dimensions of `incoming`.
 * On output, this is filled with values from `incoming` in row- or column-major format depending on `row`.
 * @param threads Number of threads to use.
 */
template <typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense(const Matrix<InputValue_, InputIndex_>* incoming, bool row, StoredValue_* store, int threads = 1) {
    InputIndex_ NR = incoming->nrow();
    InputIndex_ NC = incoming->ncol();
    bool pref_rows = incoming->prefer_rows();
    size_t primary = (pref_rows ? NR : NC);
    size_t secondary = (pref_rows ? NC : NR);

    if (row == pref_rows) {
        constexpr bool same_type = std::is_same<InputValue_, StoredValue_>::value;
        parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
            std::vector<InputValue_> temp(same_type ? 0 : secondary);
            auto store_copy = store + static_cast<size_t>(start) * secondary; // cast to size_t to avoid overflow.
            auto wrk = consecutive_extractor<false>(incoming, row, start, length);

            for (InputIndex_ x = 0; x < length; ++x) {
                if constexpr(same_type) {
                    auto ptr = wrk->fetch(store_copy);
                    copy_n(ptr, secondary, store_copy);
                } else {
                    auto ptr = wrk->fetch(temp.data());
                    std::copy_n(ptr, secondary, store_copy);
                }
                store_copy += secondary;
            }
        }, primary, threads);

    } else if (incoming->sparse()) {
        std::fill_n(store, primary * secondary, 0); // already cast to size_t to avoid overflow.

        // We iterate over the input matrix's preferred dimension but split
        // into threads along the non-preferred dimension. This aims to
        // reduce false sharing across threads during writes, as locations
        // for simultaneous writes in the transposed matrix will be
        // separated by around 'secondary * length' elements. 
        parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
            auto store_copy = store;

            auto wrk = consecutive_extractor<true, InputValue_, InputIndex_>(incoming, pref_rows, 0, primary, start, length);
            std::vector<InputValue_> vtemp(length);
            std::vector<InputIndex_> itemp(length);

            // Note that we don't use the blocked transposition strategy
            // from the dense case, because the overhead of looping is 
            // worse than the cache misses for sparse data.
            for (size_t x = 0; x < primary; ++x) {
                auto range = wrk->fetch(vtemp.data(), itemp.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    store_copy[static_cast<size_t>(range.index[i]) * primary] = range.value[i]; // again, casting here.
                }
                ++store_copy;
            }
        }, secondary, threads);

    } else {
        // Same logic as described for the sparse case; we iterate along the
        // preferred dimension but split into threads along the non-preferred
        // dimension to reduce false sharing.
        parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
            auto store_copy = store + static_cast<size_t>(start) * primary; // cast to size_t to avoid overflow.

            auto wrk = consecutive_extractor<false, InputValue_, InputIndex_>(incoming, pref_rows, 0, primary, start, length);
            const size_t length_as_size_t = length;

            // Performing a blocked transposition to be more
            // cache-friendly. This involves collecting several
            // consecutive primary dimension elements so that we can
            // transpose by blocks along the secondary dimension.
            constexpr size_t block_size = 16;
            const size_t alloc = std::min(primary, block_size);
            std::vector<InputValue_> bigbuffer(length_as_size_t * alloc);
            std::vector<const InputValue_*> ptrs(alloc);
            std::vector<InputValue_*> buf_ptrs;
            buf_ptrs.reserve(alloc);
            auto first = bigbuffer.data();
            for (size_t i = 0; i < alloc; ++i, first += length_as_size_t) {
                buf_ptrs.push_back(first);
            }

            size_t prim_i = 0;
            while (prim_i < primary) {
                size_t prim_to_process = std::min(primary - prim_i, block_size);
                for (size_t c = 0; c < prim_to_process; ++c) {
                    ptrs[c] = wrk->fetch(buf_ptrs[c]);
                }

                size_t sec_i = 0;
                while (sec_i < length_as_size_t) {
                    size_t sec_to_process = std::min(block_size, length_as_size_t - sec_i);
                    auto output = store_copy + sec_i * primary;
                    for (size_t c = 0; c < prim_to_process; ++c, ++output) {
                        auto input = ptrs[c] + sec_i;
                        auto output2 = output;
                        for (size_t r = 0; r < sec_to_process; ++r, ++input, output2 += primary) {
                            *output2 = *input;
                        }
                    }
                    sec_i += sec_to_process;
                }

                prim_i += prim_to_process;
                store_copy += prim_to_process;
            }
        }, secondary, threads);
    }

    return;
}

/**
 * @tparam row Whether to return a row-major matrix.
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam InputValue_ Type of data values in the input.
 * @tparam InputIndex_ Integer type for the indices in the input.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row = true`, the matrix is row-major, otherwise it is column-major.
 */
template <
    typename Value_ = double, 
    typename Index_ = int, 
    typename StoredValue_ = Value_, 
    typename InputValue_,
    typename InputIndex_
>
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>* incoming, bool row, int threads = 1) {
    auto NR = incoming->nrow();
    auto NC = incoming->ncol();
    std::vector<StoredValue_> buffer(static_cast<size_t>(NR) * static_cast<size_t>(NC));
    convert_to_dense(incoming, row, buffer.data(), threads);
    return std::shared_ptr<Matrix<Value_, Index_> >(new DenseMatrix<Value_, Index_, decltype(buffer)>(NR, NC, std::move(buffer), row));
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam Matrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output dense matrix - row-major (0) or column-major (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `incoming`.
 */
template<typename Value_, typename Index_, typename StoredValue_ = Value_, typename InputValue_, typename InputIndex_>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>* incoming, int order, int threads = 1) {
    if (order < 0) {
        order = static_cast<int>(incoming->prefer_rows()); 
    }
    return convert_to_dense<Value_, Index_, StoredValue_>(incoming, order == 0, threads);
}

/**
 * @cond
 */
// Backwards compatbility.
template<bool row_, typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense(const Matrix<InputValue_, InputIndex_>* incoming, StoredValue_* store, int threads = 1) {
    convert_to_dense(incoming, row_, store, threads);
}

template<bool row_, typename Value_, typename Index_, typename StoredValue_ = Value_, typename InputValue_, typename InputIndex_>
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>* incoming, int threads = 1) {
    return convert_to_dense<Value_, Index_, StoredValue_>(incoming, row_, threads);
}
/**
 * @endcond
 */

}

#endif
