#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "./DenseMatrix.hpp"
#include "./transpose.hpp"

#include "../utils/consecutive_extractor.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/copy.hpp"

#include <memory>
#include <vector>
#include <cstddef>

/**
 * @file convert_to_dense.hpp
 *
 * @brief Convert a matrix into a dense format.
 */

namespace tatami {

/**
 * @brief Options for `convert_to_dense()`.
 */
struct ConvertToDenseOptions {
    /**
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

/**
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam InputValue_ Type of data values in the input.
 * @tparam InputIndex_ Integer type for the indices in the input.
 *
 * @param matrix A `tatami::Matrix`.
 * @param row_major Whether to store the output as a row-major matrix.
 * @param[out] store Pointer to an array of length equal to the product of the dimensions of `matrix`.
 * On output, this is filled with values from `matrix` in row- or column-major format depending on `row_major`.
 * @param options Further options.
 */
template <typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense(const Matrix<InputValue_, InputIndex_>& matrix, bool row_major, StoredValue_* store, const ConvertToDenseOptions& options) {
    InputIndex_ NR = matrix.nrow();
    InputIndex_ NC = matrix.ncol();
    bool pref_rows = matrix.prefer_rows();
    std::size_t primary = (pref_rows ? NR : NC);
    std::size_t secondary = (pref_rows ? NC : NR);

    if (row_major == pref_rows) {
        constexpr bool same_type = std::is_same<InputValue_, StoredValue_>::value;
        parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
            std::vector<InputValue_> temp(same_type ? 0 : secondary);
            auto wrk = consecutive_extractor<false>(matrix, pref_rows, start, length);

            for (InputIndex_ x = 0; x < length; ++x) {
                auto store_copy = store + static_cast<std::size_t>(start + x) * secondary; // cast to avoid overflow.
                if constexpr(same_type) {
                    auto ptr = wrk->fetch(store_copy);
                    copy_n(ptr, secondary, store_copy);
                } else {
                    auto ptr = wrk->fetch(temp.data());
                    std::copy_n(ptr, secondary, store_copy);
                }
            }
        }, primary, options.num_threads);

    } else if (matrix.is_sparse()) {
        std::fill_n(store, primary * secondary, 0); // already cast to std::size_t to avoid overflow.

        // We iterate over the input matrix's preferred dimension but split
        // into threads along the non-preferred dimension. This aims to
        // reduce false sharing across threads during writes, as locations
        // for simultaneous writes in the transposed matrix will be
        // separated by around 'secondary * length' elements. 
        parallelize([&](int, std::size_t start, std::size_t length) -> void {
            auto wrk = consecutive_extractor<true, InputValue_, InputIndex_>(matrix, pref_rows, 0, primary, start, length);
            std::vector<InputValue_> vtemp(length);
            std::vector<InputIndex_> itemp(length);

            // Note that we don't use the blocked transposition strategy
            // from the dense case, because the overhead of looping is 
            // worse than the cache misses for sparse data.
            for (std::size_t x = 0; x < primary; ++x) {
                auto range = wrk->fetch(vtemp.data(), itemp.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    store[static_cast<std::size_t>(range.index[i]) * primary + x] = range.value[i]; // cast to std::size_t to avoid overflow.
                }
            }
        }, secondary, options.num_threads);

    } else {
        // Same logic as described for the sparse case; we iterate along the
        // preferred dimension but split into threads along the non-preferred
        // dimension to reduce false sharing.
        parallelize([&](int, std::size_t start, std::size_t length) -> void {
            auto wrk = consecutive_extractor<false, InputValue_, InputIndex_>(matrix, pref_rows, 0, primary, start, length);

            // Performing a blocked transposition to be more
            // cache-friendly. This involves collecting several
            // consecutive primary dimension elements so that we can
            // transpose by blocks along the secondary dimension.
            constexpr std::size_t block_size = 16;
            const std::size_t alloc = std::min(primary, block_size);
            std::vector<InputValue_> bigbuffer(length * alloc); // already size_ts, to avoid overflow.
            std::vector<const InputValue_*> ptrs(alloc);
            std::vector<InputValue_*> buf_ptrs(alloc);
            for (std::size_t i = 0; i < alloc; ++i) {
                buf_ptrs[i] = bigbuffer.data() + i * length; // already all size_t's, to avoid overflow.
            }

            std::size_t prim_i = 0;
            while (prim_i < primary) {
                std::size_t prim_to_process = std::min(static_cast<std::size_t>(primary - prim_i), block_size);
                for (std::size_t c = 0; c < prim_to_process; ++c) {
                    ptrs[c] = wrk->fetch(buf_ptrs[c]);
                }

                std::size_t sec_i = 0;
                while (sec_i < length) {
                    std::size_t sec_end = sec_i + std::min(static_cast<std::size_t>(length - sec_i), block_size);
                    for (std::size_t c = 0; c < prim_to_process; ++c) {
                        auto input = ptrs[c];
                        std::size_t offset = start * primary + (c + prim_i); // already all std::size_t's, to avoid overflow.
                        for (std::size_t r = sec_i; r < sec_end; ++r) {
                            store[r * primary + offset] = input[r];
                        }
                    }

                    sec_i = sec_end;
                }
                prim_i += prim_to_process;
            }
        }, secondary, options.num_threads);
    }

    return;
}

/**
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam InputValue_ Type of data values in the input.
 * @tparam InputIndex_ Integer type for the indices in the input.
 *
 * @param matrix A `tatami::Matrix`.
 * @param row_major Whether to return a row-major matrix.
 * @param options Further options.
 *
 * @return A pointer to a new `tatami::DenseMatrix` with the same dimensions and type as the matrix referenced by `matrix`.
 * If `row_major = true`, the matrix is row-major, otherwise it is column-major.
 */
template <
    typename Value_,
    typename Index_,
    typename StoredValue_ = Value_, 
    typename InputValue_,
    typename InputIndex_
>
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>& matrix, bool row_major, const ConvertToDenseOptions& options) {
    auto NR = matrix.nrow();
    auto NC = matrix.ncol();
    std::vector<StoredValue_> buffer(static_cast<std::size_t>(NR) * static_cast<std::size_t>(NC));
    convert_to_dense(matrix, row_major, buffer.data(), options);
    return std::shared_ptr<Matrix<Value_, Index_> >(new DenseMatrix<Value_, Index_, decltype(buffer)>(NR, NC, std::move(buffer), row_major));
}

/**
 * @cond
 */
// Backwards compatbility.
template <typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense(const Matrix<InputValue_, InputIndex_>* matrix, bool row_major, StoredValue_* store, int threads = 1) {
    convert_to_dense(
        *matrix,
        row_major,
        store,
        [&]{
            ConvertToDenseOptions options;
            options.num_threads = threads;
            return options;
        }()
    );
}

template <typename Value_ = double, typename Index_ = int, typename StoredValue_ = Value_, typename InputValue_, typename InputIndex_>
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>* matrix, bool row_major, int threads = 1) {
    ConvertToDenseOptions options;
    options.num_threads = threads;
    return convert_to_dense<Value_, Index_>(
        *matrix,
        row_major,
        [&]{
            ConvertToDenseOptions options;
            options.num_threads = threads;
            return options;
        }()
    );
}

template<bool row_, typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense(const Matrix<InputValue_, InputIndex_>* matrix, StoredValue_* store, int threads = 1) {
    convert_to_dense(matrix, row_, store, threads);
}

template<bool row_, typename Value_, typename Index_, typename StoredValue_ = Value_, typename InputValue_, typename InputIndex_>
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>* matrix, int threads = 1) {
    return convert_to_dense<Value_, Index_, StoredValue_>(matrix, row_, threads);
}
/**
 * @endcond
 */

}

#endif
