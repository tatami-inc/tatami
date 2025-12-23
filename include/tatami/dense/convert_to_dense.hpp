#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "./DenseMatrix.hpp"

#include "../utils/consecutive_extractor.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/copy.hpp"
#include "../utils/Index_to_container.hpp"

#include <memory>
#include <vector>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

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
void convert_to_dense(const Matrix<InputValue_, InputIndex_>& matrix, const bool row_major, StoredValue_* const store, const ConvertToDenseOptions& options) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const bool pref_rows = matrix.prefer_rows();
    const auto primary = (pref_rows ? NR : NC);
    const auto secondary = (pref_rows ? NC : NR);

    // We assume that 'store' was allocated correctly, in which case the product of 'primary' and 'secondary' is known to fit inside a std::size_t.
    // This saves us from various checks when computing related products (see all the product_unsafe() calls).

    if (row_major == pref_rows) {
        constexpr bool same_type = std::is_same<InputValue_, StoredValue_>::value;
        can_cast_Index_to_container_size<std::vector<InputValue_> >(secondary);

        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<false, InputValue_, InputIndex_>(matrix, pref_rows, start, length);
            auto temp = [&]{
                if constexpr(same_type) {
                    return false;
                } else {
                    return std::vector<InputValue_>(secondary);
                }
            }();

            for (InputIndex_ x = 0; x < length; ++x) {
                const auto store_copy = store + sanisizer::product_unsafe<std::size_t>(secondary, start + x);
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
        std::fill_n(store, sanisizer::product_unsafe<std::size_t>(primary, secondary), 0);

        // We iterate over the input matrix's preferred dimension but split
        // into threads along the non-preferred dimension. This aims to
        // reduce false sharing across threads during writes, as locations
        // for simultaneous writes in the transposed matrix will be
        // separated by around 'secondary * length' elements. 
        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<true, InputValue_, InputIndex_>(matrix, pref_rows, 0, primary, start, length);
            auto vtemp = create_container_of_Index_size<std::vector<InputValue_> >(length);
            auto itemp = create_container_of_Index_size<std::vector<InputIndex_> >(length);

            // Note that we don't use the blocked transposition strategy
            // from the dense case, because the overhead of looping is 
            // worse than the cache misses for sparse data.
            for (InputIndex_ x = 0; x < primary; ++x) {
                const auto range = wrk->fetch(vtemp.data(), itemp.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    store[sanisizer::nd_offset<std::size_t>(x, primary, range.index[i])] = range.value[i];
                }
            }
        }, secondary, options.num_threads);

    } else {
        // Same logic as described for the sparse case; we iterate along the
        // preferred dimension but split into threads along the non-preferred
        // dimension to reduce false sharing.
        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<false, InputValue_, InputIndex_>(matrix, pref_rows, 0, primary, start, length);

            // Performing a blocked transposition to be more
            // cache-friendly. This involves collecting several
            // consecutive primary dimension elements so that we can
            // transpose by blocks along the secondary dimension.
            constexpr InputIndex_ block_size = 16;
            const InputIndex_ alloc = std::min(primary, block_size);
            std::vector<InputValue_> bigbuffer(sanisizer::product_unsafe<typename std::vector<InputValue_>::size_type>(length, alloc));
            std::vector<const InputValue_*> ptrs(alloc); // no need for protection here, we know that alloc <= 16.
            std::vector<InputValue_*> buf_ptrs(alloc);
            for (InputIndex_ i = 0; i < alloc; ++i) {
                buf_ptrs[i] = bigbuffer.data() + sanisizer::product_unsafe<std::size_t>(length, i);
            }

            InputIndex_ prim_i = 0;
            while (prim_i < primary) {
                const InputIndex_ prim_to_process = std::min(static_cast<InputIndex_>(primary - prim_i), block_size);
                for (InputIndex_ c = 0; c < prim_to_process; ++c) {
                    ptrs[c] = wrk->fetch(buf_ptrs[c]);
                }

                InputIndex_ sec_i = 0;
                while (sec_i < length) {
                    const InputIndex_ sec_end = sec_i + std::min(static_cast<InputIndex_>(length - sec_i), block_size);
                    for (InputIndex_ c = 0; c < prim_to_process; ++c) {
                        const auto input = ptrs[c];
                        for (InputIndex_ r = sec_i; r < sec_end; ++r) {
                            store[sanisizer::nd_offset<std::size_t>(c + prim_i, primary, r + start)] = input[r];
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
std::shared_ptr<Matrix<Value_, Index_> > convert_to_dense(const Matrix<InputValue_, InputIndex_>& matrix, const bool row_major, const ConvertToDenseOptions& options) {
    const auto NR = matrix.nrow();
    const auto NC = matrix.ncol();
    const auto buffer_size = sanisizer::product<typename std::vector<StoredValue_>::size_type>(attest_for_Index(NR), attest_for_Index(NC));
    std::vector<StoredValue_> buffer(buffer_size);
    convert_to_dense(matrix, row_major, buffer.data(), options);

    return std::shared_ptr<Matrix<Value_, Index_> >(
        new DenseMatrix<Value_, Index_, I<decltype(buffer)> >(
            sanisizer::cast<Index_>(attest_for_Index(NR)),
            sanisizer::cast<Index_>(attest_for_Index(NC)),
            std::move(buffer),
            row_major
        )
    );
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
    return convert_to_dense<Value_, Index_, StoredValue_>(
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
