#ifndef TATAMI_CONVERT_TO_DENSE_H
#define TATAMI_CONVERT_TO_DENSE_H

#include "./DenseMatrix.hpp"
#include "./transpose.hpp"

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
 * @cond
 */
template <typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense_direct(const Matrix<InputValue_, InputIndex_>& matrix, const bool row, StoredValue_* const store, const ConvertToDenseOptions& options) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const auto primary = (row ? NR : NC);
    const auto secondary = (row ? NC : NR);

    constexpr bool same_type = std::is_same<InputValue_, StoredValue_>::value;
    can_cast_Index_to_container_size<std::vector<InputValue_> >(secondary);

    parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
        auto wrk = consecutive_extractor<false, InputValue_, InputIndex_>(matrix, row, start, length);
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
}

template <typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense_running_from_sparse(const Matrix<InputValue_, InputIndex_>& matrix, const bool row, StoredValue_* const store, const ConvertToDenseOptions& options) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const auto primary = (row ? NR : NC);
    const auto secondary = (row ? NC : NR);

    // Here, our parallelization strategy is to directly put values into the output store for the main thread,
    // while all other threads just store the structural non-zeros in a separate per-thread data structure. 
    // This aims to reduce false sharing at the cost of some extra allocations.
    //
    // Previously, each thread used to process a slice of each secondary dimension element and store it in the output.
    // This avoided false sharing without requiring more memory usage, but was suboptimal in terms of memory locality as adjacent memory was processed by different threads.
    // Each thread also needed to do more calculations to find its slice within each secondary dimension element.
    struct ThreadSpecificHolder {
        ThreadSpecificHolder(const InputIndex_ start, const InputIndex_ length) : 
            start(start),
            values(cast_Index_to_container_size<I<decltype(values)> >(length)),
            indices(cast_Index_to_container_size<I<decltype(indices)> >(length))
        {}
        InputIndex_ start; 
        std::vector<std::vector<InputValue_> > values;
        std::vector<std::vector<InputIndex_> > indices;
    };

    const bool do_parallel = options.num_threads > 1;
    std::optional<std::vector<std::optional<ThreadSpecificHolder> > > all_partial_contents;
    if (do_parallel) {
        // -1, as we don't need to store the results of the main thread.
        all_partial_contents.emplace(sanisizer::cast<I<decltype(all_partial_contents->size())> >(options.num_threads - 1));
    }

    // We assume that 'store' was allocated correctly, in which case the product of 'primary' and 'secondary' is known to fit inside a std::size_t.
    // This saves us from various checks when computing related products. 
    std::fill_n(store, sanisizer::product_unsafe<std::size_t>(primary, secondary), 0);

    const auto num_used = parallelize([&](const int thread, const InputIndex_ start, const InputIndex_ length) -> void {
        auto wrk = consecutive_extractor<true, InputValue_, InputIndex_>(matrix, !row, start, length);
        auto vtemp = create_container_of_Index_size<std::vector<InputValue_> >(primary);
        auto itemp = create_container_of_Index_size<std::vector<InputIndex_> >(primary);

        if (!do_parallel || thread == 0) {
            for (InputIndex_ x = 0; x < length; ++x) {
                const auto range = wrk->fetch(vtemp.data(), itemp.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    store[sanisizer::nd_offset<std::size_t>(start + x, secondary, range.index[i])] = range.value[i];
                }
            }

        } else {
            ThreadSpecificHolder tmp(start, length);
            for (InputIndex_ x = 0; x < length; ++x) {
                const auto range = wrk->fetch(vtemp.data(), itemp.data());
                tmp.values[x] = std::vector<InputValue_>(range.value, range.value + range.number);
                tmp.indices[x] = std::vector<InputIndex_>(range.index, range.index + range.number);
            }
            (*all_partial_contents)[thread - 1] = std::move(tmp);
        }
    }, secondary, options.num_threads);

    // Our reduction step takes the structural non-zeros from other threads and adds them to the output store.
    if (do_parallel) {
        for (int u = 1; u < num_used; ++u) {
            const auto start = (*all_partial_contents)[u - 1]->start;
            const auto& tmp_values = (*all_partial_contents)[u - 1]->values;
            const auto& tmp_indices = (*all_partial_contents)[u - 1]->indices;
            const auto length = tmp_indices.size();
            for (I<decltype(length)> x = 0; x < length ; ++x) {
                const auto& cur_values = tmp_values[x];
                const auto& cur_indices = tmp_indices[x];
                const auto cur_count = cur_indices.size();
                for (I<decltype(cur_count)> i = 0; i < cur_count; ++i) {
                    store[sanisizer::nd_offset<std::size_t>(start + x, secondary, cur_indices[i])] = cur_values[i];
                }
            }
        }
    }
}

template <typename StoredValue_, typename InputValue_, typename InputIndex_>
void convert_to_dense_running_from_dense(const Matrix<InputValue_, InputIndex_>& matrix, const bool row, StoredValue_* const store, const ConvertToDenseOptions& options) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const auto primary = (row ? NR : NC);
    const auto secondary = (row ? NC : NR);

    // Here, our parallelization strategy is to directly put values into the output store for the main thread,
    // while all other threads just store the submatrix values in a separate per-thread data structure. 
    // This aims to reduce false sharing at the cost of some extra allocations.
    //
    // Previously, each thread used to process a slice of each secondary dimension element and store it in the output.
    // This avoided false sharing without requiring more memory usage, but was suboptimal in terms of memory locality as adjacent memory was processed by different threads.
    // Each thread also needed to do more calculations to find its slice within each secondary dimension element.
    const bool do_parallel = options.num_threads > 1;
    struct ThreadSpecificHolder {
        ThreadSpecificHolder(const InputIndex_ start, const InputIndex_ length, const InputIndex_ primary) : 
            start(start),
            length(length),
            values(sanisizer::product<I<decltype(values.size())> >(length, primary)) // still need to check size here, as vector's size_type might be less than size_t.
        {}
        InputIndex_ start, length; 
        std::vector<InputValue_> values;
    };
    std::optional<std::vector<std::optional<ThreadSpecificHolder> > > all_partial_contents;
    if (do_parallel) {
        all_partial_contents.emplace(sanisizer::cast<I<decltype(all_partial_contents->size())> >(options.num_threads - 1));
    }

    // We assume that 'store' was allocated correctly, in which case the product of 'primary' and 'secondary' is known to fit inside a std::size_t.
    // This saves us from various checks when computing related products. 
    std::fill_n(store, sanisizer::product_unsafe<std::size_t>(primary, secondary), 0);

    const auto num_used = parallelize([&](const int thread, const InputIndex_ start, const InputIndex_ length) -> void {
        auto wrk = consecutive_extractor<false, InputValue_, InputIndex_>(matrix, !row, start, length);

        if (!do_parallel || thread == 0) {
            // Performing a blocked transposition to be more cache-friendly.
            // This involves collecting several consecutive primary dimension elements so that we can transpose by blocks along the secondary dimension.
            constexpr InputIndex_ block_size = 16;
            const InputIndex_ alloc = std::min(length, block_size);
            std::vector<InputValue_> bigbuffer(sanisizer::product_unsafe<typename std::vector<InputValue_>::size_type>(primary, alloc));
            std::vector<const InputValue_*> ptrs(alloc); // no need for protection here, we know that alloc <= 16.

            InputIndex_ sec_i = 0;
            while (sec_i < length) {
                const InputIndex_ sec_to_process = std::min(static_cast<InputIndex_>(length - sec_i), block_size);
                for (InputIndex_ x = 0; x < sec_to_process; ++x) {
                    ptrs[x] = wrk->fetch(bigbuffer.data() + sanisizer::product_unsafe<std::size_t>(primary, x));
                }

                InputIndex_ prim_i = 0;
                while (prim_i < primary) {
                    const InputIndex_ prim_end = prim_i + std::min(static_cast<InputIndex_>(primary - prim_i), block_size);
                    for (InputIndex_ x = 0; x < sec_to_process; ++x) {
                        const auto input = ptrs[x];
                        for (InputIndex_ p = prim_i; p < prim_end; ++p) {
                            store[sanisizer::nd_offset<std::size_t>(start + sec_i + x, secondary, p)] = input[p];
                        }
                    }
                    prim_i = prim_end;
                }
                sec_i += sec_to_process;
            }

        } else {
            ThreadSpecificHolder tmp(start, length, primary);
            for (InputIndex_ x = 0; x < length; ++x) {
                const auto buffer = tmp.values.data() + sanisizer::product_unsafe<std::size_t>(x, primary);
                const auto ptr = wrk->fetch(buffer);
                copy_n(ptr, primary, buffer);
            }
            (*all_partial_contents)[thread - 1] = std::move(tmp);
        }
    }, secondary, options.num_threads);

    // Our reduction step performs the transposition to the output array.
    if (do_parallel) {
        for (int u = 1; u < num_used; ++u) {
            const auto tmp_start = (*all_partial_contents)[u - 1]->start;
            const auto tmp_length = (*all_partial_contents)[u - 1]->length;
            const auto& tmp_values = (*all_partial_contents)[u - 1]->values;
            transpose(tmp_values.data(), tmp_length, primary, primary, store + tmp_start, secondary);
        }
    }
}
/**
 * @endcond
 */

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
    if (row_major == matrix.prefer_rows()) {
        convert_to_dense_direct(matrix, row_major, store, options);
    } else if (matrix.is_sparse()) {
        convert_to_dense_running_from_sparse(matrix, row_major, store, options);
    } else {
        convert_to_dense_running_from_dense(matrix, row_major, store, options);
    }
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
