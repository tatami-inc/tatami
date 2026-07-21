#ifndef TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H
#define TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H

#include <memory>
#include <vector>
#include <cstddef>
#include <optional>
#include <cassert>

#include "FragmentedSparseMatrix.hpp"
#include "convert_to_sparse_utils.hpp"

#include "../utils/parallelize.hpp"
#include "../utils/copy.hpp"
#include "../utils/consecutive_extractor.hpp"
#include "../utils/Index_to_container.hpp"

/**
 * @file convert_to_fragmented_sparse.hpp
 *
 * @brief Convert a matrix into a fragmented sparse format.
 */

namespace tatami {

/**
 * @brief Fragmented sparse contents.
 *
 * @tparam Value_ Type of value in the matrix.
 * @tparam Index_ Type of row/column index.
 *
 * The "primary" dimension is the one that is used to organize non-zero elements into vectors, while the other dimension is defined as the "secondary" dimension.
 * For example, the rows would be the primary dimension in a fragmented sparse row matrix.
 * (Check out `FragmentedSparseMatrix` for more details.)
 */
template<typename Value_, typename Index_>
struct FragmentedSparseContents {
    /**
     * @cond
     */
    FragmentedSparseContents(Index_ n) :
        value(cast_Index_to_container_size<I<decltype(value)> >(n)),
        index(cast_Index_to_container_size<I<decltype(index)> >(n))
    {}
    /**
     * @endcond
     */

    /**
     * Vector of vectors containing the values of the structural non-zero elements.
     * Each inner vector corresponds to an element of the primary dimension and contains all values for that element.
     */
    std::vector<std::vector<Value_> > value;

    /**
     * Vector of vectors containing the secondary dimension indices of the structural non-zero elements.
     * Each inner vector corresponds to an element of the primary dimension and contains all indices for that element.
     * Each inner vector is of length equal to the corresponding entry of `values` and is guaranteed to be strictly increasing.
     */
    std::vector<std::vector<Index_> > index;
};

/**
 * @brief Options for `retrieve_fragmented_sparse_contents()`.
 */
struct RetrieveFragmentedSparseContentsOptions {
    /**
     * Whether to possibly perform the retrieval in two passes.
     * Setting this to `true` allows the function to perform a preliminary pass through `matrix` to determine the size of each memory allocation.
     * This aims to reduce memory consumption at the cost of some speed.
     */
    bool two_pass = false;

    /**
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

/**
 * @cond
 */
template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_contents_consistent(
    const Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const RetrieveFragmentedSparseContentsOptions& options
) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const InputIndex_ primary = (row ? NR : NC);
    const InputIndex_ secondary = (row ? NC : NR);

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    auto& store_v = output.value;
    auto& store_i = output.index;

    if (matrix.is_sparse()) {
        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<true>(matrix, row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(secondary);
            auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(secondary);

            for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                auto& sv = store_v[p];
                auto& si = store_i[p];
                sv.reserve(range.number);
                si.reserve(range.number);

                // We don't filter out structural non-zeros that have values of zero, for consistency with convert_to_compressed_sparse().
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    sv.push_back(range.value[i]);
                    si.push_back(range.index[i]);
                }
            }
        }, primary, options.num_threads);

    } else {
        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<false>(matrix, row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(secondary);

            for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                const auto ptr = wrk->fetch(buffer_v.data());
                auto& sv = store_v[p];
                auto& si = store_i[p];

                // For dense, we do treat zero values as structural zeros and remove them, otherwise the output wouldn't actually be sparse.
                for (InputIndex_ s = 0; s < secondary; ++s) {
                    const auto val = ptr[s];

                    if (val) {
                        sv.push_back(val);
                        si.push_back(s);
                    }
                }
            }
        }, primary, options.num_threads);
    }

    return output;
}

template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_inconsistent_one_pass(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const int num_threads
) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const InputIndex_ primary = (row ? NR : NC);
    const InputIndex_ secondary = (row ? NC : NR);

    // In the one-pass strategy, we load everything in a nice format first, then we transpose it in serial.
    // This avoids messy reallocations when trying to expand vectors on an inconsistent dimension.
    std::vector<std::vector<InputValue_> > store_v;
    std::vector<std::vector<InputIndex_> > store_i;
    auto original_ranges = extract_sparse_matrix(matrix, store_v, store_i, num_threads);

    auto primary_counts = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
    for (I<decltype(secondary)> s = 0; s < secondary; ++s) {
        const auto& sec_indices = original_ranges[s].index;
        const auto num = original_ranges[s].number;
        for (I<decltype(num)> n = 0; n < num; ++n) {
            primary_counts[sec_indices[n]] += 1; // addition must be safe as this cannot exceed dimension extents.
        }
    }

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    for (InputIndex_ p = 0; p < primary; ++p) {
        output.index[p].reserve(primary_counts[p]);
        output.value[p].reserve(primary_counts[p]);
    }

    for (I<decltype(secondary)> s = 0; s < secondary; ++s) {
        const auto& sec_values = original_ranges[s].value;
        const auto& sec_indices = original_ranges[s].index;
        const auto num = original_ranges[s].number;
        for (I<decltype(num)> n = 0; n < num; ++n) {
            const auto curp = sec_indices[n];
            output.value[curp].push_back(sec_values[n]);
            output.index[curp].push_back(s);
        }
    }

    return output;
}

template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_inconsistent_two_pass(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const int num_threads
) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const InputIndex_ primary = (row ? NR : NC);
    const InputIndex_ secondary = (row ? NC : NR);

    // In the two-pass strategy, we count the number of non-zeros first, then we fill it up in the second pass.
    auto nnz_inconsistent = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
    auto per_thread = count_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, nnz_inconsistent.data(), num_threads);

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    tatami::cast_Index_to_container_size<std::vector<StoredValue_> >(secondary);
    tatami::cast_Index_to_container_size<std::vector<StoredIndex_> >(secondary);
    for (InputIndex_ p = 0; p < primary; ++p) {
        assert(nnz_inconsistent[p] <= secondary);
        output.index[p].resize(nnz_inconsistent[p]);
        output.value[p].resize(nnz_inconsistent[p]);
    }

    const bool is_sparse = matrix.is_sparse();
    if (per_thread.has_value()) {
        // Transforming the per-thread counts into per-thread starting offsets within each vector.
        auto& offsets = per_thread->counts;
        for (InputIndex_ i = 0; i < primary; ++i) {
            InputIndex_ accumulant = 0;
            static_assert(std::is_same<I<decltype(per_thread->counts[0][0])>, InputIndex_>::value); // confirm that the accumulant assignment won't overflow.
            for (auto& pt : offsets) {
                const auto count = pt[i];
                pt[i] = accumulant;
                accumulant += count;
            }
        }

        parallelize([&](const int, const int th_start, const int th_length) -> void {
            for (int t = 0; t < th_length; ++t) {
                auto& offsets = (per_thread->counts)[t + th_start];
                const auto actual_start = (per_thread->starts)[t + th_start];
                const auto actual_length = (per_thread->lengths)[t + th_start];

                // We're going to completely ignore false sharing here, see reasoning in convert_to_compressed_sparse.hpp.
                if (is_sparse) {
                    Options opt;
                    opt.sparse_ordered_index = false;
                    auto wrk = consecutive_extractor<true>(matrix, !row, actual_start, actual_length, opt);
                    auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                    auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
                    for (InputIndex_ x = 0; x < actual_length; ++x) {
                        const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                        for (InputIndex_ i = 0; i < range.number; ++i) {
                            const auto prim = range.index[i];
                            auto& pos = offsets[prim];
                            output.value[prim][pos] = range.value[i];
                            output.index[prim][pos] = x + actual_start;
                            ++pos;
                        }
                    }

                } else {
                    auto wrk = consecutive_extractor<false>(matrix, !row, actual_start, actual_length);
                    auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                    for (InputIndex_ x = 0; x < actual_length; ++x) {
                        const auto ptr = wrk->fetch(buffer_v.data());
                        for (InputIndex_ p = 0; p < primary; ++p) {
                            const auto val = ptr[p]; 
                            if (val != 0) {
                                auto& pos = offsets[p];
                                output.value[p][pos] = val;
                                output.index[p][pos] = x + actual_start; 
                                ++pos;
                            }
                        }
                    }
                }
            }
        }, per_thread->counts.size(), per_thread->counts.size());

    } else {
        auto offsets = tatami::create_container_of_Index_size<std::vector<InputIndex_> >(primary);

        if (is_sparse){ 
            Options opt;
            opt.sparse_ordered_index = false;
            auto wrk = consecutive_extractor<true>(matrix, !row, static_cast<InputIndex_>(0), secondary, opt);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
            auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
            for (InputIndex_ s = 0; s < secondary; ++s) {
                const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    const auto prim = range.index[i];
                    auto& pos = offsets[prim];
                    output.value[prim][pos] = range.value[i];
                    output.index[prim][pos] = s;
                    ++pos;
                }
            }

        } else {
            auto wrk = consecutive_extractor<false>(matrix, !row, static_cast<InputIndex_>(0), secondary);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
            for (InputIndex_ s = 0; s < secondary; ++s) {
                const auto ptr = wrk->fetch(buffer_v.data());
                for (InputIndex_ p = 0; p < primary; ++p) {
                    const auto val = ptr[p]; 
                    if (val != 0) {
                        auto& pos = offsets[p];
                        output.value[p][pos] = val;
                        output.index[p][pos] = s; 
                        ++pos;
                    }
                }
            }
        }
    }

    return output;
}
/**
 * @endcond
 */

/**
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix Pointer to a `tatami::Matrix`. 
 * @param row Whether to retrieve the contents of `matrix` by row, i.e., the output is a fragmented sparse row matrix.
 * @param options Further options.
 *
 * @return Contents of the sparse matrix in fragmented form, see `FragmentedSparseContents`.
 */
template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_contents(
    const Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const RetrieveFragmentedSparseContentsOptions& options
) {
    if (row == matrix.prefer_rows()) {
        return retrieve_fragmented_sparse_contents_consistent<StoredValue_, StoredIndex_>(matrix, row, options);
    }

    if (!options.two_pass) {
        return retrieve_fragmented_sparse_inconsistent_one_pass<StoredValue_, StoredIndex_>(matrix, row, options.num_threads);
    }

    return retrieve_fragmented_sparse_inconsistent_two_pass<StoredValue_, StoredIndex_>(matrix, row, options.num_threads);
}

/**
 * @brief Options for `convert_to_fragmented_sparse()`.
 */
struct ConvertToFragmentedSparseOptions {
    /**
     * Whether to possibly perform the conversion in two passes.
     * Setting this to `true` allows the function to perform a preliminary pass through `matrix` to determine the size of each memory allocation.
     * This aims to reduce memory consumption at the cost of some speed.
     */
    bool two_pass = false;

    /**
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

/**
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix A `tatami::Matrix`. 
 * @param row Whether to return a fragmented sparse row matrix.
 * @param options Further options.
 *
 * @return A pointer to a new `tatami::FragmentedSparseMatrix`, with the same dimensions and type as the matrix referenced by `matrix`.
 * If `row = true`, the matrix is in fragmented sparse row format, otherwise it is fragmented sparse column.
 */
template<
    typename Value_,
    typename Index_,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_fragmented_sparse(
    const Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const ConvertToFragmentedSparseOptions& options)
{
    auto frag = retrieve_fragmented_sparse_contents<StoredValue_, StoredIndex_>(
        matrix,
        row,
        [&]{
            RetrieveFragmentedSparseContentsOptions ropt;
            ropt.two_pass = options.two_pass;
            ropt.num_threads = options.num_threads;
            return ropt;
        }()
    );
    return std::shared_ptr<Matrix<Value_, Index_> >(
        new FragmentedSparseMatrix<
            Value_, 
            Index_,
            std::vector<std::vector<StoredValue_> >,
            std::vector<std::vector<StoredIndex_> >
        >(
            matrix.nrow(), 
            matrix.ncol(), 
            std::move(frag.value), 
            std::move(frag.index),
            row, 
            []{
                FragmentedSparseMatrixOptions fopt;
                fopt.check = false; // no need for checks, as we guarantee correctness.
                return fopt;
            }()
        )
    );
}

/**
 * @cond
 */
// Backwards compatbility.
template<typename Value_, typename Index_, typename StoredValue_ = Value_, typename StoredIndex_ = Index_, typename InputValue_, typename InputIndex_>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_fragmented_sparse(const Matrix<InputValue_, InputIndex_>* matrix, bool row, int threads = 1) {
    return convert_to_fragmented_sparse<Value_, Index_, StoredValue_, StoredIndex_>(
        *matrix,
        row,
        [&]{
            ConvertToFragmentedSparseOptions opt;
            opt.num_threads = threads;
            return opt;
        }()
    );
}

template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, bool row, int threads = 1) {
    return retrieve_fragmented_sparse_contents<StoredValue_, StoredIndex_>(
        *matrix,
        row,
        [&]{
            RetrieveFragmentedSparseContentsOptions opt;
            opt.num_threads = threads;
            return opt;
        }()
    );
}

template <bool row_, typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, int threads = 1) {
    return retrieve_fragmented_sparse_contents<StoredValue_, StoredIndex_>(matrix, row_, threads);
}

template <bool row_, typename Value_, typename Index_, typename StoredValue_ = Value_, typename StoredIndex_ = Index_, typename InputValue_, typename InputIndex_>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_fragmented_sparse(const Matrix<InputValue_, InputIndex_>* matrix, int threads = 1) {
    return convert_to_fragmented_sparse<Value_, Index_, StoredValue_, StoredIndex_>(matrix, row_, threads);
}
/**
 * @endcond
 */

}

#endif
