#ifndef TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H
#define TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H

#include <memory>
#include <vector>
#include <cstddef>
#include <optional>

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

            // Special conversion from dense to save ourselves from having to make
            // indices that we aren't really interested in.
            for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                const auto ptr = wrk->fetch(buffer_v.data());
                auto& sv = store_v[p];
                auto& si = store_i[p];

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

    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const InputIndex_ primary = (row ? NR : NC);
    const InputIndex_ secondary = (row ? NC : NR);

    if (!options.two_pass) {
        // In the one-pass strategy, we load everything in a nice format first, then we transpose it in serial.
        // This avoids messy reallocations when trying to expand vectors on an inconsistent dimension.
        auto tmp = retrieve_fragmented_sparse_contents_consistent<StoredValue_, StoredIndex_>(matrix, !row, options);
        auto primary_counts = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
        for (I<decltype(secondary)> s = 0; s < secondary; ++s) {
            const auto& sec_indices = tmp.index[s];
            const auto num = sec_indices.size();
            for (I<decltype(num)> n = 0; n < num; ++n) {
                primary_counts[sec_indices[n]] += 1; // addition must be space, this cannot exceed dimension extents.
            }
        }

        FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
        for (InputIndex_ p = 0; p < primary; ++p) {
            output.index[p].reserve(primary_counts[p]);
            output.value[p].reserve(primary_counts[p]);
        }

        for (I<decltype(secondary)> s = 0; s < secondary; ++s) {
            const auto& sec_values = tmp.value[s];
            const auto& sec_indices = tmp.index[s];
            const auto num = sec_indices.size();
            for (I<decltype(num)> n = 0; n < num; ++n) {
                const auto curp = sec_indices[n];
                output.value[curp].push_back(sec_values[n]);
                output.index[curp].push_back(s);
            }
        }

        return output;
    }

    // In the two-pass strategy, we count the number of non-zeros first, then we fill it up in the second pass.
    std::optional<std::vector<InputIndex_> > nnz_consistent;
    auto nnz_inconsistent = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
    count_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, nnz_inconsistent.data(), nnz_consistent, options.num_threads);

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    for (InputIndex_ p = 0; p < primary; ++p) {
        output.index[p].reserve(nnz_inconsistent[p]);
        output.value[p].reserve(nnz_inconsistent[p]);
    }

    fill_sparse_matrix_inconsistent(
        matrix,
        primary,
        secondary,
        row,
        nnz_consistent,
        /* sparse_main = */ [&](const InputIndex_ s, const SparseRange<InputValue_, InputIndex_>& range) -> void {
            for (InputIndex_ i = 0; i < range.number; ++i) {
                output.value[range.index[i]].push_back(range.value[i]);
                output.index[range.index[i]].push_back(s);
            }
        },
        /* dense_main = */ [&](const InputIndex_ s, const InputValue_* const ptr) -> void {
            for (InputIndex_ p = 0; p < primary; ++p) {
                const auto val = ptr[p]; 
                if (val != 0) {
                    output.value[p].push_back(val);
                    output.index[p].push_back(s);
                }
            }
        },
        /* reduce = */ [&](const InputIndex_ s, const std::vector<InputValue_>& cur_values, const std::vector<InputIndex_>& cur_primary_indices) {
            const auto cur_count = cur_values.size();
            for (I<decltype(cur_count)> i = 0; i < cur_count; ++i) {
                const auto primary = cur_primary_indices[i];
                output.value[primary].push_back(cur_values[i]);
                output.index[primary].push_back(s); 
            }
        },
        options.num_threads
    );

    return output;
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
