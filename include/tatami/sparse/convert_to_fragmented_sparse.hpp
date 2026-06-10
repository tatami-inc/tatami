#ifndef TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H
#define TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H

#include <memory>
#include <vector>
#include <cstddef>
#include <optional>

#include "FragmentedSparseMatrix.hpp"
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
     * Setting this to `true` allows the function to perform a preliminary pass through `matrix` to determine the memory allocation.
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

template<typename Index_>
struct CountFragmentedSparseNonZerosInconsistentResult {
    CountFragmentedSparseNonZerosInconsistentResult(const Index_ primary, const Index_ secondary, const bool is_sparse) :
        inconsistent(tatami::cast_Index_to_container_size<std::vector<Index_> >(primary))
    {
        if (!is_sparse) {
            consistent.emplace(cast_Index_to_container_size<std::vector<Index_> >(secondary));
        }
    }

    std::vector<Index_> inconsistent;
    std::optional<std::vector<Index_> > consistent;
};

template<typename Value_, typename Index_>
CountFragmentedSparseNonZerosInconsistentResult<Index_> count_fragmented_sparse_non_zeros_inconsistent(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Index_ primary,
    const Index_ secondary,
    const bool row,
    const int threads
) {
    const bool do_parallel = threads > 1;
    std::optional<std::vector<std::optional<std::vector<Index_> > > > all_partial_counts;
    if (do_parallel) {
        all_partial_counts.emplace(sanisizer::cast<I<decltype(all_partial_counts->size())> >(threads - 1));
    }

    const bool is_sparse = matrix.is_sparse();
    CountFragmentedSparseNonZerosInconsistentResult<Index_> output(primary, secondary, is_sparse);
    auto& nnz_inconsistent = output.inconsistent;
    auto& nnz_consistent = output.consistent;

    const int num_used = parallelize([&](const int thread, const Index_ start, const Index_ length) -> void {
        // To minimize false sharing, we allocate each buffer as a per-thread vector before moving it into the nnz_workers for serial use.
        // We skip the allocation for the first thread as this is allowed to use the (presumably zeroed) nnz array directly.
        Index_* cur_counts;
        std::optional<std::vector<Index_> > count_holder;
        if (!do_parallel) {
            cur_counts = nnz_inconsistent.data();
        } else {
            if (thread == 0) {
                cur_counts = nnz_inconsistent.data();
            } else {
                count_holder.emplace(cast_Index_to_container_size<std::vector<Index_> >(primary));
                cur_counts = count_holder->data();
            }
        }

        if (is_sparse) {
            Options opt;
            opt.sparse_extract_value = false;
            opt.sparse_ordered_index = false;
            auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
            auto buffer_i = create_container_of_Index_size<std::vector<Index_> >(primary);

            for (Index_ x = 0; x < length; ++x) {
                const auto range = wrk->fetch(NULL, buffer_i.data());
                for (Index_ i = 0; i < range.number; ++i) {
                    ++cur_counts[range.index[i]];
                }
            }

        } else {
            auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<Value_> >(primary);

            for (Index_ x = 0; x < length; ++x) {
                const auto ptr = wrk->fetch(buffer_v.data());
                Index_ count = 0;
                for (Index_ p = 0; p < primary; ++p) {
                    const bool is_nz = (ptr[p] != 0);
                    cur_counts[p] += is_nz;
                    count += is_nz;
                }
                (*nnz_consistent)[start + x] = count;
            }
        }

        if (do_parallel) {
            if (thread > 0) {
                (*all_partial_counts)[thread - 1] = std::move(count_holder);
            }
        }
    }, secondary, threads);

    if (do_parallel) {
        for (int t = 1; t < num_used; ++t) {
            const auto& y = *((*all_partial_counts)[t - 1]);
            for (Index_ p = 0; p < primary; ++p) {
                nnz_inconsistent[p] += y[p];
            }
        }
    }

    return output;
}

template<typename InputValue_, typename InputIndex_, typename StoredValue_, typename StoredIndex_>
void fill_fragmented_sparse_matrix_inconsistent(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const InputIndex_ primary,
    const InputIndex_ secondary,
    const bool row,
    FragmentedSparseContents<StoredValue_, StoredIndex_>& output,
    const std::optional<std::vector<InputIndex_> >& nnz_consistent, 
    const int threads
) {
    // Unlike the corresponding *compressed_sparse() function, we don't try to allocate a big block per thread and fill it up.
    // We don't have a concept of a Pointer_ type here and we can't be sure that the number of non-zeros wouldn't overflow whatever type we picked internally
    // (or indeed, the maximum size of a std::vector, for that matter).
    // So, we just allocate one vector per secondary dimension element, even though that does have some higher overhead. 
    const bool do_parallel = threads > 0;
    InputIndex_ last_secondary_added = 0;
    std::optional<std::vector<std::vector<StoredValue_> > > all_partial_values;
    std::optional<std::vector<std::vector<StoredIndex_> > > all_partial_primary_indices;
    if (do_parallel) {
        all_partial_values.emplace(sanisizer::cast<I<decltype(all_partial_values->size())> >(secondary));
        all_partial_primary_indices.emplace(sanisizer::cast<I<decltype(all_partial_primary_indices->size())> >(secondary));
    }

    const bool is_sparse = matrix.is_sparse();
    parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
        if (start == 0) {
            last_secondary_added = start + length;

            if (is_sparse){ 
                Options opt;
                opt.sparse_ordered_index = false;
                auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i) {
                        output.value[range.index[i]].push_back(range.value[i]);
                        output.index[range.index[i]].push_back(x); // start == 0, so no need to add it to 'x'.
                    }
                }

            } else {
                auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = 0; p < primary; ++p) {
                        const auto val = ptr[p]; 
                        if (val != 0) {
                            output.value[p].push_back(val);
                            output.index[p].push_back(x); // start == 0, so no need to add it to 'x'.
                        }
                    }
                }
            }

        } else {
            if (is_sparse){ 
                Options opt;
                opt.sparse_ordered_index = false;
                auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    (*all_partial_values)[start + x] = std::vector<StoredValue_>(range.value, range.value + range.number);
                    (*all_partial_primary_indices)[start + x] = std::vector<StoredIndex_>(range.index, range.index + range.number);
                }

            } else {
                auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);

                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto nnz = (*nnz_consistent)[start + x]; 
                    std::vector<StoredValue_> out_v;
                    std::vector<StoredIndex_> out_i;
                    out_v.reserve(nnz);
                    out_i.reserve(nnz);

                    const auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = 0; p < primary; ++p) {
                        const auto val = ptr[p]; 
                        if (val != 0) {
                            out_v.push_back(val);
                            out_i.push_back(p);
                        }
                    }

                    (*all_partial_values)[start + x] = std::move(out_v);
                    (*all_partial_primary_indices)[start + x] = std::move(out_i);
                }
            }
        }
    }, secondary, threads);

    if (do_parallel) {
        for (InputIndex_ s = last_secondary_added; s < secondary; ++s) {
            const auto& cur_primary_indices = (*all_partial_primary_indices)[s];
            const auto& cur_values = (*all_partial_values)[s];
            const auto cur_count = cur_values.size();
            for (I<decltype(cur_count)> i = 0; i < cur_count; ++i) {
                const auto primary = cur_primary_indices[i];
                output.value[primary].push_back(cur_values[i]);
                output.index[primary].push_back(s); 
            }
        }
    }
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
    auto nnz_counts = count_fragmented_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, options.num_threads);
    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    for (InputIndex_ p = 0; p < primary; ++p) {
        output.index[p].reserve(nnz_counts.inconsistent[p]);
        output.value[p].reserve(nnz_counts.inconsistent[p]);
    }
    fill_fragmented_sparse_matrix_inconsistent(matrix, primary, secondary, row, output, nnz_counts.consistent, options.num_threads);
    return output;
}

/**
 * @brief Options for `convert_to_fragmented_sparse()`.
 */
struct ConvertToFragmentedSparseOptions {
    /**
     * Whether to possibly perform the retrieval in two passes.
     * Setting this to `true` allows the function to perform a preliminary pass through `matrix` to determine the memory allocation.
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
