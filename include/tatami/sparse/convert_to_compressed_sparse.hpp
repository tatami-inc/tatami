#ifndef TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H
#define TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H

#include <memory>
#include <vector>
#include <cstddef>
#include <optional>

#include "CompressedSparseMatrix.hpp"
#include "convert_to_fragmented_sparse.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/consecutive_extractor.hpp"
#include "../utils/Index_to_container.hpp"
#include "../utils/copy.hpp"

/**
 * @file convert_to_compressed_sparse.hpp
 *
 * @brief Convert a matrix into a compressed sparse format.
 */

namespace tatami {

/**
 * @cond
 */
template<typename Value_, typename Index_, typename Count_>
void count_compressed_sparse_non_zeros_consistent(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Index_ primary,
    const Index_ secondary,
    const bool row,
    Count_* const output,
    const int threads
) {
    sanisizer::cast<Count_>(secondary); // confirm that the counts don't overflow the Count_.

    if (matrix.is_sparse()) {
        Options opt;
        opt.sparse_extract_value = false;
        opt.sparse_extract_index = false;
        opt.sparse_ordered_index = false;

        parallelize([&](const int, const Index_ start, const Index_ length) -> void {
            auto wrk = consecutive_extractor<true>(matrix, row, start, length, opt);
            for (Index_ x = 0; x < length; ++x) {
                const auto range = wrk->fetch(NULL, NULL);
                output[start + x] = range.number;
            }
        }, primary, threads);

    } else {
        parallelize([&](const int, const Index_ start, const Index_ length) -> void {
            auto buffer_v = create_container_of_Index_size<std::vector<Value_> >(secondary);
            auto wrk = consecutive_extractor<false>(matrix, row, start, length);
            for (Index_ p = start, pe = start + length; p < pe; ++p) {
                const auto ptr = wrk->fetch(buffer_v.data());
                Count_ count = 0;
                for (Index_ s = 0; s < secondary; ++s) {
                    count += (ptr[s] != 0);
                }
                output[p] = count;
            }
        }, primary, threads);
    }
}

template<typename Value_, typename Index_, typename Count_>
void count_compressed_sparse_non_zeros_inconsistent(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Index_ primary,
    const Index_ secondary,
    const bool row,
    Count_* const nnz_inconsistent,
    Index_* const nnz_consistent,
    const int threads
) {
    sanisizer::cast<Count_>(secondary); // confirm that the counts don't overflow the Count_.
    const bool is_sparse = matrix.is_sparse();

    const bool do_parallel = threads > 1;
    std::optional<std::vector<std::optional<std::vector<Count_> > > > all_partial_counts;
    if (do_parallel) {
        all_partial_counts.emplace(sanisizer::cast<I<decltype(all_partial_counts->size())> >(threads - 1));
    }

    std::fill_n(nnz_inconsistent, primary, 0);

    const int num_used = parallelize([&](const int thread, const Index_ start, const Index_ length) -> void {
        // To minimize false sharing, we allocate each buffer as a per-thread vector before moving it into the nnz_workers for serial use.
        // We skip the allocation for the first thread as this is allowed to use the (presumably zeroed) nnz array directly.
        Count_* cur_counts;
        std::optional<std::vector<Count_> > count_holder;
        if (!do_parallel) {
            cur_counts = nnz_inconsistent;
        } else {
            if (thread == 0) {
                cur_counts = nnz_inconsistent;
            } else {
                count_holder.emplace(cast_Index_to_container_size<std::vector<Count_> >(primary));
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
                if (nnz_consistent != NULL) {
                    nnz_consistent[start + x] = range.number;
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
                if (nnz_consistent != NULL) {
                    nnz_consistent[start + x] = count;
                }
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
}

// For back-compatiblity only, this functionality should probably not have been exported.
// It's not even entirely correct as we only count structural non-zeros in the sparse case,
// so it's hard to think of a case where someone would want to use this.
struct CountCompressedSparseNonZerosOptions {
    int num_threads = 1;
};

// For back-compatiblity only, see above.
template<typename Value_, typename Index_, typename Count_>
void count_compressed_sparse_non_zeros(
    const tatami::Matrix<Value_, Index_>& matrix,
    const bool row,
    Count_* const output,
    const CountCompressedSparseNonZerosOptions& options
) {
    const Index_ NR = matrix.nrow();
    const Index_ NC = matrix.ncol();
    const Index_ primary = (row ? NR : NC);
    const Index_ secondary = (row ? NC : NR);

    if (row == matrix.prefer_rows()) {
        count_compressed_sparse_non_zeros_consistent(matrix, primary, secondary, row, output, options.num_threads);
    } else {
        count_compressed_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, output, static_cast<Index_*>(NULL), options.num_threads);
    }
}

template<typename InputValue_, typename InputIndex_, typename Pointer_, typename StoredValue_, typename StoredIndex_>
void fill_compressed_sparse_matrix_consistent(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const InputIndex_ primary,
    const InputIndex_ secondary,
    const bool row,
    const Pointer_* const pointers,
    StoredValue_* const output_value,
    StoredIndex_* const output_index,
    const int threads
) {
    if (matrix.is_sparse()) {
        Options opt;
        opt.sparse_ordered_index = false;

        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<true>(matrix, row, start, length, opt);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(secondary);
            auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(secondary);

            for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                // Resist the urge to `fetch()` straight into 'output_v'
                // and 'output_i', as implementations may assume that they
                // have the entire 'length' length to play with, and the
                // output vectors only have whatever is allocated from the
                // first pass (which might be nothing for an all-zero matrix).
                const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                const auto offset = pointers[p];
                std::copy_n(range.value, range.number, output_value + offset);
                std::copy_n(range.index, range.number, output_index + offset);
            }
        }, primary, threads);

    } else {
        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(secondary);
            auto wrk = consecutive_extractor<false>(matrix, row, start, length);

            for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                const auto ptr = wrk->fetch(buffer_v.data());
                auto offset = pointers[p];
                for (InputIndex_ s = 0; s < secondary; ++s) {
                    const auto val = ptr[s];
                    if (val != 0) {
                        output_value[offset] = val;
                        output_index[offset] = s;
                        ++offset;
                    }
                }
            }
        }, primary, threads);
    }
}

// Here, the general strategy is to store all non-zero elements in thread-specific containers, and then stitch them together into the output buffers in the serial section.
// This avoids false sharing at the cost of almost doubling the memory usage as the number of threads increases - I think this is mostly acceptable.
// (There is no penalty for single-threaded use as we just use the output buffers directly.) 
//
// Contrast this to the previous approach, where we take a different slice of each primary dimension element in each thread.
// This tends to be a suboptimal extraction strategy as each thread would need to do some indexing to find the slice.
// Depending on the representation, the work done in each thread might be significant, e.g., a full read of the indices from file in a HDF5-backed sparse matrix.
template<typename InputValue_, typename InputIndex_, typename Pointer_, typename StoredValue_, typename StoredIndex_>
void fill_compressed_sparse_matrix_inconsistent(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const InputIndex_ primary,
    const InputIndex_ secondary,
    const bool row,
    std::vector<Pointer_> offsets, // indptrs of the entire matrix. We'll need to modify this, so we take a copy.
    const std::vector<InputIndex_>& secondary_counts, // count of the number of non-zero elements in each secondary dimension element, for buffer allocations.
    StoredValue_* const output_value,
    StoredIndex_* const output_index,
    const int threads
) {
    const bool do_parallel = threads > 0;
    InputIndex_ last_secondary_added = 0;
    std::optional<std::vector<int> > secondary_thread_assignment;
    std::optional<std::vector<std::optional<std::vector<StoredValue_> > > > all_partial_values;
    std::optional<std::vector<std::optional<std::vector<StoredIndex_> > > > all_partial_primary_indices;
    if (do_parallel) {
        // Don't subtract 1 for the threads, thread 0 might not be processing start == 0.
        // We aren't guaranteed that thread numbers are positive when storing the intermediate results. 
        secondary_thread_assignment.emplace(cast_Index_to_container_size<std::vector<int> >(secondary));
        all_partial_values.emplace(sanisizer::cast<I<decltype(all_partial_values->size())> >(threads));
        all_partial_primary_indices.emplace(sanisizer::cast<I<decltype(all_partial_primary_indices->size())> >(threads));
    }

    const bool is_sparse = matrix.is_sparse();
    parallelize([&](const int thread, const InputIndex_ start, const InputIndex_ length) -> void {
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
                        auto& pos = offsets[range.index[i]];
                        output_value[pos] = range.value[i];
                        output_index[pos] = x; // start == 0, so no need to add it.
                        ++pos;
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
                            auto& pos = offsets[p];
                            output_value[pos] = val;
                            output_index[pos] = x; // start == 0, so no need to add it.
                            ++pos;
                        }
                    }
                }
            }

        } else {
            Pointer_ cur_nnz = 0;
            for (InputIndex_ s = start, end = start + length; s < end; ++s) {
                // This summation is known to be safe as Pointer_ can hold the total number of non-zeros in the entire matrix.
                cur_nnz += secondary_counts[s];
            }
            auto cur_values = sanisizer::create<std::vector<StoredValue_> >(cur_nnz);
            auto cur_primary_indices = sanisizer::create<std::vector<StoredIndex_> >(cur_nnz);

            if (is_sparse){ 
                Options opt;
                opt.sparse_ordered_index = false;
                auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
                Pointer_ position = 0;
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    std::copy_n(range.index, range.number, cur_primary_indices.begin() + position);
                    std::copy_n(range.value, range.number, cur_values.begin() + position);
                    position += range.number;
                    (*secondary_thread_assignment)[start + x] = thread;
                }

            } else {
                auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                Pointer_ position = 0;
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = 0; p < primary; ++p) {
                        const auto val = ptr[p]; 
                        if (val != 0) {
                            cur_primary_indices[position] = p;
                            cur_values[position] = val;
                            ++position;
                        }
                    }
                    (*secondary_thread_assignment)[start + x] = thread;
                }
            }

            (*all_partial_primary_indices)[thread] = std::move(cur_primary_indices);
            (*all_partial_values)[thread] = std::move(cur_values);
        }
    }, secondary, threads);

    if (do_parallel) {
        auto positions = sanisizer::create<std::vector<Pointer_> >(threads);
        for (InputIndex_ s = last_secondary_added; s < secondary; ++s) {
            const auto t = (*secondary_thread_assignment)[s];
            const auto& cur_primary_indices = *((*all_partial_primary_indices)[t]);
            const auto& cur_values = *((*all_partial_values)[t]);
            auto& src_pos = positions[t];

            const auto cur_count = secondary_counts[s];
            for (InputIndex_ i = 0; i < cur_count; ++i) {
                auto& dest_pos = offsets[cur_primary_indices[i + src_pos]];
                output_value[dest_pos] = cur_values[i + src_pos];
                output_index[dest_pos] = s; 
                ++dest_pos;
            }

            src_pos += cur_count;;
        }
    }
}

// For back-compatiblity only, this functionality should probably not have been exported.
// This is only useful in the context of retrieve_compressed_sparse_contents, so why would someone use this when they could just use retrieve?
struct FillCompressedSparseContentsOptions {
    int num_threads = 1;
};

// For back-compatiblity only, see above.
template<typename InputValue_, typename InputIndex_, typename Pointer_, typename StoredValue_, typename StoredIndex_>
void fill_compressed_sparse_contents(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const Pointer_* const pointers,
    StoredValue_* const output_value,
    StoredIndex_* const output_index,
    const FillCompressedSparseContentsOptions& options
) {
    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const InputIndex_ primary = (row ? NR : NC);
    const InputIndex_ secondary = (row ? NC : NR);

    if (row == matrix.prefer_rows()) {
        fill_compressed_sparse_matrix_consistent(matrix, primary, secondary, row, pointers, output_value, output_index, options.num_threads);
    } else {
        // Force it to be single-threaded as we don't have the per-worker pointers to parallelize effectively.
        fill_compressed_sparse_matrix_inconsistent(
            matrix,
            primary,
            secondary,
            row,
            std::vector<Pointer_>(pointers, pointers + sanisizer::sum<std::size_t>(primary, 1)), // making a copy because this is modified on input.
            {},
            output_value,
            output_index,
            1
        );
    }
}
/**
 * @endcond
 */

/**
 * @brief Compressed sparse contents.
 *
 * @tparam Value_ Type of value in the matrix.
 * @tparam Index_ Type of row/column index.
 * @tparam Pointer_ Integer type for the row/column pointers.
 *
 * The "primary" dimension is the one that is used to create the pointers for the compressed sparse format, while the other dimension is defined as the "secondary" dimension.
 * For example, the rows would be the primary dimension in a compressed sparse row matrix.
 */
template<typename Value_, typename Index_, typename Pointer_>
struct CompressedSparseContents {
    /**
     * Vector containing values of the structural non-zero elements in a compressed sparse format.
     */
    std::vector<Value_> value;

    /**
     * Vector containing the secondary dimension indices of the structural non-zero elements in a compressed sparse format.
     */
    std::vector<Index_> index;

    /**
     * Vector containing the pointers for each primary dimension element in a compressed sparse format.
     */
    std::vector<Pointer_> pointers;
};

/**
 * @brief Options for `retrieve_compressed_sparse_contents()`.
 */
struct RetrieveCompressedSparseContentsOptions {
    /**
     * Whether to perform the retrieval in two passes.
     * This requires another pass through `matrix` but is more memory-efficient.
     */
    bool two_pass = false;

    /**
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

/**
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the row/column indices in the output. 
 * @tparam StoredPointer_ Integer type for the row/column pointers in the output.
 * This should be large enough to hold the number of non-zero elements in `matrix`.
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix A `tatami::Matrix`. 
 * @param row Whether to retrieve the contents of `matrix` by row, i.e., the output is a compressed sparse row matrix.
 * @param options Further options.
 *
 * @return Contents of the sparse matrix in compressed form, see `CompressedSparseContents`.
 *
 * The behavior of this function can be replicated by manually calling `count_compressed_sparse_non_zeros()` followed by `fill_compressed_sparse_contents()`.
 * This may be desirable for users who want to put the compressed sparse contents into pre-existing memory allocations.
 */
template<typename StoredValue_, typename StoredIndex_, typename StoredPointer_ = std::size_t, typename InputValue_, typename InputIndex_>
CompressedSparseContents<StoredValue_, StoredIndex_, StoredPointer_> retrieve_compressed_sparse_contents(
    const Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const RetrieveCompressedSparseContentsOptions& options
) {
    // We use size_t as the default pointer type here, as our output consists of vectors
    // with the default allocator, for which the size_type is unlikely to be bigger than size_t. 
    CompressedSparseContents<StoredValue_, StoredIndex_, StoredPointer_> output;
    auto& output_v = output.value;
    auto& output_i = output.index;
    auto& output_p = output.pointers;

    const InputIndex_ NR = matrix.nrow();
    const InputIndex_ NC = matrix.ncol();
    const InputIndex_ primary = (row ? NR : NC);
    const InputIndex_ secondary = (row ? NC : NR);

    const auto primary_p1 = sanisizer::sum<I<decltype(output_p.size())> >(attest_for_Index(primary), 1);
    output_p.resize(primary_p1);

    if (!options.two_pass) {
        // Doing a single fragmented run and then concatenating everything together.
        const auto frag = retrieve_fragmented_sparse_contents<InputValue_, InputIndex_>(
            matrix,
            row,
            [&]{
                RetrieveFragmentedSparseContentsOptions roptions;
                roptions.num_threads = options.num_threads;
                return roptions;
            }()
        );
        const auto& store_v = frag.value;
        const auto& store_i = frag.index;

        for (InputIndex_ p = 0; p < primary; ++p) {
            output_p[p + 1] = sanisizer::sum<StoredPointer_>(output_p[p], store_v[p].size());
        }

        output_v.reserve(output_p.back());
        output_i.reserve(output_p.back());
        for (InputIndex_ p = 0; p < primary; ++p) {
            output_v.insert(output_v.end(), store_v[p].begin(), store_v[p].end());
            output_i.insert(output_i.end(), store_i[p].begin(), store_i[p].end());
        }

    } else if (row == matrix.prefer_rows()) {
        // First pass to figure out how many non-zeros there are.
        count_compressed_sparse_non_zeros_consistent(matrix, primary, secondary, row, output_p.data() + 1, options.num_threads);
        for (InputIndex_ i = 1; i <= primary; ++i) {
            output_p[i] = sanisizer::sum<StoredPointer_>(output_p[i], output_p[i - 1]);
        }

        // Second pass to actually fill our vectors.
        sanisizer::resize(output_v, output_p.back());
        sanisizer::resize(output_i, output_p.back());
        fill_compressed_sparse_matrix_consistent(
            matrix,
            primary,
            secondary,
            row,
            output_p.data(),
            output_v.data(),
            output_i.data(),
            options.num_threads
        );

    } else {
        // First pass to figure out how many non-zeros there are.
        auto secondary_counts = create_container_of_Index_size<std::vector<InputIndex_> >(secondary);
        count_compressed_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, output_p.data() + 1, secondary_counts.data(), options.num_threads);
        for (InputIndex_ i = 1; i <= primary; ++i) {
            output_p[i] = sanisizer::sum<StoredPointer_>(output_p[i], output_p[i - 1]);
        }

        // Second pass to actually fill our vectors.
        sanisizer::resize(output_v, output_p.back());
        sanisizer::resize(output_i, output_p.back());
        fill_compressed_sparse_matrix_inconsistent(
            matrix,
            primary,
            secondary,
            row,
            output_p,
            secondary_counts,
            output_v.data(),
            output_i.data(),
            options.num_threads
        );
    }

    return output;
}

/**
 * @brief Options for `convert_to_compressed_sparse()`.
 */
struct ConvertToCompressedSparseOptions {
    /**
     * Whether to perform the retrieval in two passes.
     * This requires another pass through `matrix` but is more memory-efficient.
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
 * @tparam StoredPointer_ Integer type for the row/column pointers in the output.
 * This should be large enough to hold the number of non-zero elements in `matrix`.
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix A `tatami::Matrix`. 
 * @param row Whether to return a compressed sparse row matrix.
 * @param options Further options.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `matrix`.
 * If `row = true`, the matrix is in compressed sparse row format, otherwise it is compressed sparse column.
 */
template<
    typename Value_,
    typename Index_,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename StoredPointer_ = std::size_t,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_compressed_sparse(
    const Matrix<InputValue_, InputIndex_>& matrix,
    const bool row,
    const ConvertToCompressedSparseOptions& options
) {
    auto comp = retrieve_compressed_sparse_contents<StoredValue_, StoredIndex_, StoredPointer_>(
        matrix,
        row, 
        [&]{
            RetrieveCompressedSparseContentsOptions ropt;
            ropt.two_pass = options.two_pass;
            ropt.num_threads = options.num_threads;
            return ropt;
        }()
    );
    return std::shared_ptr<Matrix<Value_, Index_> >(
        new CompressedSparseMatrix<
            Value_, 
            Index_,
            std::vector<StoredValue_>,
            std::vector<StoredIndex_>,
            std::vector<StoredPointer_>
        >(
            matrix.nrow(), 
            matrix.ncol(), 
            std::move(comp.value), 
            std::move(comp.index), 
            std::move(comp.pointers),
            row,
            []{
                CompressedSparseMatrixOptions copt;
                copt.check = false; // no need for checks, as we guarantee correctness.
                return copt;
            }()
        )
    );
}

/**
 * @cond
 */
// Backwards compatbility.
template<typename Value_, typename Index_, typename Count_>
void count_compressed_sparse_non_zeros(const tatami::Matrix<Value_, Index_>* matrix, bool row, Count_* output, int threads) {
    return count_compressed_sparse_non_zeros(
        *matrix,
        row,
        output,
        [&]{
            CountCompressedSparseNonZerosOptions copt;
            copt.num_threads = threads;
            return copt;
        }()
    );
}

template<typename InputValue_, typename InputIndex_, typename Pointer_, typename StoredValue_, typename StoredIndex_>
void fill_compressed_sparse_contents(const tatami::Matrix<InputValue_, InputIndex_>* matrix,
    bool row,
    const Pointer_* pointers,
    StoredValue_* output_value,
    StoredIndex_* output_index,
    int threads)
{
    fill_compressed_sparse_contents(
        *matrix,
        row,
        pointers,
        output_value,
        output_index,
        [&]{
            FillCompressedSparseContentsOptions fopt;
            fopt.num_threads = threads;
            return fopt;
        }()
    );
}

template<typename StoredValue_, typename StoredIndex_, typename StoredPointer_ = std::size_t, typename InputValue_, typename InputIndex_>
CompressedSparseContents<StoredValue_, StoredIndex_, StoredPointer_> retrieve_compressed_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, bool row, bool two_pass, int threads = 1) {
    return retrieve_compressed_sparse_contents<StoredValue_, StoredIndex_>(
        *matrix,
        row,
        [&]{
            RetrieveCompressedSparseContentsOptions opt;
            opt.two_pass = two_pass;
            opt.num_threads = threads;
            return opt;
        }()
    );
}

template<typename Value_ = double, typename Index_ = int, typename StoredValue_ = Value_, typename StoredIndex_ = Index_, typename InputValue_, typename InputIndex_>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_compressed_sparse(const Matrix<InputValue_, InputIndex_>* matrix, bool row, bool two_pass = false, int threads = 1) {
    return convert_to_compressed_sparse<Value_, Index_, StoredValue_, StoredIndex_>(
        *matrix,
        row,
        [&]{
            ConvertToCompressedSparseOptions opt;
            opt.two_pass = two_pass;
            opt.num_threads = threads;
            return opt;
        }()
    );
}

template <bool row_, typename Value_, typename Index_, typename InputValue_, typename InputIndex_>
CompressedSparseContents<Value_, Index_, std::size_t> retrieve_compressed_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, bool two_pass, int threads = 1) {
    return retrieve_compressed_sparse_contents<Value_, Index_>(matrix, row_, two_pass, threads);
}

template <bool row_, typename Value_, typename Index_, typename StoredValue_ = Value_, typename StoredIndex_ = Index_, typename InputValue_, typename InputIndex_>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_compressed_sparse(const Matrix<InputValue_, InputIndex_>* matrix, bool two_pass = false, int threads = 1) {
    return convert_to_compressed_sparse<Value_, Index_, StoredValue_, StoredIndex_>(matrix, row_, two_pass, threads);
}
/**
 * @endcond
 */

}

#endif
