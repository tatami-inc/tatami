#ifndef TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H
#define TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H

#include <memory>
#include <vector>
#include <cstddef>
#include <optional>

#include "CompressedSparseMatrix.hpp"
#include "convert_to_fragmented_sparse.hpp"
#include "convert_to_sparse_utils.hpp"

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
        std::fill_n(output, primary, 0);
        count_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, output, options.num_threads);
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

template<typename InputValue_, typename InputIndex_, typename Pointer_, typename StoredValue_, typename StoredIndex_>
void fill_compressed_sparse_matrix_inconsistent(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const InputIndex_ primary,
    const InputIndex_ secondary,
    const bool row,
    const Pointer_* const output_ptrs, 
    StoredValue_* const output_value,
    StoredIndex_* const output_index,
    std::optional<CountNonZerosPerThread<InputIndex_, Pointer_> >& per_thread // see count_sparse_non_zeros_inconsistent() in convert_to_sparse_utils.hpp.
) {
    const bool is_sparse = matrix.is_sparse();
    if (per_thread.has_value()) {
        // Transforming the per-thread counts into per-thread starting offsets.
        auto& offsets = per_thread->counts;
        Pointer_ accumulant = 0;
        static_assert(std::is_same<I<decltype(per_thread->counts[0][0])>, Pointer_>::value); // confirm that the accumulant assignment won't overflow.
        for (InputIndex_ i = 0; i < primary; ++i) {
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

                // We're going to completely ignore the potential for false sharing here.
                // False sharing would only be a risk for very fat/thin matrices (depending on row= and num_threads=),
                // and if the input matrix is sparse, this further lowers the chance of contention between threads.
                // The alternative would be to allocate a per-thread buffer to store all of the values,
                // but in that case, we might as well use a one-pass algorithm.
                if (is_sparse) {
                    Options opt;
                    opt.sparse_ordered_index = false;
                    auto wrk = consecutive_extractor<true>(matrix, !row, actual_start, actual_length, opt);
                    auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                    auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
                    for (InputIndex_ x = 0; x < actual_length; ++x) {
                        const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                        for (InputIndex_ i = 0; i < range.number; ++i) {
                            auto& pos = offsets[range.index[i]];
                            output_value[pos] = range.value[i];
                            output_index[pos] = x + actual_start;
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
                                output_value[pos] = val;
                                output_index[pos] = x + actual_start; 
                                ++pos;
                            }
                        }
                    }
                }
            }
        }, per_thread->counts.size(), per_thread->counts.size());

    } else {
        std::vector<Pointer_> offsets(output_ptrs, output_ptrs + primary);

        if (is_sparse){ 
            Options opt;
            opt.sparse_ordered_index = false;
            auto wrk = consecutive_extractor<true>(matrix, !row, static_cast<InputIndex_>(0), secondary, opt);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
            auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
            for (InputIndex_ s = 0; s < secondary; ++s) {
                const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    auto& pos = offsets[range.index[i]];
                    output_value[pos] = range.value[i];
                    output_index[pos] = s;
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
                        output_value[pos] = val;
                        output_index[pos] = s; 
                        ++pos;
                    }
                }
            }
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
        std::optional<CountNonZerosPerThread<InputIndex_, Pointer_> > empty; // Force it to be single-threaded as we don't have the per-worker pointers to parallelize effectively.
        fill_compressed_sparse_matrix_inconsistent(
            matrix,
            primary,
            secondary,
            row,
            pointers,
            output_value,
            output_index,
            empty
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

    output_p.resize(sanisizer::sum<I<decltype(output_p.size())> >(attest_for_Index(primary), 1));

    if (!options.two_pass) {
        // In the one-pass strategy, we load matrix contents along the preferred dimension first, then we transform it in serial.
        std::vector<std::vector<InputValue_> > store_v;
        std::vector<std::vector<InputIndex_> > store_i;
        auto original_ranges = extract_sparse_matrix(matrix, store_v, store_i, options.num_threads);

        const bool use_rows = matrix.prefer_rows();
        if (use_rows == row) {
            // Now concatenating everything together, if we're fortunate enough that the dimensions are consistent.
            for (InputIndex_ p = 0; p < primary; ++p) {
                output_p[p + 1] = sanisizer::sum<StoredPointer_>(output_p[p], original_ranges[p].number);
            }

            output_v.reserve(output_p.back());
            output_i.reserve(output_p.back());
            for (InputIndex_ p = 0; p < primary; ++p) {
                output_v.insert(output_v.end(), original_ranges[p].value, original_ranges[p].value + original_ranges[p].number);
                output_i.insert(output_i.end(), original_ranges[p].index, original_ranges[p].index + original_ranges[p].number);
            }

        } else {
            // Otherwise we need to compute the non-zeros on the inconsistent dimension before populating the output vectors.
            for (InputIndex_ s = 0; s < secondary; ++s) {
                const auto& range = original_ranges[s];
                for (InputIndex_ x = 0; x < range.number; ++x) {
                    output_p[range.index[x] + 1] += 1; // increments are safe at this point: p < primary and the total count must be less than 'secondary'.
                }
            }
            for (InputIndex_ p = 0; p < primary; ++p) {
                output_p[p + 1] = sanisizer::sum<StoredPointer_>(output_p[p + 1], output_p[p]);
            }

            sanisizer::resize(output_v, output_p.back());
            sanisizer::resize(output_i, output_p.back());
            std::vector<StoredPointer_> offsets(output_p.begin(), output_p.begin() + primary);
            for (InputIndex_ s = 0; s < secondary; ++s) {
                const auto& range = original_ranges[s];
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    auto& pos = offsets[range.index[i]];
                    output_v[pos] = range.value[i];
                    output_i[pos] = s;
                    ++pos;
                }
            }
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
        auto per_thread = count_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, output_p.data() + 1, options.num_threads);
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
            output_p.data(),
            output_v.data(),
            output_i.data(),
            per_thread
        );
    }

    return output;
}

/**
 * @brief Options for `convert_to_compressed_sparse()`.
 */
struct ConvertToCompressedSparseOptions {
    /**
     * Whether to perform the conversion in two passes.
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
