#ifndef TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H
#define TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H

#include <memory>
#include <vector>
#include <cstddef>

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
namespace convert_to_compressed_sparse_internal {

template<typename Value_, typename Index_, typename Count_>
void count_compressed_sparse_non_zeros_consistent(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Index_ primary,
    const Index_ secondary,
    const bool row,
    Count_* const output,
    const int threads)
{
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
    Count_* const output,
    const int threads
) {
    auto nz_counts = sanisizer::create<std::vector<std::vector<Count_> > >(threads - 1);
    for (auto& x : nz_counts) {
        x.resize(primary);
    }

    if (matrix.is_sparse()) {
        Options opt;
        opt.sparse_extract_value = false;
        opt.sparse_ordered_index = false;

        parallelize([&](const int t, const Index_ start, const Index_ length) -> void {
            auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
            auto buffer_i = create_container_of_Index_size<std::vector<Index_> >(primary);
            const auto my_counts = (t > 0 ? nz_counts[t - 1].data() : output);

            for (Index_ x = 0; x < length; ++x) {
                const auto range = wrk->fetch(NULL, buffer_i.data());
                for (Index_ i = 0; i < range.number; ++i) {
                    ++my_counts[range.index[i]];
                }
            }
        }, secondary, threads);

    } else {
        parallelize([&](const int t, const Index_ start, const Index_ length) -> void {
            auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<Value_> >(primary);
            const auto my_counts = (t > 0 ? nz_counts[t - 1].data() : output);

            for (Index_ x = 0; x < length; ++x) {
                const auto ptr = wrk->fetch(buffer_v.data());
                for (Index_ p = 0; p < primary; ++p) {
                    my_counts[p] += (ptr[p] != 0);
                }
            }
        }, secondary, threads);
    }

    for (auto& y : nz_counts) {
        for (Index_ p = 0; p < primary; ++p) {
            output[p] += y[p];
        }
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
    const int threads)
{
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
    const Pointer_* const pointers,
    StoredValue_* const output_value,
    StoredIndex_* const output_index,
    const int threads)
{
    if (matrix.is_sparse()) {
        Options opt;
        opt.sparse_ordered_index = false;

        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<true>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length, opt);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(length);
            auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(length);
            std::vector<Pointer_> offset_copy(pointers + start, pointers + start + length);

            for (InputIndex_ x = 0; x < secondary; ++x) {
                const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                for (InputIndex_ i = 0; i < range.number; ++i) {
                    auto& pos = offset_copy[range.index[i] - start];
                    output_value[pos] = range.value[i];
                    output_index[pos] = x; 
                    ++pos;
                }
            }
        }, primary, threads);

    } else {
        parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
            auto wrk = consecutive_extractor<false>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(length);
            std::vector<Pointer_> offset_copy(pointers + start, pointers + start + length);

            for (InputIndex_ x = 0; x < secondary; ++x) {
                const auto ptr = wrk->fetch(buffer_v.data());
                for (InputIndex_ p = 0; p < length; ++p) {
                    const auto val = ptr[p]; 
                    if (val != 0) {
                        auto& pos = offset_copy[p];
                        output_value[pos] = val;
                        output_index[pos] = x;
                        ++pos;
                    }
                }
            }
        }, primary, threads);
    }
}

}
/**
 * @endcond
 */

/**
 * @brief Options for `count_compressed_sparse_non_zeros()`.
 */
struct CountCompressedSparseNonZerosOptions {
    /**
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

/**
 * @tparam Value_ Type of value in the matrix.
 * @tparam Index_ Integer type of row/column index.
 * @tparam Count_ Integer type for the non-zero count.
 *
 * @param matrix A `tatami::Matrix`. 
 * @param row Whether to count structural non-zeros by row.
 * @param[out] output Pointer to an array of length equal to the number of rows (if `row = true`) or columns (otherwise) of `matrix`.
 * On output, this stores the number of structural non-zeros in each row (if `row = true`) or column (otherwise).
 * @param options Further options.
 *
 * For sparse `matrix`, all structural non-zero elements are reported, even if they have actual values of zero.
 * In contrast, for dense `matrix`, only the non-zero values are counted;
 * these are considered to be structural non-zeros upon conversion to a sparse matrix (e.g., in `fill_compressed_sparse_contents()`).
 */
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
    std::fill_n(output, primary, 0);

    if (row == matrix.prefer_rows()) {
        convert_to_compressed_sparse_internal::count_compressed_sparse_non_zeros_consistent(matrix, primary, secondary, row, output, options.num_threads);
    } else {
        convert_to_compressed_sparse_internal::count_compressed_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, output, options.num_threads);
    }
}

/**
 * @brief Options for `fill_compressed_sparse_contents()`.
 */
struct FillCompressedSparseContentsOptions {
    /**
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

/**
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam Pointer_ Integer type for the row/column pointers.
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix A `tatami::Matrix`. 
 * @param row Whether to fill `output_value` and `output_index` by row, i.e., the output represents a compressed sparse row matrix.
 * @param[in] pointers Pointer to an array of length greater than or equal to the number of rows (if `row = true`) or columns (otherwise) of `matrix`. 
 * Each entry contains the position of the start of each row/column in `output_value` and `output_index`.
 * This argument is equivalent to the array of pointers for the compressed sparse format (e.g., `CompressedSparseContents::pointers`),
 * and can be obtained by taking the cumulative sum of the per-row/column counts from `count_compressed_sparse_non_zeros()`.
 * @param[out] output_value Pointer to an array of length equal to the total number of structural non-zero elements.
 * On output, this is used to store the values of those elements in a compressed sparse format (e.g., `CompressedSparseContents::value`).
 * @param[out] output_index Pointer to an array of length equal to the total number of structural non-zero elements.
 * On output, this is used to store the row/column indices of those elements in a compressed sparse format (e.g., `CompressedSparseContents::index`).
 * @param options Further options.
 */
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
        convert_to_compressed_sparse_internal::fill_compressed_sparse_matrix_consistent(matrix, primary, secondary, row, pointers, output_value, output_index, options.num_threads);
    } else {
        convert_to_compressed_sparse_internal::fill_compressed_sparse_matrix_inconsistent(matrix, primary, secondary, row, pointers, output_value, output_index, options.num_threads);
    }
}

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
    const RetrieveCompressedSparseContentsOptions& options)
{
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
            output_p[p + 1] = output_p[p] + store_v[p].size();
        }

        output_v.reserve(output_p.back());
        output_i.reserve(output_p.back());
        for (InputIndex_ p = 0; p < primary; ++p) {
            output_v.insert(output_v.end(), store_v[p].begin(), store_v[p].end());
            output_i.insert(output_i.end(), store_i[p].begin(), store_i[p].end());
        }

    } else if (row == matrix.prefer_rows()) {
        // First pass to figure out how many non-zeros there are.
        convert_to_compressed_sparse_internal::count_compressed_sparse_non_zeros_consistent(matrix, primary, secondary, row, output_p.data() + 1, options.num_threads);
        for (InputIndex_ i = 1; i <= primary; ++i) {
            output_p[i] += output_p[i - 1];
        }

        // Second pass to actually fill our vectors.
        output_v.resize(output_p.back());
        output_i.resize(output_p.back());
        convert_to_compressed_sparse_internal::fill_compressed_sparse_matrix_consistent(
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
        convert_to_compressed_sparse_internal::count_compressed_sparse_non_zeros_inconsistent(matrix, primary, secondary, row, output_p.data() + 1, options.num_threads);
        for (InputIndex_ i = 1; i <= primary; ++i) {
            output_p[i] += output_p[i - 1];
        }

        // Second pass to actually fill our vectors.
        output_v.resize(output_p.back());
        output_i.resize(output_p.back());
        convert_to_compressed_sparse_internal::fill_compressed_sparse_matrix_inconsistent(
            matrix,
            primary,
            secondary,
            row,
            output_p.data(),
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
