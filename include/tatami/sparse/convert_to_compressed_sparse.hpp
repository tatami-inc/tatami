#ifndef TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H
#define TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H

#include "CompressedSparseMatrix.hpp"
#include "convert_to_fragmented_sparse.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/consecutive_extractor.hpp"

#include <memory>
#include <vector>

/**
 * @file convert_to_compressed_sparse.hpp
 *
 * @brief Convert a matrix into a compressed sparse format.
 */

namespace tatami {

/**
 * @brief Compressed sparse contents.
 *
 * @tparam Value_ Type of value in the matrix.
 * @tparam Index_ Type of row/column index.
 *
 * The "primary" dimension is the one that is used to create the pointers for the compressed sparse format, while the other dimension is defined as the "secondary" dimension.
 * For example, the rows would be the primary dimension in a compressed sparse row matrix.
 */
template<typename Value_, typename Index_>
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
    std::vector<size_t> pointers;
};

/**
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix Pointer to a `tatami::Matrix`. 
 * @param row Whether to retrieve the contents of `matrix` by row, i.e., the output is a compressed sparse row matrix.
 * @param two_pass Whether to perform the retrieval in two passes.
 * This requires another pass through `matrix` but is more memory-efficient.
 * @param threads Number of threads to use.
 *
 * @return Contents of the sparse matrix in compressed form, see `CompressedSparseContents`.
 */
template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
CompressedSparseContents<StoredValue_, StoredIndex_> retrieve_compressed_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, bool row, bool two_pass, int threads = 1) {
    CompressedSparseContents<StoredValue_, StoredIndex_> output;
    auto& output_v = output.value;
    auto& output_i = output.index;
    auto& output_p = output.pointers;

    InputIndex_ NR = matrix->nrow();
    InputIndex_ NC = matrix->ncol();
    InputIndex_ primary = (row ? NR : NC);
    InputIndex_ secondary = (row ? NC : NR);

    if (!two_pass) {
        // Doing a single fragmented run and then concatenating everything together.
        auto frag = retrieve_fragmented_sparse_contents<InputValue_, InputIndex_>(matrix, row, threads);
        const auto& store_v = frag.value;
        const auto& store_i = frag.index;

        output_p.resize(static_cast<size_t>(primary) + 1);
        for (InputIndex_ p = 0; p < primary; ++p) {
            output_p[p + 1] = output_p[p] + store_v[p].size();
        }

        output_v.reserve(output_p.back());
        output_i.reserve(output_p.back());
        for (InputIndex_ p = 0; p < primary; ++p) {
            output_v.insert(output_v.end(), store_v[p].begin(), store_v[p].end());
            output_i.insert(output_i.end(), store_i[p].begin(), store_i[p].end());
        }

    } else if (row == matrix->prefer_rows()) {
        // First pass to figure out how many non-zeros there are.
        output_p.resize(static_cast<size_t>(primary) + 1);

        if (matrix->is_sparse()) {
            Options opt;
            opt.sparse_extract_value = false;
            opt.sparse_extract_index = false;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<true>(matrix, row, start, length, opt);
                for (InputIndex_ x = 0; x < length; ++x) {
                    auto range = wrk->fetch(NULL, NULL);
                    output_p[start + x + 1] = range.number;
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                auto wrk = consecutive_extractor<false>(matrix, row, start, length);
                for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                    auto ptr = wrk->fetch(buffer_v.data());
                    size_t count = 0;
                    for (InputIndex_ s = 0; s < secondary; ++s, ++ptr) {
                        count += (*ptr != 0);
                    }
                    output_p[p + 1] = count;
                }
            }, primary, threads);
        }

        for (InputIndex_ i = 1; i <= primary; ++i) {
            output_p[i] += output_p[i - 1];
        }
        output_v.resize(output_p.back());
        output_i.resize(output_p.back());

        // Second pass to actually fill our vectors.
        if (matrix->is_sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<true>(matrix, row, start, length, opt);

                for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                    // Resist the urge to `fetch()` straight into 'output_p'
                    // and 'output_i', as implementations may assume that they
                    // have the entire 'length' length to play with, and the
                    // output vectors only have whatever is allocated from the
                    // first pass (which might be nothing for an all-zero matrix).
                    auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    auto offset = output_p[p];
                    std::copy_n(range.value, range.number, output_v.data() + offset);
                    std::copy_n(range.index, range.number, output_i.data() + offset);
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                auto wrk = consecutive_extractor<false>(matrix, row, start, length);

                for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                    auto ptr = wrk->fetch(buffer_v.data());
                    auto offset = output_p[p];
                    for (InputIndex_ s = 0; s < secondary; ++s, ++ptr) {
                        if (*ptr != 0) {
                            output_v[offset] = *ptr;
                            output_i[offset] = s;
                            ++offset;
                        }
                    }
                }
            }, primary, threads);
        }

    } else {
        // First pass to figure out how many non-zeros there are.
        std::vector<std::vector<size_t> > nz_counts(threads);
        for (auto& x : nz_counts) {
            x.resize(static_cast<size_t>(primary) + 1);
        }

        if (matrix->is_sparse()) {
            Options opt;
            opt.sparse_extract_value = false;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputIndex_> buffer_i(primary);
                auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
                auto& my_counts = nz_counts[t];

                for (InputIndex_ x = 0; x < length; ++x) {
                    auto range = wrk->fetch(NULL, buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.index) {
                        ++my_counts[*range.index + 1];
                    }
                }
            }, secondary, threads);

        } else {
            parallelize([&](size_t t, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
                std::vector<InputValue_> buffer_v(primary);
                auto& my_counts = nz_counts[t];

                for (InputIndex_ x = 0; x < length; ++x) {
                    auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = 0; p < primary; ++p, ++ptr) {
                        if (*ptr) {
                            ++my_counts[p + 1];
                        }
                    }
                }
            }, secondary, threads);
        }

        output_p.swap(nz_counts.front()); // There better be at least 1 thread!
        for (int i = 1; i < threads; ++i) {
            auto& y = nz_counts[i];
            for (InputIndex_ p = 0; p <= primary; ++p) {
                output_p[p] += y[p];
            }
            y.clear();
            y.shrink_to_fit();
        }

        for (InputIndex_ i = 1; i <= primary; ++i) {
            output_p[i] += output_p[i - 1];
        }
        output_v.resize(output_p.back());
        output_i.resize(output_p.back());

        // Second pass to actually fill our vectors.
        if (matrix->is_sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(length);
                std::vector<InputIndex_> buffer_i(length);
                auto wrk = consecutive_extractor<true>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length, opt);
                std::vector<size_t> offset_copy(output_p.begin() + start, output_p.begin() + start + length);

                for (InputIndex_ x = 0; x < secondary; ++x) {
                    auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        auto& pos = offset_copy[*(range.index) - start];
                        output_v[pos] = *(range.value);
                        output_i[pos] = x; 
                        ++pos;
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(length);
                auto wrk = consecutive_extractor<false>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length);
                std::vector<size_t> offset_copy(output_p.begin() + start, output_p.begin() + start + length);

                for (InputIndex_ x = 0; x < secondary; ++x) {
                    auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = 0; p < length; ++p, ++ptr) {
                        if (*ptr != 0) {
                            auto& pos = offset_copy[p];
                            output_v[pos] = *ptr;
                            output_i[pos] = x;
                            ++pos;
                        }
                    }
                }
            }, primary, threads);
        }
    }

    return output;
}

/**
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix Pointer to a `tatami::Matrix`, possibly containing delayed operations.
 * @param row Whether to return a compressed sparse row matrix.
 * @param two_pass Whether to use a two-pass strategy that reduces peak memory usage at the cost of speed.
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `matrix`.
 * If `row = true`, the matrix is in compressed sparse row format, otherwise it is compressed sparse column.
 */
template<
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_compressed_sparse(const Matrix<InputValue_, InputIndex_>* matrix, bool row, bool two_pass = false, int threads = 1) {
    auto comp = retrieve_compressed_sparse_contents<StoredValue_, StoredIndex_>(matrix, row, two_pass, threads);
    return std::shared_ptr<Matrix<Value_, Index_> >(
        new CompressedSparseMatrix<
            Value_, 
            Index_,
            std::vector<StoredValue_>,
            std::vector<StoredIndex_>,
            std::vector<size_t>
        >(
            matrix->nrow(), 
            matrix->ncol(), 
            std::move(comp.value), 
            std::move(comp.index), 
            std::move(comp.pointers),
            row,
            false // no need for checks, as we guarantee correctness.
        )
    );
}

/**
 * @cond
 */
// Backwards compatbility.
template <bool row_, typename Value_, typename Index_, typename InputValue_, typename InputIndex_>
CompressedSparseContents<Value_, Index_> retrieve_compressed_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, bool two_pass, int threads = 1) {
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
