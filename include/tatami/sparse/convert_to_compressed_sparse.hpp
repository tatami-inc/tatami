#ifndef TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H
#define TATAMI_CONVERT_TO_COMPRESSED_SPARSE_H

#include "CompressedSparseMatrix.hpp"
#include "convert_to_fragmented_sparse.hpp"
#include "../stats/utils.hpp"

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
 * Check out `CompressedSparseMatrix` for more details.
 */
template<typename Value_, typename Index_>
struct CompressedSparseContents {
    /**
     * Vector containing values for the structural non-zero elements in a compressed sparse format.
     */
    std::vector<Value_> value;

    /**
     * Vector containing the secondary dimension indices for the structural non-zero elements in a compressed sparse format.
     */
    std::vector<Index_> index;

    /**
     * Vector containing the pointers for each primary dimension element in a compressed sparse format.
     */
    std::vector<size_t> pointers;
};

/**
 * @tparam row_ Whether to retrieve contents by row.
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`. 
 * @param threads Number of threads to use.
 *
 * @return Contents of the sparse matrix in compressed form, see `FragmentedSparseContents`.
 */
template <bool row_, typename Value_, typename Index_, typename InputValue_, typename InputIndex_>
CompressedSparseContents<Value_, Index_> retrieve_compressed_sparse_contents(const Matrix<InputValue_, InputIndex_>* incoming, bool two_pass, int threads = 1) {
    CompressedSparseContents<Value_, Index_> output;
    auto& output_v = output.value;
    auto& output_i = output.index;
    auto& output_p = output.pointers;

    InputIndex_ NR = incoming->nrow();
    InputIndex_ NC = incoming->ncol();
    InputIndex_ primary = (row_ ? NR : NC);
    InputIndex_ secondary = (row_ ? NC : NR);

    if (!two_pass) {
        // Doing a single fragmented run and then concatenating everything together.
        auto frag = retrieve_fragmented_sparse_contents<row_, InputValue_, InputIndex_>(incoming, threads);
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

    } else if (row_ == incoming->prefer_rows()) {
        // First pass to figure out how many non-zeros there are.
        output_p.resize(static_cast<size_t>(primary) + 1);

        if (incoming->sparse()) {
            Options opt;
            opt.sparse_extract_value = false;
            opt.sparse_extract_index = false;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<row_, true>(incoming, start, length, opt);
                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    auto range = wrk->fetch(p, NULL, NULL);
                    output_p[p + 1] = range.number;
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false>(incoming, start, length);
                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
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
        if (incoming->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<row_, true>(incoming, start, length, opt);

                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    // Resist the urge to `fetch_copy()` straight into
                    // store_v, as implementations may assume that they
                    // have the entire 'secondary' length to play with.
                    auto range = wrk->fetch(p, buffer_v.data(), buffer_i.data());
                    std::copy(range.value, range.value + range.number, output_v.data() + output_p[p]);
                    std::copy(range.index, range.index + range.number, output_i.data() + output_p[p]);
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false>(incoming, start, length);

                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    size_t offset = output_p[p];

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

        if (incoming->sparse()) {
            Options opt;
            opt.sparse_extract_value = false;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputIndex_> buffer_i(primary);
                auto wrk = consecutive_extractor<!row_, true>(incoming, start, length, opt);
                auto& my_counts = nz_counts[t];

                for (InputIndex_ s = start, end = start + length; s < end; ++s) {
                    auto range = wrk->fetch(s, NULL, buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.index) {
                        ++my_counts[*range.index + 1];
                    }
                }
            }, secondary, threads);

        } else {
            parallelize([&](size_t t, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<!row_, false>(incoming, start, length);
                std::vector<InputValue_> buffer_v(primary);
                auto& my_counts = nz_counts[t];

                for (InputIndex_ s = start, end = start + length; s < end; ++s) {
                    auto ptr = wrk->fetch(s, buffer_v.data());
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
        if (incoming->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(length);
                std::vector<InputIndex_> buffer_i(length);
                auto wrk = consecutive_extractor<!row_, true>(incoming, static_cast<InputIndex_>(0), secondary, start, length, opt);
                std::vector<size_t> offset_copy(output_p.begin() + start, output_p.begin() + start + length);

                for (InputIndex_ s = 0; s < secondary; ++s) {
                    auto range = wrk->fetch(s, buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        auto& pos = offset_copy[*(range.index) - start];
                        output_v[pos] = *(range.value);
                        output_i[pos] = s; 
                        ++pos;
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(length);
                auto wrk = consecutive_extractor<!row_, false>(incoming, static_cast<InputIndex_>(0), secondary, start, length);
                std::vector<size_t> offset_copy(output_p.begin() + start, output_p.begin() + start + length);

                for (InputIndex_ s = 0; s < secondary; ++s) {
                    auto ptr = wrk->fetch(s, buffer_v.data());
                    for (InputIndex_ p = 0; p < length; ++p, ++ptr) {
                        if (*ptr != 0) {
                            auto& pos = offset_copy[p];
                            output_v[pos] = *ptr;
                            output_i[pos] = s;
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
 * @tparam row_ Whether to return a compressed sparse row matrix.
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`, possibly containing delayed operations.
 * @param two_pass Whether to use a two-pass strategy that reduces peak memory usage at the cost of speed.
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row_ = true`, the matrix is compressed sparse row, otherwise it is compressed sparse column.
 */
template <bool row_, 
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_compressed_sparse(const Matrix<InputValue_, InputIndex_>* incoming, bool two_pass = false, int threads = 1) {
    auto comp = retrieve_compressed_sparse_contents<row_, StoredValue_, StoredIndex_>(incoming, two_pass, threads);
    return std::shared_ptr<Matrix<Value_, Index_> >(
        new CompressedSparseMatrix<
            row_, 
            Value_, 
            Index_,
            std::vector<StoredValue_>,
            std::vector<StoredIndex_>,
            std::vector<size_t>
        >(
            incoming->nrow(), 
            incoming->ncol(), 
            std::move(comp.value), 
            std::move(comp.index), 
            std::move(comp.pointers),
            false // no need for checks, as we guarantee correctness.
        )
    );
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output matrix - compressed sparse row (0) or column (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 * @param two_pass Whether to use a two-pass strategy that reduces peak memory usage at the cost of speed.
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_compressed_sparse(const Matrix<InputValue_, InputIndex_>* incoming, int order, bool two_pass = false, int threads = 1) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows());
    }
    if (order == 0) {
        return convert_to_compressed_sparse<true, Value_, Index_, StoredValue_, StoredIndex_>(incoming, two_pass, threads);
    } else {
        return convert_to_compressed_sparse<false, Value_, Index_, StoredValue_, StoredIndex_>(incoming, two_pass, threads);
    }
}

}

#endif
