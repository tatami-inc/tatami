#ifndef TATAMI_CONVERT_TO_SPARSE_H
#define TATAMI_CONVERT_TO_SPARSE_H

#include "../sparse/CompressedSparseMatrix.hpp"
#include "../sparse/FragmentedSparseMatrix.hpp"
#include "../stats/utils.hpp"

#include <memory>
#include <vector>
#include <deque>

/**
 * @file convert_to_sparse.hpp
 *
 * @brief Convert a matrix into a compressed sparse format.
 */

namespace tatami {

/**
 * @brief Fragmented sparse contents.
 *
 * @tparam Value_ Type of value in the matrix.
 * @tparam Index_ Type of row/column index.
 *
 * Check out `FragmentedSparseMatrix` for more details.
 */
template<typename Value_, typename Index_>
struct FragmentedSparseContents {
    /**
     * @cond
     */
    FragmentedSparseContents(size_t n) : value(n), index(n) {}
    /**
     * @endcond
     */

    /**
     * Vector of vectors containing the values for the structural non-zero elements.
     * Each inner vector corresponds to an element of the primary dimension and contains all values for that element.
     */
    std::vector<std::vector<Value_> > value;

    /**
     * Vector of vectors containing the secondary dimension indices for the structural non-zero elements.
     * Each inner vector corresponds to an element of the primary dimension and contains all indices for that element.
     * Each inner vector is of length equal to the corresponding entry of `values` and is guaranteed to be strictly increasing.
     */
    std::vector<std::vector<Index_> > index;
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
 * @return Contents of the sparse matrix in fragmented form, see `FragmentedSparseContents`.
 */
template <bool row_, typename Value_, typename Index_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_contents(const Matrix<InputValue_, InputIndex_>* incoming, int threads = 1) {
    Index_ NR = incoming->nrow();
    Index_ NC = incoming->ncol();
    Index_ primary = (row_ ? NR : NC);
    Index_ secondary = (row_ ? NC : NR);

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    auto& store_v = output.value;
    auto& store_i = output.index;

    if (row_ == incoming->prefer_rows()) {
        if (incoming->sparse()) {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<row_, true>(incoming, start, length);

                for (Index_ p = start, e = start + length; p < e; ++p) {
                    auto range = wrk->fetch(p, buffer_v.data(), buffer_i.data());
                    auto& sv = store_v[p];
                    auto& si = store_i[p];
                    sv.reserve(range.number);
                    si.reserve(range.number);

                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        if (*range.value) {
                            sv.push_back(*range.value);
                            si.push_back(*range.index);
                        }
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false>(incoming, start, length);

                // Special conversion from dense to save ourselves from having to make
                // indices that we aren't really interested in.
                for (Index_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    auto& sv = store_v[p];
                    auto& si = store_i[p];

                    for (Index_ s = 0; s < secondary; ++s, ++ptr) {
                        if (*ptr) {
                            sv.push_back(*ptr);
                            si.push_back(s);
                        }
                    }
                }
            }, primary, threads);
        }

    } else {
        // We iterate on the incoming matrix's preferred dimension, under the
        // assumption that it may be arbitrarily costly to extract in the
        // non-preferred dim; it is thus cheaper to do cache-unfriendly inserts
        // into the output buffers. 

        if (incoming->sparse()) {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(primary);
                std::vector<InputIndex_> buffer_i(primary);
                auto wrk = consecutive_extractor<!row_, true>(incoming, 0, secondary, start, length);

                for (Index_ s = 0; s < secondary; ++s) {
                    auto range = wrk->fetch(s, buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        if (*range.value) {
                            store_v[*range.index].push_back(*range.value);
                            store_i[*range.index].push_back(s);
                        }
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                auto wrk = consecutive_extractor<!row_, false>(incoming, 0, secondary, start, length);
                auto len = wrk->block_length;
                std::vector<InputData_> buffer_v(len);

                for (Index_ s = 0; s < secondary; ++s) {
                    auto ptr = wrk->fetch(s, buffer_v.data());
                    for (InputIndex_ p = 0; p < len; ++p, ++ptr) {
                        if (*ptr) {
                            store_v[p + start].push_back(*ptr);
                            store_i[p + start].push_back(s);
                        }
                    }
                }
            }, primary, threads);
        }
    }

    return output;
}

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
CompressedSparseContents<StoredValue_, StoredIndex_> retrieve_compressed_sparse_contents(const Matrix<InputValue_, InputIndex_>* incoming, bool two_pass, int threads = 1) {
    CompressedSparseContents<StoredValue_, StoredIndex_> output;
    auto& output_v = output.value;
    auto& output_i = output.index;

    Index_ primary = (row_ ? incoming->nrow() : incoming->ncol());
    Index_ secondary = (row_ ? incoming->ncol() : incoming->nrow());

    if (!two_pass) {
        auto frag = retrieve_fragmented_sparse_contents<row_, Value_, Index_>(incoming, threads);
        auto& store_v = frag.value;
        auto& store_i = frag.index;

        // Concatenating everything together.
        size_t total_size = 0;
        auto& indptrs = output.pointers;
        indptrs.reserve(static_cast<size_t>(primary) + 1);
        for (Index_ p = 0; p < primary; ++p) {
            total_size += store_v[p].size();
            indptrs[p + 1] = total_size;
        }

        output_v.reserve(total_size);
        output_i.reserve(total_size);
        for (Index_ p = 0; p < primary; ++p) {
            output_v.insert(output_v.end(), store_v[p].begin(), store_v[p].end());
            output_i.insert(output_i.end(), store_i[p].begin(), store_i[p].end());
        }

    } else if (row_ == incoming->prefer_rows()) {
        // First pass to figure out how many non-zeros there are.
        auto& offsets = output.pointers;
        offsets.resize(primary + 1);

        if (incoming->sparse()) {
            Options opt;
            opt.extract_sparse_value = false;
            opt.extract_sparse_index = false;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t t, Index_ start, Index_ length) -> void {
                auto wrk = consecutive_extractor<row_, true>(incoming, start, length, opt);
                for (Index_ p = start, e = start + length; p < e; ++p) {
                    auto range = wrk->fetch(p, NULL, NULL);
                    offsets[p + 1] = range.number;
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false>(incoming, start, length);
                for (Index_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    size_t count = 0;
                    for (Index_ s = 0; s < secondary; ++s, ++ptr) {
                        count += (*ptr != 0);
                    }
                    offsets[p + 1] = count;
                }
            }, primary, threads);
        }

        for (Index_ i = 1; i <= primary; ++i) {
            offsets[i] += offsets[i - 1];
        }
        output_v.resize(offsets.back());
        output_i.resize(offsets.back());

        // Second pass to actually fill our vectors.
        if (incoming->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<row_, true>(incoming, start, length, opt);

                for (Index_ p = start, e = start + length; p < e; ++p) {
                    // Resist the urge to `fetch_copy()` straight into
                    // store_v, as implementations may assume that they
                    // have the entire 'secondary' length to play with.
                    auto range = wrk->fetch(p, buffer_v.data(), buffer_i.data());
                    std::copy(range.value, range.value + range.number, store_v.data() + offsets[p]);
                    std::copy(range.index, range.index + range.number, store_i.data() + offsets[p]);
                }
            }, primary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false>(incoming, start, length);

                for (Index_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    size_t offset = offsets[p];

                    for (Index_ s = 0; s < secondary; ++s, ++ptr) {
                        if (*ptr != 0) {
                            store_v[offset] = *ptr;
                            store_i[offset] = s;
                            ++offset;
                        }
                    }
                }
            }, primary, threads);
        }

    } else {
        // First pass to figure out how many non-zeros there are.
        std::vector<std::vector<size_t> > nz_counts(nthreads);
        for (auto& x : nz_counts) {
            x.resize(static_cast<size_t>(primary));
        }

        if (incoming->sparse()) {
            Options opt;
            opt.extract_sparse_value = false;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t t, Index_ start, Index_ length) -> void {
                std::vector<InputIndex_> buffer_i(length);
                auto wrk = consecutive_extractor<!row_, true>(incoming, start, length, opt);
                auto& my_counts = nz_counts[t];
                for (Index_ s = start, end = start + length; s < end; ++s) {
                    auto range = wrk->fetch(s, NULL, buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.index) {
                        ++my_counts[*range.index + 1];
                    }
                }
            }, secondary, threads);

        } else {
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                auto len = wrk->block_length;
                std::vector<InputData_> buffer_v(len);
                auto wrk = consecutive_extractor<!row_, false>(incoming, start, length);
                auto& my_counts = nz_counts[t];
                for (Index_ s = start, end = start + length; s < end; ++s) {
                    auto ptr = wrk->fetch(s, buffer_v.data());
                    for (InputIndex_ p = 0; p < len; ++p, ++ptr) {
                        if (*ptr) {
                            ++my_counts[p + start + 1];
                        }
                    }
                }
            }, secondary, threads);
        }

        auto& offsets = output.pointers;
        offsets.swap(nz_counts.front()); // There better be at least 1 thread!
        for (int i = 1; i < threads; ++i) {
            auto& y = nz_counts[i];
            for (Index_ p = 0; p < primary; ++p) {
                offsets[p] += y[p];
            }
            y.clear();
            y.shrink_to_fit();
        }

        for (Index_ i = 1; i <= primary; ++i) {
            offsets[i] += offsets[i - 1];
        }
        output_v.resize(offsets.back());
        output_i.resize(offsets.back());

        // Second pass to actually fill our vectors.
        if (incoming->sparse()) {
            Options opt;
            opt.sparse_ordered_index = false;

            parallelize([&](size_t t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<!row_, true>(incoming, 0, secondary, start, length, opt);
                std::vector<size_t> offset_copy(offsets.begin() + start, offsets.begin() + start + length);

                for (size_t s = 0; s < secondary; ++p) {
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
            parallelize([&](size_t, Index_ start, Index_ length) -> void {
                std::vector<InputData_> buffer_v(secondary);
                auto wrk = consecutive_extractor<!row_, false>(incoming, 0, secondary, start, length);
                std::vector<size_t> offset_copy(offsets.begin() + start, offsets.begin() + start + length);

                for (size_t s = 0; s < secondary; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    size_t offset = offsets[p];

                    for (Index_ p = 0; p < length; ++p, ++ptr) {
                        if (*ptr != 0) {
                            auto& pos = offset_copy[p];
                            store_v[pos] = *ptr;
                            store_i[pos] = p + start;
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
 * @tparam InputMatrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`, possibly containing delayed operations.
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
inline std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const Matrix<InputValue_, InputIndex_>* incoming, int threads = 1) {
    if (threads == 0) {
        threads = 1; // for back-compatibility with calls using the old 'reserve' argument.
    }

        auto comp = retrieve_compressed_sparse_contents(incoming, threads);
        return std::shared_ptr<Matrix<Value_, Index_> >(
            new CompressedSparseMatrix<
                row_, 
                Value_, 
                Index_,
                decltype(output_v),
                decltype(output_i),
                decltype(indptrs)
            >(
                NR, 
                NC, 
                std::move(output_v), 
                std::move(output_i), 
                std::move(indptrs)
            )
        );
    }
}

/**
 * This overload makes it easier to control the desired output order when it is not known at compile time.
 *
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputMatrix_ Input matrix class, most typically a `tatami::Matrix`.
 *
 * @param incoming Pointer to a `tatami::Matrix`.
 * @param order Ordering of values in the output matrix - compressed sparse row (0) or column (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::CompressedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    class InputMatrix_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_sparse(const InputMatrix_* incoming, int order, int threads = 1) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows());
    }
    if (order == 0) {
        return convert_to_sparse<true, Value_, Index_, StoredValue_, StoredIndex_, InputMatrix_>(incoming, threads);
    } else {
        return convert_to_sparse<false, Value_, Index_, StoredValue_, StoredIndex_, InputMatrix_>(incoming, threads);
    }
}

}

#endif
