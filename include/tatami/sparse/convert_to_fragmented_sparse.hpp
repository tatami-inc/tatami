#ifndef TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H
#define TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H

#include "FragmentedSparseMatrix.hpp"
#include "../stats/utils.hpp"

#include <memory>
#include <vector>

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
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`. 
 * @param threads Number of threads to use.
 *
 * @return Contents of the sparse matrix in fragmented form, see `FragmentedSparseContents`.
 */
template <bool row_, typename Value_, typename Index_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<Value_, Index_> retrieve_fragmented_sparse_contents(const Matrix<InputValue_, InputIndex_>* incoming, int threads = 1) {
    InputIndex_ NR = incoming->nrow();
    InputIndex_ NC = incoming->ncol();
    InputIndex_ primary = (row_ ? NR : NC);
    InputIndex_ secondary = (row_ ? NC : NR);

    FragmentedSparseContents<Value_, Index_> output(primary);
    auto& store_v = output.value;
    auto& store_i = output.index;

    if (row_ == incoming->prefer_rows()) {
        if (incoming->sparse()) {
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<row_, true>(incoming, start, length);

                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
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
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                auto wrk = consecutive_extractor<row_, false>(incoming, start, length);

                // Special conversion from dense to save ourselves from having to make
                // indices that we aren't really interested in.
                for (InputIndex_ p = start, e = start + length; p < e; ++p) {
                    auto ptr = wrk->fetch(p, buffer_v.data());
                    auto& sv = store_v[p];
                    auto& si = store_i[p];

                    for (InputIndex_ s = 0; s < secondary; ++s, ++ptr) {
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
            parallelize([&](size_t, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(primary);
                std::vector<InputIndex_> buffer_i(primary);
                auto wrk = consecutive_extractor<!row_, true>(incoming, static_cast<InputIndex_>(0), secondary, start, length);

                for (InputIndex_ s = 0; s < secondary; ++s) {
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
                auto wrk = consecutive_extractor<!row_, false>(incoming, static_cast<InputIndex_>(0), secondary, start, length);
                auto len = wrk->block_length;
                std::vector<InputValue_> buffer_v(len);

                for (InputIndex_ s = 0; s < secondary; ++s) {
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
 * @tparam row_ Whether to return a fragmented sparse row matrix.
 * @tparam Value_ Type of data values in the output interface.
 * @tparam Index_ Integer type for the indices in the output interface.
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param incoming Pointer to a `tatami::Matrix`, possibly containing delayed operations.
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::FragmentedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 * If `row_ = true`, the matrix is fragmented sparse row, otherwise it is fragmented sparse column.
 */
template <bool row_, 
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_fragmented_sparse(const Matrix<InputValue_, InputIndex_>* incoming, int threads = 1) {
    auto frag = retrieve_fragmented_sparse_contents<row_, StoredValue_, StoredIndex_>(incoming, threads);
    return std::shared_ptr<Matrix<Value_, Index_> >(
        new FragmentedSparseMatrix<
            row_, 
            Value_, 
            Index_,
            std::vector<std::vector<StoredValue_> >,
            std::vector<std::vector<StoredIndex_> >
        >(
            incoming->nrow(), 
            incoming->ncol(), 
            std::move(frag.value), 
            std::move(frag.index),
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
 * @param order Ordering of values in the output matrix - fragmented sparse row (0) or column (1).
 * If set to -1, the ordering is chosen based on `tatami::Matrix::prefer_rows()`. 
 * @param threads Number of threads to use.
 *
 * @return A pointer to a new `tatami::FragmentedSparseMatrix`, with the same dimensions and type as the matrix referenced by `incoming`.
 */
template <
    typename Value_ = double,
    typename Index_ = int,
    typename StoredValue_ = Value_,
    typename StoredIndex_ = Index_,
    typename InputValue_,
    typename InputIndex_
>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_fragmented_sparse(const Matrix<InputValue_, InputIndex_>* incoming, int order, int threads = 1) {
    if (order < 0) {
        order = static_cast<int>(!incoming->prefer_rows());
    }
    if (order == 0) {
        return convert_to_fragmented_sparse<true, Value_, Index_, StoredValue_, StoredIndex_, InputValue_, InputIndex_>(incoming, threads);
    } else {
        return convert_to_fragmented_sparse<false, Value_, Index_, StoredValue_, StoredIndex_, InputValue_, InputIndex_>(incoming, threads);
    }
}

}

#endif
