#ifndef TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H
#define TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H

#include "FragmentedSparseMatrix.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/consecutive_extractor.hpp"

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
 * The "primary" dimension is the one that is used to organize non-zero elements into vectors, while the other dimension is defined as the "secondary" dimension.
 * For example, the rows would be the primary dimension in a fragmented sparse row matrix.
 * (Check out `FragmentedSparseMatrix` for more details.)
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
 * @tparam StoredValue_ Type of data values to be stored in the output.
 * @tparam StoredIndex_ Integer type for storing the indices in the output. 
 * @tparam InputValue_ Type of data values in the input interface.
 * @tparam InputIndex_ Integer type for indices in the input interface.
 *
 * @param matrix Pointer to a `tatami::Matrix`. 
 * @param row Whether to retrieve the contents of `matrix` by row, i.e., the output is a fragmented sparse row matrix.
 * @param threads Number of threads to use, for parallelization with `parallelize()`.
 *
 * @return Contents of the sparse matrix in fragmented form, see `FragmentedSparseContents`.
 */
template<typename StoredValue_, typename StoredIndex_, typename InputValue_, typename InputIndex_>
FragmentedSparseContents<StoredValue_, StoredIndex_> retrieve_fragmented_sparse_contents(const Matrix<InputValue_, InputIndex_>* matrix, bool row, int threads = 1) {
    InputIndex_ NR = matrix->nrow();
    InputIndex_ NC = matrix->ncol();
    InputIndex_ primary = (row ? NR : NC);
    InputIndex_ secondary = (row ? NC : NR);

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    auto& store_v = output.value;
    auto& store_i = output.index;

    if (row == matrix->prefer_rows()) {
        if (matrix->is_sparse()) {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                std::vector<InputIndex_> buffer_i(secondary);
                auto wrk = consecutive_extractor<true>(matrix, row, start, length);

                for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                    auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
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
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(secondary);
                auto wrk = consecutive_extractor<false>(matrix, row, start, length);

                // Special conversion from dense to save ourselves from having to make
                // indices that we aren't really interested in.
                for (InputIndex_ p = start, pe = start + length; p < pe; ++p) {
                    auto ptr = wrk->fetch(buffer_v.data());
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
        // We iterate on the matrix matrix's preferred dimension, under the
        // assumption that it may be arbitrarily costly to extract in the
        // non-preferred dim; it is thus cheaper to do cache-unfriendly inserts
        // into the output buffers. 

        if (matrix->is_sparse()) {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                std::vector<InputValue_> buffer_v(primary);
                std::vector<InputIndex_> buffer_i(primary);
                auto wrk = consecutive_extractor<true>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length);

                for (InputIndex_ x = 0; x < secondary; ++x) {
                    auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        if (*range.value) {
                            store_v[*range.index].push_back(*range.value);
                            store_i[*range.index].push_back(x);
                        }
                    }
                }
            }, primary, threads);

        } else {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<false>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length);
                std::vector<InputValue_> buffer_v(length);

                for (InputIndex_ x = 0; x < secondary; ++x) {
                    auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = start, pe = start + length; p < pe; ++p, ++ptr) {
                        if (*ptr) {
                            store_v[p].push_back(*ptr);
                            store_i[p].push_back(x);
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
 * @param row Whether to return a fragmented sparse row matrix.
 * @param threads Number of threads to use, for parallelization with `parallelize()`.
 *
 * @return A pointer to a new `tatami::FragmentedSparseMatrix`, with the same dimensions and type as the matrix referenced by `matrix`.
 * If `row = true`, the matrix is in fragmented sparse row format, otherwise it is fragmented sparse column.
 */
template<typename Value_, typename Index_, typename StoredValue_ = Value_, typename StoredIndex_ = Index_, typename InputValue_, typename InputIndex_>
std::shared_ptr<Matrix<Value_, Index_> > convert_to_fragmented_sparse(const Matrix<InputValue_, InputIndex_>* matrix, bool row, int threads = 1) {
    auto frag = retrieve_fragmented_sparse_contents<StoredValue_, StoredIndex_>(matrix, row, threads);
    return std::shared_ptr<Matrix<Value_, Index_> >(
        new FragmentedSparseMatrix<
            Value_, 
            Index_,
            std::vector<std::vector<StoredValue_> >,
            std::vector<std::vector<StoredIndex_> >
        >(
            matrix->nrow(), 
            matrix->ncol(), 
            std::move(frag.value), 
            std::move(frag.index),
            row, 
            false // no need for checks, as we guarantee correctness.
        )
    );
}

/**
 * @cond
 */
// Backwards compatbility.
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
