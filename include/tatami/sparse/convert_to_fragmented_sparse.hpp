#ifndef TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H
#define TATAMI_CONVERT_TO_FRAGMENTED_SPARSE_H

#include "FragmentedSparseMatrix.hpp"
#include "../utils/parallelize.hpp"
#include "../utils/consecutive_extractor.hpp"
#include "../utils/Index_to_container.hpp"

#include <memory>
#include <vector>
#include <cstddef>

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
        value(cast_Index_to_container_size<decltype(value)>(n)),
        index(cast_Index_to_container_size<decltype(index)>(n))
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
     * Number of threads to use, for parallelization with `parallelize()`.
     */
    int num_threads = 1;
};

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
    bool row,
    const RetrieveFragmentedSparseContentsOptions& options)
{
    InputIndex_ NR = matrix.nrow();
    InputIndex_ NC = matrix.ncol();
    InputIndex_ primary = (row ? NR : NC);
    InputIndex_ secondary = (row ? NC : NR);

    FragmentedSparseContents<StoredValue_, StoredIndex_> output(primary);
    auto& store_v = output.value;
    auto& store_i = output.index;

    if (row == matrix.prefer_rows()) {
        if (matrix.is_sparse()) {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<true>(matrix, row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(secondary);
                auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(secondary);

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
            }, primary, options.num_threads);

        } else {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<false>(matrix, row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(secondary);

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
            }, primary, options.num_threads);
        }

    } else {
        // We iterate on the matrix matrix's preferred dimension, under the
        // assumption that it may be arbitrarily costly to extract in the
        // non-preferred dim; it is thus cheaper to do cache-unfriendly inserts
        // into the output buffers. 

        if (matrix.is_sparse()) {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<true>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);

                for (InputIndex_ x = 0; x < secondary; ++x) {
                    auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    for (InputIndex_ i = 0; i < range.number; ++i, ++range.value, ++range.index) {
                        if (*range.value) {
                            store_v[*range.index].push_back(*range.value);
                            store_i[*range.index].push_back(x);
                        }
                    }
                }
            }, primary, options.num_threads);

        } else {
            parallelize([&](int, InputIndex_ start, InputIndex_ length) -> void {
                auto wrk = consecutive_extractor<false>(matrix, !row, static_cast<InputIndex_>(0), secondary, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(length);

                for (InputIndex_ x = 0; x < secondary; ++x) {
                    auto ptr = wrk->fetch(buffer_v.data());
                    for (InputIndex_ p = start, pe = start + length; p < pe; ++p, ++ptr) {
                        if (*ptr) {
                            store_v[p].push_back(*ptr);
                            store_i[p].push_back(x);
                        }
                    }
                }
            }, primary, options.num_threads);
        }
    }

    return output;
}

/**
 * @brief Options for `convert_to_fragmented_sparse()`.
 */
struct ConvertToFragmentedSparseOptions {
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
    bool row,
    const ConvertToFragmentedSparseOptions& options)
{
    auto frag = retrieve_fragmented_sparse_contents<StoredValue_, StoredIndex_>(
        matrix,
        row,
        [&]{
            RetrieveFragmentedSparseContentsOptions ropt;
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
