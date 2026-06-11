#ifndef TATAMI_CONVERT_TO_SPARSE_UTILS_HPP
#define TATAMI_CONVERT_TO_SPARSE_UTILS_HPP

#include <vector>
#include <optional>

#include "../utils/consecutive_extractor.hpp"
#include "../utils/Index_to_container.hpp"

namespace tatami {

template<typename Value_, typename Index_, typename Count_>
void count_sparse_non_zeros_inconsistent(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Index_ primary,
    const Index_ secondary,
    const bool row,
    Count_* const nnz_inconsistent,
    std::optional<std::vector<Index_> >& nnz_consistent,
    const int threads
) {
    // First, we confirm that the counts don't overflow the Count_.
    sanisizer::cast<Count_>(secondary); 

    const bool do_parallel = threads > 1;
    std::optional<std::vector<std::optional<std::vector<Count_> > > > all_partial_counts;
    if (do_parallel) {
        all_partial_counts.emplace(sanisizer::cast<I<decltype(all_partial_counts->size())> >(threads - 1));
    }

    const bool is_sparse = matrix.is_sparse();
    if (!is_sparse) {
        if (do_parallel) {
            // We only report the secondary non-zero counts in the dense case, just so that downstream functions can easily preallocate memory.
            // This is not necessary in the sparse case as we get the number of non-zeros directly from the query.
            nnz_consistent.emplace(cast_Index_to_container_size<std::vector<Index_> >(secondary));
        }
    }

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
            }

        } else {
            auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<Value_> >(primary);

            if (do_parallel) {
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

            } else {
                for (Index_ x = 0; x < length; ++x) {
                    const auto ptr = wrk->fetch(buffer_v.data());
                    for (Index_ p = 0; p < primary; ++p) {
                        cur_counts[p] += (ptr[p] != 0);
                    }
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

// Here, the general strategy is to store all non-zero elements in thread-specific containers, and then stitch them together into the output buffers in the serial section.
// This avoids false sharing at the cost of almost doubling the memory usage as the number of threads increases - I think this is mostly acceptable.
// (There is no penalty for single-threaded use as we just use the output buffers directly.) 
//
// Contrast this to the previous approach, where we take a different slice of each primary dimension element in each thread.
// This tends to be a suboptimal extraction strategy as each thread would need to do some indexing to find the slice.
// Depending on the representation, the work done in each thread might be significant, e.g., a full read of the indices from file in a HDF5-backed sparse matrix.
template<typename InputValue_, typename InputIndex_, class SparseMain_, class DenseMain_, class Reduce_>
void fill_sparse_matrix_inconsistent(
    const tatami::Matrix<InputValue_, InputIndex_>& matrix,
    const InputIndex_ primary,
    const InputIndex_ secondary,
    const bool row,
    const std::optional<std::vector<InputIndex_> > & nnz_consistent, // count of the number of non-zero elements in each secondary dimension element, for buffer allocations.
    SparseMain_ sparse_main,
    DenseMain_ dense_main,
    Reduce_ reduce,
    const int threads
) {
    const bool do_parallel = threads > 0;
    InputIndex_ last_secondary_added = 0;
    std::optional<std::vector<std::vector<InputValue_> > > all_partial_values;
    std::optional<std::vector<std::vector<InputIndex_> > > all_partial_primary_indices;
    if (do_parallel) {
        all_partial_values.emplace(sanisizer::cast<I<decltype(all_partial_values->size())> >(secondary));
        all_partial_primary_indices.emplace(sanisizer::cast<I<decltype(all_partial_primary_indices->size())> >(secondary));
    }

    const bool is_sparse = matrix.is_sparse();
    parallelize([&](const int, const InputIndex_ start, const InputIndex_ length) -> void {
        if (!do_parallel || start == 0) {
            last_secondary_added = start + length;

            if (is_sparse){ 
                Options opt;
                opt.sparse_ordered_index = false;
                auto wrk = consecutive_extractor<true>(matrix, !row, start, length, opt);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                auto buffer_i = create_container_of_Index_size<std::vector<InputIndex_> >(primary);
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                    sparse_main(x, range); // start == 0, so no need to add start to get the actual secondary dimension index.
                }

            } else {
                auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);
                for (InputIndex_ x = 0; x < length; ++x) {
                    const auto ptr = wrk->fetch(buffer_v.data());
                    dense_main(x, ptr); // again, start == 0.
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
                    (*all_partial_values)[start + x] = std::vector<InputValue_>(range.value, range.value + range.number);
                    (*all_partial_primary_indices)[start + x] = std::vector<InputIndex_>(range.index, range.index + range.number);
                }

            } else {
                auto wrk = consecutive_extractor<false>(matrix, !row, start, length);
                auto buffer_v = create_container_of_Index_size<std::vector<InputValue_> >(primary);

                for (InputIndex_ x = 0; x < length; ++x) {
                    // Note, nnz_consistent is only populated if !is_sparse && do_parallel.
                    const auto nnz = (*nnz_consistent)[start + x]; 
                    std::vector<InputValue_> out_v;
                    std::vector<InputIndex_> out_i;
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
            const auto& cur_values = (*all_partial_values)[s];
            const auto& cur_primary_indices = (*all_partial_primary_indices)[s];
            reduce(s, cur_values, cur_primary_indices);
        }
    }
}

}

#endif
