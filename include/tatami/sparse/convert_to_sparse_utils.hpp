#ifndef TATAMI_CONVERT_TO_SPARSE_UTILS_HPP
#define TATAMI_CONVERT_TO_SPARSE_UTILS_HPP

#include <vector>
#include <optional>

#include "../utils/consecutive_extractor.hpp"
#include "../utils/Index_to_container.hpp"

namespace tatami {

template<typename Index_, typename Count_>
struct CountNonZerosPerThread {
    std::vector<Index_> starts, lengths;
    std::vector<std::vector<Count_> > counts;
};

template<typename Value_, typename Index_, typename Count_>
std::optional<CountNonZerosPerThread<Index_, Count_> > count_sparse_non_zeros_inconsistent(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Index_ primary,
    const Index_ secondary,
    const bool row,
    Count_* const counts, // assume that this is already zeroed.
    const int threads
) {
    // First, we confirm that the counts don't overflow the Count_.
    sanisizer::cast<Count_>(secondary); 

    const bool is_sparse = matrix.is_sparse();
    const bool do_parallel = threads > 1;
    std::optional<std::vector<Index_> > all_partial_starts, all_partial_lengths;
    std::optional<std::vector<std::optional<std::vector<Count_> > > > all_partial_counts;
    if (do_parallel) {
        all_partial_starts.emplace(sanisizer::cast<I<decltype(all_partial_starts->size())> >(threads));
        all_partial_lengths.emplace(sanisizer::cast<I<decltype(all_partial_lengths->size())> >(threads));
        all_partial_counts.emplace(sanisizer::cast<I<decltype(all_partial_counts->size())> >(threads));
    }

    const int num_used = parallelize([&](const int thread, const Index_ start, const Index_ length) -> void {
        // To minimize false sharing, we allocate each buffer as a per-thread vector before moving it into the nnz_workers for serial use.
        // We skip the allocation for the first thread as this is allowed to use the (presumably zeroed) nnz array directly.
        Count_* cur_counts;
        std::optional<std::vector<Count_> > count_holder;
        if (!do_parallel) {
            cur_counts = counts;
        } else {
            count_holder.emplace(cast_Index_to_container_size<std::vector<Count_> >(primary));
            cur_counts = count_holder->data();
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
            for (Index_ x = 0; x < length; ++x) {
                const auto ptr = wrk->fetch(buffer_v.data());
                for (Index_ p = 0; p < primary; ++p) {
                    cur_counts[p] += (ptr[p] != 0);
                }
            }
        }

        if (do_parallel) {
            (*all_partial_counts)[thread] = std::move(count_holder);
            (*all_partial_starts)[thread] = start;
            (*all_partial_lengths)[thread] = length;
        }
    }, secondary, threads);

    std::optional<CountNonZerosPerThread<Index_, Count_> > output;
    if (do_parallel) {
        output.emplace();
        output->counts.reserve(num_used);
        output->starts = std::move(*all_partial_starts);
        output->lengths = std::move(*all_partial_lengths);
        for (int t = 0; t < num_used; ++t) {
            auto& y = *((*all_partial_counts)[t]);
            if (t == 0) {
                std::copy(y.begin(), y.end(), counts);
            } else {
                for (Index_ p = 0; p < primary; ++p) {
                    counts[p] += y[p];
                }
            }
            output->counts.push_back(std::move(y));
        }
    }

    return output;
}

template<typename Value_, typename Index_> 
std::vector<SparseRange<Value_, Index_> > extract_sparse_matrix(
    const tatami::Matrix<Value_, Index_>& matrix,
    std::vector<std::vector<Value_> >& store_v,
    std::vector<std::vector<Index_> >& store_i,
    const int num_threads
) {
    const bool row = matrix.prefer_rows();
    const Index_ NR = matrix.nrow();
    const Index_ NC = matrix.ncol();
    const Index_ primary = (row ? NR : NC);
    const Index_ secondary = (row ? NC : NR);

    resize_container_to_Index_size(store_v, primary);
    resize_container_to_Index_size(store_i, primary);
    auto original_ranges = create_container_of_Index_size<std::vector<SparseRange<Value_, Index_> > >(primary);

    if (matrix.is_sparse()) {
        parallelize([&](const int, const Index_ start, const Index_ length) -> void {
            auto wrk = consecutive_extractor<true>(matrix, row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<Value_> >(secondary);
            auto buffer_i = create_container_of_Index_size<std::vector<Index_> >(secondary);

            // Only make a copy of the underlying buffers if we don't have any choice.
            for (Index_ p = start, pend = start + length; p < pend; ++p) {
                auto range = wrk->fetch(buffer_v.data(), buffer_i.data());
                if (range.value == buffer_v.data()) {
                    auto& sv = store_v[p];
                    sv.insert(sv.end(), range.value, range.value + range.number);
                    range.value = sv.data();
                }
                if (range.index == buffer_i.data()) {
                    auto& si = store_i[p];
                    si.insert(si.end(), range.index, range.index + range.number);
                    range.index = si.data();
                }
                original_ranges[p] = std::move(range);
            }
        }, primary, num_threads);

    } else {
        parallelize([&](const int, const Index_ start, const Index_ length) -> void {
            auto wrk = consecutive_extractor<false>(matrix, row, start, length);
            auto buffer_v = create_container_of_Index_size<std::vector<Value_> >(secondary);

            for (Index_ p = start, pend = start + length; p < pend; ++p) {
                const auto ptr = wrk->fetch(buffer_v.data());
                auto& sv = store_v[p];
                auto& si = store_i[p];

                // For dense, we do treat zero values as structural zeros and remove them, otherwise the output wouldn't actually be sparse.
                Index_ nnz = 0;
                for (Index_ s = 0; s < secondary; ++s) {
                    nnz += (ptr[s] != 0);
                }
                sv.reserve(nnz);
                si.reserve(nnz);
                for (Index_ s = 0; s < secondary; ++s) {
                    const auto val = ptr[s];
                    if (val) {
                        sv.push_back(val);
                        si.push_back(s);
                    }
                }

                original_ranges[p] = SparseRange<Value_, Index_>(sv.size(), sv.data(), si.data());
            }
        }, primary, num_threads);
    }

    return original_ranges;
}

}

#endif
