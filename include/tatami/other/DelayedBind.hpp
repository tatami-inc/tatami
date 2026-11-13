#ifndef TATAMI_DELAYED_BIND_HPP
#define TATAMI_DELAYED_BIND_HPP

#include "../base/Matrix.hpp"
#include "../utils/new_extractor.hpp"
#include "../utils/ConsecutiveOracle.hpp"
#include "../utils/FixedOracle.hpp"
#include "../utils/PseudoOracularExtractor.hpp"
#include "../utils/copy.hpp"
#include "../utils/Index_to_container.hpp"

#include <numeric>
#include <algorithm>
#include <memory>
#include <array>
#include <type_traits>
#include <cstddef>
#include <limits>
#include <cassert>

/**
 * @file DelayedBind.hpp
 *
 * @brief Delayed combining of multiple `tatami::Matrix` objects.
 *
 * This is equivalent to the `DelayedAbind` class in the **DelayedArray** package.
 */

namespace tatami {

/**
 * @cond
 */
namespace DelayedBind_internal {

/*
 * Consider the problem of mapping an index 'i' of the combined dimension of the DelayedBind object to the component matrix.
 * The simplest approach is to do a binary search on 'cumulative' but that has O(log(n)) scaling for 'n' component matrices.
 * Another approach is to create a look-up table that can be indexed directly by 'i', which has O(1) runtime but consumes a lot of memory.
 *
 * The CompressedMapping class is a compromise that breaks up the combined dimensions into bins of regular size no greater than the extent of the smallest component matrix.
 * Each bin stores the identity the component matrix in which it starts, plus the length of that component matrix covered by that bin.
 * Given 'i', we can divide by the bin size to find its matching bin; we then compare the remainder to the covered length to assign it to the starting or next matrix.
 * This has O(1) runtime while being much more memory efficient for a minimum bin size >> 1.
 *
 * In practice, the bin size is set to the largest power of 2 that is no greater than the small matrix extent.
 * This allows for (i) some very fast division via bit shifting, and (ii) bitpacking to store both the matrix identity and the covered length in a single integer.
 */
template<typename Index_>
class CompressedMapping {
public:
    CompressedMapping() = default;

    CompressedMapping(Index_ min_size, const std::vector<Index_>& cumulative) {
        // Define bin size as the smallest power of 2 that is no greater than 'min_size'.
        // This means that 'my_shift' must be less than the number of binary digits of 'Index_'.
        while (min_size > 1) {
            ++my_shift;
            min_size /= 2;
        }
        assert(sanisizer::is_less_than(min_size, std::numeric_limits<Index_>::digits));
        const auto bin_size = static_cast<Index_>(1) << my_shift;
        my_mask = (bin_size - 1);

        const auto total_dim = cumulative.back();
        assert(min_size <= total_dim);
        const auto total_bins = (total_dim >> my_shift) + ((total_dim & my_mask) > 0);
        my_bins.reserve(total_bins);

        Index_ pos = 1;
        for (I<decltype(total_bins)> b = 0; b < total_bins; ++b) {
            // Find the first cumulative boundary that exceeds the start of the bin.
            // 'pos - 1' is the index of the matrix covering the start of the bin.
            Index_ lower = b << my_shift;
            while (cumulative[pos] <= lower) {
                ++pos;
            }

            // We compute the gap between the start of the bin and the cumulative boundary immediately to its right.
            // This is the number of elements in the bin that are assigned to the 'pos - 1'-th matrix.
            const auto until_next = sanisizer::min(cumulative[pos] - lower, bin_size);

            // Now we pack both of these values into a single integer to save some space.
            // This is possible based on the following logic:
            //
            // - The total cumulative dimension is less than '2^X' where 'X' is the number of digits of 'Index_'.
            // - The minimum size of each matrix is greater than or equal to '2^Y', where 'Y' is 'my_shift'.
            //   We note that 'X > Y' as the minimum size cannot be greater than the total cumulative dimension.
            // - Thus, the maximum possible value of the number of matrices must be less than '2^(X - Y)'.
            //   This implies that the matrix index must fit into a '(X - Y)'-bit unsigned integer.
            // - We know that 'until_next' must be less than or equal to '2^Y'.
            //   We also know it must be positive as the cumulative boundary exceeds the bin start.
            //   This implies that 'until_next - 1' must fit into a Y-bit unsigned integer.
            //
            // Thus, we can shift the matrix index by 'Y' and use the cleared bits to store 'until_next - 1'.
            const Index_ encoding = (static_cast<Index_>(pos - 1) << my_shift) | static_cast<Index_>(until_next - 1);
            my_bins.push_back(encoding);
        }
    }

    Index_ operator[](const Index_ i) const {
        const auto encoded = my_bins[i >> my_shift];
        // 'encoded >> my_shift' recovers the matrix index.
        // 'encoded & my_mask' recovers 'until_next - 1'.
        // We want to know if 'i & my_mask >= until_next', which is the same as 'i & my_mask > until_next - 1', saving us an addition.
        return sanisizer::sum_unsafe<Index_>(encoded >> my_shift, sanisizer::is_greater_than(i & my_mask, encoded & my_mask)); 
    }

public:
    std::vector<Index_> my_bins;
    Index_ my_shift = 0; 
    Index_ my_mask;
};

/**********************
 *** Dense parallel ***
 **********************/

template<typename Index_, class Initialize_>
Index_ initialize_parallel_block(
    const std::vector<Index_>& cumulative, 
    const CompressedMapping<Index_>& mapping,
    const Index_ block_start, 
    const Index_ block_length, 
    const Initialize_ init
) {
    const Index_ nmats = cumulative.size() - 1; // Number of matrices is guaranteed to fit in Index_, see reasoning in the DelayedBind constructor.
    if (nmats == 0) {
        return 0;
    }

    const Index_ start_index = mapping[block_start];
    Index_ actual_start = block_start - cumulative[start_index];
    const Index_ block_end = block_start + block_length;

    for (Index_ index = start_index; index < nmats; ++index) {
        const Index_ submat_end = cumulative[index + 1]; 
        bool not_final = (block_end > submat_end);
        const Index_ actual_end = (not_final ? submat_end : block_end) - cumulative[index];
        init(index, actual_start, actual_end - actual_start);
        if (!not_final) {
            break;
        }
        actual_start = 0;
    }

    return start_index;
}

template<typename Index_, class Initialize_>
void initialize_parallel_index(
    const std::vector<Index_>& cumulative, 
    const CompressedMapping<Index_>& mapping,
    const std::vector<Index_>& indices, 
    const Initialize_ init
) {
    Index_ counter = 0;
    const Index_ il = indices.size();
    while (counter < il) {
        const Index_ first_index = indices[counter];
        const Index_ bind_index = mapping[first_index];
        const Index_ lower = cumulative[bind_index];
        const Index_ upper = cumulative[bind_index + 1];

        // Creating the slice with one element already.
        auto slice_ptr = std::make_shared<std::vector<Index_> >(1, first_index - lower);
        ++counter;

        while (counter < il && indices[counter]  < upper) {
            slice_ptr->push_back(indices[counter] - lower);
            ++counter;
        }

        init(bind_index, std::move(slice_ptr));
    }
}

template<bool oracle_, typename Value_, typename Index_>
class  ParallelDense final : public DenseExtractor<oracle_, Value_, Index_> {
public:
    ParallelDense(
        const std::vector<Index_>&, // Not used, just provided for consistency with other constructors.
        const CompressedMapping<Index_>&, // Ditto.
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) {
        my_exts.reserve(matrices.size());
        my_count.reserve(matrices.size());
        for (const auto& m : matrices) {
            my_count.emplace_back(row ? m->ncol() : m->nrow());
            my_exts.emplace_back(new_extractor<false, oracle_>(m.get(), row, oracle, opt));
        }
    }

    ParallelDense(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping,
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start, 
        const Index_ block_length, 
        const Options& opt
    ) {
        my_exts.reserve(matrices.size());
        my_count.reserve(matrices.size());
        initialize_parallel_block(
            cumulative, 
            mapping,
            block_start, 
            block_length,
            [&](const Index_ i, const Index_ sub_block_start, const Index_ sub_block_length) -> void {
                my_count.emplace_back(sub_block_length);
                my_exts.emplace_back(new_extractor<false, oracle_>(matrices[i].get(), row, oracle, sub_block_start, sub_block_length, opt));
            }
        );
    }

    ParallelDense(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping,
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) {
        my_exts.reserve(matrices.size());
        my_count.reserve(matrices.size());
        initialize_parallel_index(
            cumulative, 
            mapping,
            *indices_ptr,
            [&](const Index_ i, VectorPtr<Index_> sub_indices_ptr) -> void {
                my_count.emplace_back(sub_indices_ptr->size());
                my_exts.emplace_back(new_extractor<false, oracle_>(matrices[i].get(), row, oracle, std::move(sub_indices_ptr), opt));
            }
        );
    }

public:
    const Value_* fetch(const Index_ i, Value_* const buffer) {
        auto copy = buffer;
        const Index_ nmats = my_count.size();
        for (Index_ x = 0; x < nmats; ++x) {
            const auto ptr = my_exts[x]->fetch(i, copy); 
            const auto num = my_count[x];
            copy_n(ptr, num, copy);
            copy += num;
        }
        return buffer;
    }

private:
    std::vector<std::unique_ptr<DenseExtractor<oracle_, Value_, Index_> > > my_exts;
    std::vector<Index_> my_count;
};

/***********************
 *** Sparse parallel ***
 ***********************/

template<bool oracle_, typename Value_, typename Index_>
class ParallelFullSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    ParallelFullSparse(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>&, // not actually used, just provided for consistency with the other constructors.
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Options& opt
    ) : 
        my_cumulative(cumulative),
        my_needs_value(opt.sparse_extract_value),
        my_needs_index(opt.sparse_extract_index)
    {
        my_exts.reserve(matrices.size());
        for (const auto& m : matrices) {
            my_exts.emplace_back(new_extractor<true, oracle_>(m.get(), row, oracle, opt));
        }
    }

    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        auto vcopy = value_buffer;
        auto icopy = index_buffer;
        Index_ accumulated = 0;

        const Index_ nmats = my_exts.size();
        for (Index_ x = 0; x < nmats; ++x) {
            const auto range = my_exts[x]->fetch(i, vcopy, icopy); 
            accumulated += range.number;
            if (my_needs_value) {
                copy_n(range.value, range.number, vcopy);
                vcopy += range.number;
            }
            if (my_needs_index) {
                const auto offset = my_cumulative[x];
                for (Index_ y = 0; y < range.number; ++y) {
                    icopy[y] = range.index[y] + offset;
                }
                icopy += range.number;
            }
        }

        return SparseRange<Value_, Index_>(accumulated, (my_needs_value ? value_buffer : NULL), (my_needs_index ? index_buffer : NULL));
    }

private:
    const std::vector<Index_>& my_cumulative;
    bool my_needs_value, my_needs_index;
    std::vector<std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > > my_exts;
};

template<bool oracle_, typename Value_, typename Index_>
class ParallelBlockSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    ParallelBlockSparse(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping,
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        const Index_ block_start, 
        const Index_ block_length, 
        const Options& opt
    ) : 
        my_cumulative(cumulative), 
        my_needs_value(opt.sparse_extract_value), 
        my_needs_index(opt.sparse_extract_index) 
    {
        my_exts.reserve(matrices.size());
        my_start_matrix = initialize_parallel_block(
            my_cumulative, 
            mapping,
            block_start, 
            block_length,
            [&](const Index_ i, const Index_ sub_block_start, const Index_ sub_block_length) -> void {
                my_exts.emplace_back(new_extractor<true, oracle_>(matrices[i].get(), row, oracle, sub_block_start, sub_block_length, opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        auto vcopy = value_buffer;
        auto icopy = index_buffer;
        Index_ count = 0;

        const Index_ nmats = my_exts.size();
        for (Index_ x = 0; x < nmats; ++x) {
            const auto range = my_exts[x]->fetch(i, vcopy, icopy);
            count += range.number;
            if (my_needs_value) {
                copy_n(range.value, range.number, vcopy);
                vcopy += range.number;
            }
            if (my_needs_index) {
                const Index_ offset = my_cumulative[x + my_start_matrix];
                for (Index_ y = 0; y < range.number; ++y) {
                    icopy[y] = range.index[y] + offset;
                }
                icopy += range.number;
            }
        }

        return SparseRange<Value_, Index_>(count, (my_needs_value ? value_buffer : NULL), (my_needs_index ? index_buffer : NULL));
    }

private:
    const std::vector<Index_>& my_cumulative;
    bool my_needs_value, my_needs_index;
    std::vector<std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > > my_exts;
    Index_ my_start_matrix;
};

template<bool oracle_, typename Value_, typename Index_>
class ParallelIndexSparse final : public SparseExtractor<oracle_, Value_, Index_> {
public:
    ParallelIndexSparse(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping,
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        MaybeOracle<oracle_, Index_> oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) : 
        my_cumulative(cumulative),
        my_needs_value(opt.sparse_extract_value), 
        my_needs_index(opt.sparse_extract_index) 
    {
        my_exts.reserve(matrices.size());
        my_which_matrix.reserve(matrices.size());
        initialize_parallel_index(
            my_cumulative, 
            mapping,
            *indices_ptr,
            [&](const Index_ i, VectorPtr<Index_> sub_indices_ptr) -> void {
                my_which_matrix.emplace_back(i);
                my_exts.emplace_back(new_extractor<true, oracle_>(matrices[i].get(), row, oracle, std::move(sub_indices_ptr), opt));
            }
        );
    }

    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const value_buffer, Index_* const index_buffer) {
        auto vcopy = value_buffer;
        auto icopy = index_buffer;
        Index_ count = 0;

        const Index_ nmats = my_which_matrix.size();
        for (Index_ x = 0; x < nmats; ++x) {
            const auto range = my_exts[x]->fetch(i, vcopy, icopy);
            count += range.number;
            if (my_needs_value) {
                copy_n(range.value, range.number, vcopy);
                vcopy += range.number;
            }

            if (my_needs_index) {
                const Index_ offset = my_cumulative[my_which_matrix[x]];
                for (Index_ y = 0; y < range.number; ++y) {
                    icopy[y] = range.index[y] + offset;
                }
                icopy += range.number;
            }
        }

        return SparseRange<Value_, Index_>(count, (my_needs_value ? value_buffer : NULL), (my_needs_index ? index_buffer : NULL));
    }

private:
    const std::vector<Index_>& my_cumulative;
    bool my_needs_value, my_needs_index;
    std::vector<std::unique_ptr<SparseExtractor<oracle_, Value_, Index_> > > my_exts;
    std::vector<Index_> my_which_matrix;
};

/****************************
 *** Myopic Perpendicular ***
 ****************************/

template<typename Value_, typename Index_>
class MyopicPerpendicularDense final : public MyopicDenseExtractor<Value_, Index_> {
public:
    template<typename ... Args_>
    MyopicPerpendicularDense(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row, 
        const Args_& ... args
    ) : 
        my_cumulative(cumulative),
        my_mapping(mapping)
    {
        my_exts.reserve(matrices.size());
        for (const auto& m : matrices) {
            my_exts.emplace_back(m->dense(row, args...));
        }
    }

    const Value_* fetch(const Index_ i, Value_* const buffer) {
        const Index_ chosen = my_mapping[i];
        return my_exts[chosen]->fetch(i - my_cumulative[chosen], buffer);
    }

private:
    const std::vector<Index_>& my_cumulative;
    const CompressedMapping<Index_>& my_mapping;
    std::vector<std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > > my_exts;
};

template<typename Value_, typename Index_>
class MyopicPerpendicularSparse final : public MyopicSparseExtractor<Value_, Index_> {
public:
    template<typename ... Args_>
    MyopicPerpendicularSparse(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        const Args_& ... args
    ) : 
        my_cumulative(cumulative),
        my_mapping(mapping)
    {
        my_exts.reserve(matrices.size());
        for (const auto& m : matrices) {
            my_exts.emplace_back(m->sparse(row, args...));
        }
    }

    SparseRange<Value_, Index_> fetch(const Index_ i, Value_* const vbuffer, Index_* const ibuffer) {
        const Index_ chosen = my_mapping[i];
        return my_exts[chosen]->fetch(i - my_cumulative[chosen], vbuffer, ibuffer);
    }

private:
    const std::vector<Index_>& my_cumulative;
    const CompressedMapping<Index_>& my_mapping;
    std::vector<std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > > my_exts;
};

/******************************
 *** Oracular Perpendicular ***
 ******************************/

template<typename Index_>
class OffsetOracleWrapper : public Oracle<Index_> {
public:
    OffsetOracleWrapper(std::shared_ptr<const Oracle<Index_> > ptr, PredictionIndex start, PredictionIndex number, Index_ offset) :
        my_ptr(std::move(ptr)),
        my_start(start),
        my_number(number),
        my_offset(offset)
    {}

    PredictionIndex total() const {
        return my_number;
    }

    Index_ get(PredictionIndex i) const {
        return my_ptr->get(i + my_start) - my_offset;
    }
private:
    std::shared_ptr<const Oracle<Index_> > my_ptr;
    PredictionIndex my_start, my_number;
    Index_ my_offset;
};

template<typename Index_>
class OffsetOracleWrapperWithExtras : public Oracle<Index_> {
public:
    OffsetOracleWrapperWithExtras(std::shared_ptr<const Oracle<Index_> > ptr, PredictionIndex start, PredictionIndex number, Index_ offset, std::vector<Index_> extras) :
        my_ptr(std::move(ptr)),
        my_start(start),
        my_number(number),
        my_offset(offset),
        my_extras(std::move(extras))
    {
        assert(sanisizer::sum<PredictionIndex>(my_number, my_extras.size()));
    }

    PredictionIndex total() const {
        // This addition should be safe under the assumption that this is just a subset of predictions from the parent oracle.
        return my_number + my_extras.size();
    }

    Index_ get(PredictionIndex i) const {
        if (sanisizer::is_greater_than_or_equal(i, my_number)) {
            return my_extras[i - my_number];
        } else {
            return my_ptr->get(i + my_start) - my_offset;
        }
    }
private:
    std::shared_ptr<const Oracle<Index_> > my_ptr;
    PredictionIndex my_start, my_number;
    Index_ my_offset;
    std::vector<Index_> my_extras;
};

template<typename Index_>
class ConsecutiveOracleWithExtras : public Oracle<Index_> {
public:
    ConsecutiveOracleWithExtras(Index_ start, Index_ number, std::vector<Index_> extras) :
        my_start(start),
        my_number(number),
        my_extras(std::move(extras))
    {
        assert(sanisizer::sum<PredictionIndex>(my_number, my_extras.size()));
    }

    PredictionIndex total() const {
        // This addition should be safe under the assumption that this is just a subset of predictions from the parent oracle.
        return my_number + my_extras.size();
    }

    Index_ get(PredictionIndex i) const {
        if (sanisizer::is_greater_than_or_equal(i, my_number)) {
            return my_extras[i - my_number];
        } else {
            return my_start + i;
        }
    }
private:
    Index_ my_start, my_number;
    std::vector<Index_> my_extras;
};

template<typename Index_>
struct PredictionProperties {
    bool consecutive_predictions = true;
    Index_ start_prediction = 0;
    Index_ number_of_predictions = 0;

    bool consecutive_calls = true;
    PredictionIndex start_call = 0;
    PredictionIndex number_of_calls = 0;

    std::vector<Index_> realized_predictions;
};

template<typename Index_>
std::vector<PredictionProperties<Index_> > compute_prediction_properties(
    const std::vector<Index_>& cumulative,
    const CompressedMapping<Index_>& mapping, 
    const Oracle<Index_>& oracle
) {
    const auto nmats = cumulative.size() - 1;
    auto predictions = create_container_of_Index_size<std::vector<PredictionProperties<Index_> > >(nmats); // nmats should fit in an Index_, so this call is legal.

    const auto ntotal = oracle.total();
    for (I<decltype(ntotal)> i = 0; i < ntotal; ++i) {
        const auto raw_prediction = oracle.get(i);
        const Index_ choice = mapping[raw_prediction];
        const auto pred = raw_prediction - cumulative[choice];

        auto& current = predictions[choice];
        if (current.consecutive_predictions) {
            if (current.number_of_predictions == 0) {
                current.start_prediction = pred;
                current.number_of_predictions = 1;
            } else if (current.number_of_predictions + current.start_prediction == pred) {
                ++current.number_of_predictions;
            } else {
                current.consecutive_predictions = false;
            }
        }

        if (current.consecutive_calls) {
            if (current.number_of_calls == 0) {
                current.start_call = i;
                current.number_of_calls = 1;
            } else if (current.number_of_calls + current.start_call == i) {
                ++current.number_of_calls;
            } else {
                current.consecutive_calls = false;
            }
        }

        // We only need to start recording the realized predictions when both runs have failed,
        // as we will always use the longer run to initialize an Oracle.
        if (!current.consecutive_predictions && !current.consecutive_calls) {
            current.realized_predictions.push_back(pred);
        }
    }

    return predictions;
}

template<typename Index_, class Initialize_>
void initialize_perp_oracular(
    const std::vector<Index_>& cumulative, 
    const std::shared_ptr<const Oracle<Index_> >& ora, 
    std::vector<PredictionProperties<Index_> > properties,
    Initialize_ init
) {
    const auto nmats = cumulative.size() - 1;
    assert(sanisizer::is_equal(nmats, properties.size()));

    for (I<decltype(nmats)> x = 0; x < nmats; ++x) {
        auto& prop = properties[x];
        if (prop.consecutive_predictions) {
            if (prop.number_of_predictions) {
                init(x, std::make_shared<ConsecutiveOracle<Index_> >(prop.start_prediction, prop.number_of_predictions));
            }
        } else if (prop.consecutive_calls) {
            if (prop.number_of_calls) {
                init(x, std::make_shared<OffsetOracleWrapper<Index_> >(ora, prop.start_call, prop.number_of_calls, cumulative[x]));
            }
        } else if (sanisizer::is_greater_than(prop.number_of_predictions, prop.number_of_calls)) {
            init(x, std::make_shared<ConsecutiveOracleWithExtras<Index_> >(prop.start_prediction, prop.number_of_predictions, std::move(prop.realized_predictions)));
        } else if (sanisizer::is_less_than(prop.number_of_predictions, prop.number_of_calls)) {
            init(x, std::make_shared<OffsetOracleWrapperWithExtras<Index_> >(ora, prop.start_call, prop.number_of_calls, cumulative[x], std::move(prop.realized_predictions)));
        } else if (prop.number_of_predictions > 1) {
            // i.e., if the runs are equal and are not of length 1, we prefer the ConsecutiveOracleWithExtras as that doesn't need to copy 'ora'.
            init(x, std::make_shared<ConsecutiveOracleWithExtras<Index_> >(prop.start_prediction, prop.number_of_predictions, std::move(prop.realized_predictions)));
        } else {
            // Need to add back the first element from a length-1 run.
            prop.realized_predictions.insert(prop.realized_predictions.begin(), 1, prop.start_prediction);
            init(x, std::make_shared<FixedVectorOracle<Index_> >(std::move(prop.realized_predictions)));
        }
    }
}

template<typename Value_, typename Index_>
class OracularPerpendicularDense final : public OracularDenseExtractor<Value_, Index_> {
public:
    template<typename ... Args_>
    OracularPerpendicularDense(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        std::shared_ptr<const Oracle<Index_> > ora, 
        const Args_& ... args
    ) :
        my_mapping(mapping),
        my_oracle(std::move(ora))
    {
        resize_container_to_Index_size(my_exts, matrices.size()); // number of matrices should fit in an Index_, so this call is allowed.
        initialize_perp_oracular(
            cumulative,
            my_oracle,
            compute_prediction_properties(cumulative, mapping, *my_oracle),
            [&](const Index_ x, std::shared_ptr<const Oracle<Index_> > subora) -> void {
                my_exts[x] = matrices[x]->dense(row, std::move(subora), args...);
            }
        );
    }

    const Value_* fetch(const Index_ i, Value_* const buffer) {
        const auto chosen = my_mapping[my_oracle->get(my_used++)];
        const auto output = my_exts[chosen]->fetch(i, buffer);
        return output;
    }

private:
    const CompressedMapping<Index_>& my_mapping;
    std::shared_ptr<const Oracle<Index_> > my_oracle;
    PredictionIndex my_used = 0;
    std::vector<std::unique_ptr<OracularDenseExtractor<Value_, Index_> > > my_exts;
};

template<typename Value_, typename Index_>
class OracularPerpendicularSparse final : public OracularSparseExtractor<Value_, Index_> {
public:
    template<typename ... Args_>
    OracularPerpendicularSparse(
        const std::vector<Index_>& cumulative, 
        const CompressedMapping<Index_>& mapping, 
        const std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >& matrices, 
        const bool row,
        std::shared_ptr<const Oracle<Index_> > ora, 
        const Args_& ... args
    ) :
        my_mapping(mapping),
        my_oracle(std::move(ora))
    {
        resize_container_to_Index_size(my_exts, matrices.size()); // number of matrices should fit in an Index_, so this call is legal.
        initialize_perp_oracular(
            cumulative,
            my_oracle,
            compute_prediction_properties(cumulative, mapping, *my_oracle),
            [&](const Index_ x, std::shared_ptr<const Oracle<Index_> > subora) -> void {
                my_exts[x] = matrices[x]->sparse(row, std::move(subora), args...);
            }
        );
    }

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        const auto chosen = my_mapping[my_oracle->get(my_used++)];
        const auto output = my_exts[chosen]->fetch(i, vbuffer, ibuffer);
        return output;
    }

private:
    const CompressedMapping<Index_>& my_mapping;
    std::shared_ptr<const Oracle<Index_> > my_oracle;
    PredictionIndex my_used = 0;
    std::vector<std::unique_ptr<OracularSparseExtractor<Value_, Index_> > > my_exts;
};

}
/**
 * @endcond
 */

/**
 * @brief Delayed combining of a matrix.
 *
 * Implements delayed combining by rows or columns of a matrix.
 * This operation is "delayed" in that it is only performed on a row-by-row or column-by-column basis during data extraction.
 *
 * @tparam Value_ Type of matrix value.
 * @tparam Index_ Type of index value.
 */
template<typename Value_, typename Index_>
class DelayedBind final : public Matrix<Value_, Index_> {
public:
    /**
     * @param matrices Pointers to the matrices to be combined.
     * All matrices to be combined should have the same number of columns (if `row = true`) or rows (otherwise).
     * @param by_row Whether to combine matrices by the rows (i.e., the output matrix has number of rows equal to the sum of the number of rows in `matrices`).
     * If false, combining is applied by the columns.
     */
    DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > matrices, const bool by_row) : 
        my_matrices(std::move(matrices)), my_by_row(by_row)
    {
        auto nmats = my_matrices.size();
        my_cumulative.reserve(sanisizer::sum<I<decltype(my_cumulative.size())> >(nmats, 1));
        I<decltype(nmats)> sofar = 0;
        my_cumulative.push_back(0);
        Index_ min_extent = std::numeric_limits<Index_>::max();

        for (I<decltype(nmats)> i = 0; i < nmats; ++i) {
            auto& current = my_matrices[i];
            Index_ primary, secondary;
            if (my_by_row) {
                primary = current->nrow();
                secondary = current->ncol();
            } else {
                primary = current->ncol();
                secondary = current->nrow();
            }

            if (i == 0) {
                my_otherdim = secondary;
            } else if (my_otherdim != secondary) {
                throw std::runtime_error("all 'matrices' should have the same number of " + (my_by_row ? std::string("columns") : std::string("rows")));
            }

            // Removing the matrices that don't contribute anything,
            // so we don't have to deal with their overhead.
            if (primary > 0) {
                if (sofar != i) {
                    my_matrices[sofar] = std::move(current);
                }
                my_cumulative.push_back(sanisizer::sum<Index_>(attest_for_Index(my_cumulative.back()), attest_for_Index(primary)));
                min_extent = std::min(min_extent, primary);
                ++sofar;
            }
        }

        my_matrices.resize(sofar);
        nmats = sofar;

        if (nmats == 0) {
            min_extent = 0;
        }
        my_mapping = DelayedBind_internal::CompressedMapping(min_extent, my_cumulative);

        double denom = 0;
        for (const auto& x : my_matrices) {
            const double total = static_cast<double>(x->nrow()) * static_cast<double>(x->ncol());
            denom += total;
            my_sparse_prop += total * x->is_sparse_proportion();
            my_by_row_prop += total * x->prefer_rows_proportion();
        }
        if (denom) {
            my_sparse_prop /= denom;
            my_by_row_prop /= denom;
        }

        for (int d = 0; d < 2; ++d) {
            my_uses_oracle[d] = false;
            for (const auto& x : my_matrices) {
                if (x->uses_oracle(d)) {
                    my_uses_oracle[d] = true;
                    break;
                }
            }
        }
    }

    /**
     * @cond
     */
    // Soft-deprecated.
    DelayedBind(const std::vector<std::shared_ptr<Matrix<Value_, Index_> > >& matrices, const bool by_row) : 
        DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > >(matrices.begin(), matrices.end()), by_row) {}
    /**
     * @endcond
     */

private:
    std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > my_matrices;
    bool my_by_row;

    Index_ my_otherdim = 0;
    std::vector<Index_> my_cumulative;
    DelayedBind_internal::CompressedMapping<Index_> my_mapping;

    double my_sparse_prop = 0, my_by_row_prop = 0;
    std::array<bool, 2> my_uses_oracle;

public:
    Index_ nrow() const {
        if (my_by_row) {
            return my_cumulative.back();
        } else {
            return my_otherdim;
        }
    }

    Index_ ncol() const {
        if (my_by_row) {
            return my_otherdim;
        } else {
            return my_cumulative.back();
        }
    }

    bool is_sparse() const {
        return my_sparse_prop > 0.5;
    }

    double is_sparse_proportion() const {
        return my_sparse_prop;
    }

    bool prefer_rows() const {
        return my_by_row_prop > 0.5;
    }

    double prefer_rows_proportion() const {
        return my_by_row_prop;
    }

    bool uses_oracle(const bool row) const {
        return my_uses_oracle[row];
    }

    using Matrix<Value_, Index_>::dense;

    using Matrix<Value_, Index_>::sparse;

    /**********************************
     ********** Myopic dense **********
     **********************************/
public:
    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->dense(row, opt);
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<false, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                false,
                opt
            );
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->dense(row, block_start, block_length, opt);
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                block_start,
                block_length,
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<false, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                false,
                block_start,
                block_length,
                opt
            );
        }
    }

    std::unique_ptr<MyopicDenseExtractor<Value_, Index_> > dense(
        const bool row,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->dense(row, std::move(indices_ptr), opt);
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularDense<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(indices_ptr),
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<false, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                false,
                std::move(indices_ptr),
                opt
            );
        }
    }

    /***********************************
     ********** Myopic sparse **********
     ***********************************/
private:
    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->sparse(row, opt);
        } else  if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelFullSparse<false, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                false,
                opt
            );
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->sparse(row, block_start, block_length, opt);
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                block_start,
                block_length,
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelBlockSparse<false, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                false,
                block_start,
                block_length,
                opt
            );
        }
    }

    std::unique_ptr<MyopicSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->sparse(row, std::move(indices_ptr), opt);
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::MyopicPerpendicularSparse<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(indices_ptr),
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelIndexSparse<false, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                false,
                std::move(indices_ptr),
                opt
            );
        }
    }

    /************************************
     ********** Oracular dense **********
     ************************************/
public:
    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->dense(row, std::move(oracle), opt);
        } else if (!my_uses_oracle[row]) {
            return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(
                std::move(oracle),
                dense(row, opt)
            );
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<true, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                opt
            );
        }
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->dense(row, std::move(oracle), block_start, block_length, opt);
        } else if (!my_uses_oracle[row]) {
            return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(
                std::move(oracle),
                dense(row, block_start, block_length, opt)
            );
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                block_start,
                block_length,
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<true, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                block_start,
                block_length,
                opt
            );
        }
    }

    std::unique_ptr<OracularDenseExtractor<Value_, Index_> > dense(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->dense(row, std::move(oracle), std::move(indices_ptr), opt);
        } else if (!my_uses_oracle[row]) {
            return std::make_unique<PseudoOracularDenseExtractor<Value_, Index_> >(
                std::move(oracle),
                dense(row, std::move(indices_ptr), opt)
            );
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularDense<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                std::move(indices_ptr),
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelDense<true, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                std::move(indices_ptr),
                opt
            );
        }
    }

    /*************************************
     ********** Oracular sparse **********
     *************************************/
private:
    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->sparse(row, std::move(oracle), opt);
        } else if (!my_uses_oracle[row]) {
            return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(
                std::move(oracle),
                sparse(row, opt)
            );
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelFullSparse<true, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                opt
            );
        }
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        const Index_ block_start,
        const Index_ block_length,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->sparse(row, std::move(oracle), block_start, block_length, opt);
        } else if (!my_uses_oracle[row]) {
            return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(
                std::move(oracle),
                sparse(row, block_start, block_length, opt)
            );
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                block_start,
                block_length,
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelBlockSparse<true, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                block_start,
                block_length,
                opt
            );
        }
    }

    std::unique_ptr<OracularSparseExtractor<Value_, Index_> > sparse(
        const bool row,
        std::shared_ptr<const Oracle<Index_> > oracle,
        VectorPtr<Index_> indices_ptr,
        const Options& opt
    ) const {
        if (my_matrices.size() == 1) {
            return my_matrices[0]->sparse(row, std::move(oracle), std::move(indices_ptr), opt);
        } else if (!my_uses_oracle[row]) {
            return std::make_unique<PseudoOracularSparseExtractor<Value_, Index_> >(
                std::move(oracle),
                sparse(row, std::move(indices_ptr), opt)
            );
        } else if (row == my_by_row) {
            return std::make_unique<DelayedBind_internal::OracularPerpendicularSparse<Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                std::move(indices_ptr),
                opt
            );
        } else {
            return std::make_unique<DelayedBind_internal::ParallelIndexSparse<true, Value_, Index_> >(
                my_cumulative,
                my_mapping,
                my_matrices,
                row,
                std::move(oracle),
                std::move(indices_ptr),
                opt
            );
        }
    }
};

/**
 * @cond
 */
// These methods are soft-deprecated: kept around for back-compatibility only.
template<typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > matrices, bool row) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBind<Value_, Index_>(std::move(matrices), row));
}

template<typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<Matrix<Value_, Index_> > > matrices, bool row) {
    return std::shared_ptr<Matrix<Value_, Index_> >(new DelayedBind<Value_, Index_>(std::move(matrices), row));
}

template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<const Matrix<Value_, Index_> > > matrices) {
    return make_DelayedBind(std::move(matrices), margin_ == 0);
}

template<int margin_, typename Value_, typename Index_>
std::shared_ptr<Matrix<Value_, Index_> > make_DelayedBind(std::vector<std::shared_ptr<Matrix<Value_, Index_> > > matrices) {
    return make_DelayedBind(std::move(matrices), margin_ == 0);
}
/**
 * @endcond
 */

}

#endif
