#ifndef TATAMI_PROCESS_CONSECUTIVE_INDICES_HPP
#define TATAMI_PROCESS_CONSECUTIVE_INDICES_HPP

/**
 * @file process_consecutive_indices.hpp
 * @brief Utility to process consecutive indices.
 */

namespace tatami {

/**
 * Process runs of consecutive indices that are used in the index-aware `dense_row()`, `sparse_column()`, etc. methods.
 * This provides some opportunities for optimization when the indices contain contiguous stretches.
 * For example, third-party libraries can be asked to process blocks of observations rather than handling them one at a time.
 *
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Function_ Function to apply to each contiguous run.
 *
 * @param indices Pointer to an array of sorted and unique indices for row/column elements.
 * @param length Length of the array pointed to by `indices`.
 * @param fun Function to apply to each contiguous run of indices.
 * This should take two arguments - the start index of each run, and the length of the run.
 * Calls to `fun` are guaranteed to contain increasing start indices with non-overlapping runs.
 * The return value of this function is ignored.
 */
template<typename Index_, class Function_>
void process_consecutive_indices(const Index_* const indices, const Index_ length, const Function_ fun) {
    if (length > 0) {
        Index_ start = indices[0], last = start + 1;
        for (Index_ i = 1; i < length; ++i) {
            if (indices[i] > last) {
                fun(start, last - start);
                start = indices[i];
            }
            last = indices[i] + 1;
        }
        fun(start, last - start);
    }
}

}

#endif
