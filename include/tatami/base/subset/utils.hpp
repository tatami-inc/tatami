#ifndef TATAMI_DELAYED_SUBSET_UTILS_HPP
#define TATAMI_DELAYED_SUBSET_UTILS_HPP

#include "../Matrix.hpp"
#include "../utils.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>

namespace tatami {

namespace subset_utils {

struct DenseSupplement {
    std::vector<size_t> reverse_mapping;
    static constexpr bool sparse = false;
};

template<bool WORKROW, typename T, typename IDX, class InputWorkspace>
const T* remap_dense(const Matrix<T, IDX>* mat, size_t i, T* buffer, InputWorkspace* work, const std::vector<size_t>& rmapping) {
    const T* dump = extract_dense<WORKROW>(mat, i, work->buffer.data(), work->internal.get());
    auto temp = buffer;
    for (auto i : rmapping) {
        *temp = dump[i];
        ++temp;
    } 
    return buffer;
}

template<typename IDX>
struct SparseSupplement {
    std::vector<std::pair<size_t, size_t> > mapping_duplicates; 
    std::vector<IDX> mapping_duplicates_pool; 
    static constexpr bool sparse = true;
};

template<bool WORKROW, typename T, typename IDX, class InputWorkspace>
SparseRange<T, IDX> remap_sparse_duplicates(
    const Matrix<T, IDX>* mat, 
    size_t i, 
    T* vbuffer, 
    IDX* ibuffer, 
    InputWorkspace* work, 
    const std::vector<std::pair<size_t, size_t> >& dups, 
    const std::vector<IDX>& pool)
{
    // Allocation status of work->vbuffer depends on the extraction mode used to construct work->internal.
    T* vin = work->vbuffer.data();

    // work->ibuffer should always be allocated, as we need this to get the expanded counts.
    IDX* iin = work->ibuffer.data();
    auto raw = extract_sparse<WORKROW>(mat, i, vin, iin, work->internal.get());

    if (!raw.value) {
        vbuffer = NULL;
    }
    if (!(work->report_index)) {
        ibuffer = NULL;
    }

    auto vcopy = vbuffer;
    auto icopy = ibuffer;
    size_t counter = 0;

    for (size_t i = 0; i < raw.number; ++i) {
        const auto& pool_pos = dups[raw.index[i]];
        counter += pool_pos.second;

        if (vcopy) {
            std::fill(vcopy, vcopy + pool_pos.second, raw.value[i]);
            vcopy += pool_pos.second;
        }

        if (icopy) {
            auto istart = pool.begin() + pool_pos.first;
            std::copy(istart, istart + pool_pos.second, icopy);
            icopy += pool_pos.second;
        }
    }

    return SparseRange<T, IDX>(counter, vbuffer, ibuffer);
}

template<typename IDX, bool SPARSE>
using ConditionalSupplement = typename std::conditional<SPARSE, SparseSupplement<IDX>, DenseSupplement>::type;

}

}

#endif
