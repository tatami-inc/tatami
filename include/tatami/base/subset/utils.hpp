#ifndef TATAMI_DELAYED_SUBSET_UTILS_HPP
#define TATAMI_DELAYED_SUBSET_UTILS_HPP

#include "../Matrix.hpp"
#include "../utils.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>

namespace tatami {

namespace subset_utils {

template<typename Value_, typename Index_>
const Value_* remap_dense(const Value_* input, Value_* buffer, const std::vector<Index_>& rmapping) {
    auto temp = buffer;
    for (auto i : rmapping) {
        *temp = input[i];
        ++temp;
    } 
    return buffer;
}

template<DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_, class IndexStorage_>
struct PerpendicularExtractor : public Extractor<selection_, sparse_, Value_, Index_> {
    PerpendicularExtractor(std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > i, const IndexStorage_& in) : 
        internal(std::move(i)), indices(&in)
    {
        if constexpr(selection_ == DimensionSelectionType::FULL) {
            this->full_length = internal->full_length;
        } else if constexpr(selection_ == DimensionSelectionType::BLOCK) {
            this->block_start = internal->block_start;
            this->block_length = internal->block_length;
        } else {
            this->index_length = internal->index_length;
        }
    }

    const Index_* index_start() const {
        if constexpr(selection_ == DimensionSelectionType::INDEX) {
            return internal->index_start();
        } else {
            return NULL;
        }
    }

protected:
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > internal;
    const IndexStorage_* indices;
};

template<DimensionSelectionType selection_, typename Value_, typename Index_, class IndexStorage_>
struct DensePerpendicularExtractor : public PerpendicularExtractor<selection_, false, Value_, Index_, IndexStorage_> {
    DensePerpendicularExtractor(std::unique_ptr<Extractor<selection_, false, Value_, Index_> > i, const IndexStorage_& p) : 
        PerpendicularExtractor<selection_, false, Value_, Index_, IndexStorage_>(std::move(i), p) {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        return this->internal->fetch((*(this->indices))[i], buffer);
    }
};

template<DimensionSelectionType selection_, typename Value_, typename Index_, class IndexStorage_>
struct SparsePerpendicularExtractor : public PerpendicularExtractor<selection_, true, Value_, Index_, IndexStorage_> {
    SparsePerpendicularExtractor(std::unique_ptr<Extractor<selection_, true, Value_, Index_> > i, const IndexStorage_& p) : 
        PerpendicularExtractor<selection_, true, Value_, Index_, IndexStorage_>(std::move(i), p) {}

    SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        return this->internal->fetch((*(this->indices))[i], vbuffer, ibuffer);
    }
};

template<bool accrow_, DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_, class IndexStorage_, typename ... Args_>
std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > populate_perpendicular(
    const Matrix<Value_, Index_>* mat, 
    const IndexStorage_& indices, 
    const Options<Index_>& options, 
    Args_... args)
{
    // TODO: handle variable access patterns here.
    std::unique_ptr<Extractor<selection_, sparse_, Value_, Index_> > output;

    if constexpr(sparse_) {
        output.reset(new SparsePerpendicularExtractor<selection_, Value_, Index_, IndexStorage_>(new_extractor<accrow_, sparse_>(mat, options, std::move(args)...), indices));
    } else {
        output.reset(new DensePerpendicularExtractor<selection_, Value_, Index_, IndexStorage_>(new_extractor<accrow_, sparse_>(mat, options, std::move(args)...), indices));
    }
    return output;
}

}

}

#endif
