#ifndef TATAMI_COMPRESSED_SPARSE_SECONDARY_EXTRACTOR_BASIC_HPP
#define TATAMI_COMPRESSED_SPARSE_SECONDARY_EXTRACTOR_BASIC_HPP

#include <vector>
#include "SparseSecondaryExtractorCore.hpp"

namespace tatami {

template<typename Index_, typename StoredIndex_, typename CustomPointer_, class CustomPointerModifier_> 
struct CompressedSparseSecondaryExtractorBasic : public SparseSecondaryExtractorCore<Index_, StoredIndex_, CustomPointer_, CustomPointerModifier_> {
public:
    CompressedSparseSecondaryExtractorBasic() = default;

    template<class IndexStorage_, class PointerStorage_>
    CompressedSparseSecondaryExtractorBasic(StoredIndex_ max_index, const IndexStorage_& idx, const PointerStorage_& idp, Index_ start, Index_ length) :
        SparseSecondaryExtractorCore<Index_, StoredIndex_, CustomPointer_, CustomPointerModifier_>(max_index, length)
    {
        auto idpIt = idp.begin() + start;
        for (Index_ i = 0; i < length; ++i, ++idpIt) {
            this->current_indptrs[i] = *idpIt;
            this->current_indices[i] = (*idpIt < *(idpIt + 1) ? idx[*idpIt] : max_index);
        }
        this->closest_current_index = (length ? *std::min_element(this->current_indices.begin(), this->current_indices.end()) : max_index);
        return;
    } 

    template<class IndexStorage_, class PointerStorage_>
    CompressedSparseSecondaryExtractorBasic(StoredIndex_ max_index, const IndexStorage_& idx, const PointerStorage_& idp) :
        CompressedSparseSecondaryExtractorBasic(max_index, idx, idp, static_cast<Index_>(0), static_cast<Index_>(idp.size() - 1)) {}

    template<class IndexStorage_, class PointerStorage_>
    CompressedSparseSecondaryExtractorBasic(StoredIndex_ max_index, const IndexStorage_& idx, const PointerStorage_& idp, const Index_* subset, Index_ length) :
        SparseSecondaryExtractorCore<Index_, StoredIndex_, CustomPointer_, CustomPointerModifier_>(max_index, length)
    {
        for (Index_ i0 = 0; i0 < length; ++i0) {
            auto i = subset[i0];
            this->current_indptrs[i0] = idp[i];
            this->current_indices[i0] = (idp[i] < idp[i + 1] ? idx[idp[i]] : max_index);
        }
        this->closest_current_index = (length ? *std::min_element(this->current_indices.begin(), this->current_indices.end()) : max_index);
        return;
    }

public:
    template<class IndexStorage_, class PointerStorage_, class PrimaryFunction_, class StoreFunction_, class SkipFunction_>
    bool search(
        StoredIndex_ secondary,
        Index_ primary_length,
        PrimaryFunction_&& to_primary, 
        const IndexStorage_& indices,
        const PointerStorage_& indptrs,
        StoreFunction_&& store,
        SkipFunction_&& skip
    ) {
        return this->search_base(
            secondary, 
            primary_length, 
            std::forward<PrimaryFunction_>(to_primary), 
            indices, 
            indptrs, 
            std::forward<StoreFunction_>(store), 
            std::forward<SkipFunction_>(skip)
        );
    }
};

}

#endif
