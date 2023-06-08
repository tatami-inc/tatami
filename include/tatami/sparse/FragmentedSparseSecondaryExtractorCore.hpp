#ifndef TATAMI_FRAGMENTED_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP
#define TATAMI_FRAGMENTED_SPARSE_SECONDARY_EXTRACTOR_CORE_HPP

#include <vector>

namespace tatami {

template<typename Index_, typename StoredIndex_, typename CustomPointer_, class CustomPointerModifier_> 
struct FragmentedSparseSecondaryExtractorCore : public SparseSecondaryExtractorCore<Index_, StoredIndex_, CustomPointer_, CustomPointerModifier_> {
public:
    FragmentedSparseSecondaryExtractorCore() = default;

    FragmentedSparseSecondaryExtractorCore(StoredIndex_ mi, Index_ length) : SparseSecondaryExtractorCore<Index_, StoredIndex_, CustomPointer_, CustomPointerModifier_>(mi, length) {}

public:
    template<class IndexVectorStorage_, class PrimaryFunction_, class StoreFunction_, class SkipFunction_>
    bool search(
        StoredIndex_ secondary,
        Index_ primary_length,
        PrimaryFunction_ to_primary, 
        const IndexVectorStorage_& indices,
        StoreFunction_ store,
        SkipFunction_ skip
    ) {
        return this->search_base(
            secondary, 
            primary_length, 
            std::forward<PrimaryFunction_>(to_primary), 
            indices, 
            false, 
            std::forward<StoreFunction_>(store), 
            std::forward<SkipFunction_>(skip)
        );
    }
};

}

#endif
