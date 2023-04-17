#ifndef TATAMI_STANDARD_EXTRACTOR_HPP
#define TATAMI_STANDARD_EXTRACTOR_HPP

#include "Extractor.hpp"
#include <vector>
#include <type_traits>

namespace tatami {

template<bool sparse_, typename Value_, typename Index_>
struct IndexStartExtractor : public Extractor<sparse_, Value_, Index_> {
protected:
    const Index_* extracted_index_start = NULL;

    const Index_* quick_extracted_index() const { 
        return extracted_index_start;
    }
};

template<bool sparse_, typename Value_, typename Index_>
struct IndicesExtractor : public Extractor<sparse_, Value_, Index_> {
protected:
    std::vector<Index_> extracted_indices;

    const Index_* quick_extracted_index() const { 
        return extracted_indices.data();
    }
};

template<bool use_start_, bool sparse_, typename Value_, typename Index_>
struct IndexExtractor : public std::conditional<use_start_, IndexStartExtractor<sparse_, Value_, Index_>, IndicesExtractor<sparse_, Value_, Index_> >::type {
public:
    IndexExtractor(ExtractionOptions<Index_>& options) {
        this->extracted_selection = DimensionSelectionType::INDEX;
        if constexpr(use_start_) {
            this->extracted_length = options.selection.index_length;
            this->extracted_index_start = options.selection.index_start;
        } else {
            this->extracted_length = options.selection.indices.size();
            this->extracted_indices = std::move(options.selection.indices);
        }
    }

    const Index_* extracted_index() const { 
        return this->quick_extracted_index();
    }
};

template<DimensionSelectionType selection_, bool sparse_, typename Value_, typename Index_>
struct SimpleExtractor : public Extractor<sparse_, Value_, Index_> {
    SimpleExtractor(ExtractionOptions<Index_>& options) {
        this->extracted_selection = options.selection.type;
        if constexpr(selection_ == DimensionSelectionType::BLOCK) {
            this->extracted_length = options.selection.block_length;
            this->extracted_block = options.selection.block_start;
        }
    }

    const Index_* extracted_index() const { return NULL; }
};

template<DimensionSelectionType selection_, bool use_start_, bool sparse_, typename Value_, typename Index_>
using StandardExtractor = typename std::conditional<
        selection_ == DimensionSelectionType::INDEX, 
        IndexExtractor<use_start_, sparse_, Value_, Index_>,
        SimpleExtractor<selection_, sparse_, Value_, Index_>
    >::type;

}

#endif
