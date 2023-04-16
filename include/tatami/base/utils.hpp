#ifndef TATAMI_BASE_UTILS_HPP
#define TATAMI_BASE_UTILS_HPP

#include "Matrix.hpp"
#include "Options.hpp"
#include <memory>

namespace tatami {

/**
 * @cond
 */
// Default to false.
template<typename T, class V, typename = int>
struct has_data {
    static const bool value = false;
};

// Specialization is only run if it _has_ a data method.
template<typename T, class V>
struct has_data<T, V, decltype((void) V().data(), 0)> { 
    static const bool value = std::is_same<T*, decltype(V().data())>::value;
};
/**
 * @endcond
 */

template<bool accrow_, bool sparse_, typename Value_, typename Index_>
struct StandardExtractor : public Extractor<sparse_, Value_, Index_> {
    StandardExtractor(const Matrix<Value_, Index_>* matrix, ExtractionOptions<Index_>& options) {
        this->extracted_selection = options.selection.type;
        switch (this->extracted_selection) {
            case DimensionSelectionType::FULL:
                this->extracted_length = (accrow_ ? matrix->ncol() : matrix->nrow());
                break;
            case DimensionSelectionType::BLOCK:
                this->extracted_length = options.selection.block_length;
                this->extracted_block = options.selection.block_start;
                break;
            case DimensionSelectionType::INDEX:
                if (options.selection.index_start) {
                    this->extracted_length = options.selection.index_length;
                    this->extracted_index_start = options.selection.index_start;
                } else {
                    this->extracted_length = options.selection.indices.size();
                    this->extracted_indices = std::move(options.selection.indices);
                }
                break;
        }
    }

    const Index_* extracted_index() const { 
        return quick_extracted_index();
    }

protected:
    const Index_* extracted_index_start = NULL;

    std::vector<Index_> extracted_indices;

    const Index_* quick_extracted_index() const { 
        return (extracted_index_start ? extracted_index_start : extracted_indices.data()); 
    }
};

/**
 * @tparam row_ Whether to iterate over rows.
 * @tparam sparse_ Whether to perform sparse retrieval.
 * @tparam Value_ Data value type, should be numeric.
 * @tparam Index_ Row/column index type, should be integer.
 * @tparam Args_ Further arguments.
 *
 * @param[in] ptr Pointer to a `Matrix` object to iterate over.
 * @param args Zero or more additional arguments to pass to methods like `Matrix::dense_row()`, e.g., `IterationOptions`, `ExtractionOptions`.
 *
 * @return A `DimensionAccess` object to access the requested dimension of `ptr`.
 */
template<bool row_, bool sparse_, typename Value_, typename Index_, typename ... Args_>
std::unique_ptr<Extractor<sparse_, Value_, Index_> > new_extractor(const Matrix<Value_, Index_>* ptr, Args_... args) {
    if constexpr(sparse_) {
        if constexpr(row_) {
            return ptr->sparse_row(std::move(args)...);
        } else {
            return ptr->sparse_column(std::move(args)...);
        }
    } else {
        if constexpr(row_) {
            return ptr->dense_row(std::move(args)...);
        } else {
            return ptr->dense_column(std::move(args)...);
        }
    }
}

}

#endif
