#ifndef TATAMI_BIND_INTERSECTION_HPP
#define TATAMI_BIND_INTERSECTION_HPP

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <memory>
#include "../subset/make_DelayedSubset.hpp"
#include "../other/DelayedBind.hpp"

/**
 * @file bind_intersection.hpp
 *
 * @brief Combine matrices based on the intersection of identifiers.
 */

namespace tatami {

/**
 * Combine multiple matrices along the specified dimension while accounting for non-identical IDs along the other dimension. 
 * This function will identify the intersection of common identifiers across all matrices, 
 * subset each individual matrix to the intersection in the other dimension, 
 * and then combine them along the specified dimension.
 *
 * @tparam margin The dimension to combine along - by rows (0) or columns (1).
 * @tparam Matrix The **tatami** matrix class.
 * @tparam Id Integer type representing the identifiers.
 *
 * @param inputs Vector of pointers to matrices to be combined along the dimension specified by `margin`.
 * @param ids Vector of length equal to that of `inputs`, containing pointers to the identifiers for the matrices.
 * Each pointer should refer to an array equal to the number of columns (if `margin = 0`) or rows (if `margin = 1`) in the corresponding entry of `inputs`.
 *
 * @return Pair containing the combined matrix (first) and a vector of indices for the intersection (second).
 * In the combined matrix, each column (if `margin = 0`) or row (if `margin = 1`) will correspond to an entry in the intersection vector. 
 * Indices in the intersection vector should be applied to `ids[0]` to obtain the actual identifiers of the rows/columns.
 */
template<int margin, class Matrix, typename Id>
std::pair<std::shared_ptr<Matrix>, std::vector<size_t> > bind_intersection(const std::vector<std::shared_ptr<Matrix> >& inputs, const std::vector<const Id*>& ids) {
    size_t n = inputs.size();
    if (n == 0) {
        throw std::runtime_error("at least one matrix must be supplied for binding");
    }
    if (n != ids.size()) {
        throw std::runtime_error("'inputs' and 'ids' should have the same length");
    }

    auto otherlen = [](const auto& ptr) -> size_t {
        if constexpr(margin) {
            return ptr->nrow();
        } else {
            return ptr->ncol();
        }
    };

    // Finding the intersection of all names and creating a mapping.
    const auto& first_ptr = inputs[0];
    auto first_ids = ids[0];
    std::unordered_set<int> in_use(first_ids, first_ids + otherlen(first_ptr));

    for (int i = 1; i < n; ++i) {
        std::vector<int> intersection;
        intersection.reserve(in_use.size());

        const auto& current = inputs[i];
        size_t current_other = otherlen(current);
        auto current_ids = ids[i];

        for (size_t j = 0; j < current_other; ++j) {
            if (in_use.find(current_ids[j]) != in_use.end()) {
                intersection.push_back(current_ids[j]);            
            }
        }

        in_use = std::unordered_set<int>(intersection.begin(), intersection.end());
    }

    std::unordered_map<int, int> mapping;
    std::vector<int> as_vec(in_use.begin(), in_use.end());
    std::sort(as_vec.begin(), as_vec.end());
    for (int s = 0; s < as_vec.size(); ++s) {
        mapping[as_vec[s]] = s;
    }

    // Applying the mapping. 
    std::vector<std::shared_ptr<Matrix> > collected;
    collected.reserve(inputs.size());
    std::vector<size_t> indices;

    for (int i = 0; i < n; ++i) {
        const auto& current = inputs[i];
        size_t current_other = otherlen(current);
        auto current_ids = ids[i];

        std::vector<size_t> reorder(mapping.size());
        for (size_t j = 0; j < current_other; ++j) {
            auto it = mapping.find(current_ids[j]);
            if (it != mapping.end()) {
                reorder[it->second] = j;
            }
        }

        if (i == 0) {
            indices = reorder;
        }

        collected.push_back(make_DelayedSubset<1 - margin>(current, std::move(reorder)));
    }

    return std::make_pair(make_DelayedBind<margin>(std::move(collected)), std::move(indices));
}

}

#endif
