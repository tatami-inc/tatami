#ifndef TATAMI_LAYERED_UTILS_HPP
#define TATAMI_LAYERED_UTILS_HPP

#include <limits>
#include <vector>
#include <numeric>

#include "../base/Matrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../base/CompressedSparseMatrix.hpp"

namespace tatami {

/**
 * @cond
 */
namespace layered_utils {

template<typename T>
int categorize_row(T v) {
    constexpr uint8_t max8 = std::numeric_limits<uint8_t>::max();
    constexpr uint16_t max16 = std::numeric_limits<uint16_t>::max();
    constexpr T maxv = std::numeric_limits<T>::max();

    if constexpr(maxv <= max8) {
        return 0;
    } else {
        if (v <= max8) {
            return 0;
        }
    }

    if constexpr(maxv <= max16) {
        return 1;
    } else {
        if (v <= max16) {
            return 1;
        }
    }

    return 2;
}

template<typename O>
struct RowRemapping {
    RowRemapping() {}
    RowRemapping(size_t n) : per_category(3), new_indices(n), permutation(n) {}

    std::vector<O> per_category;
    std::vector<O> new_indices;
    std::vector<size_t> permutation;
};

template<typename O>
RowRemapping<O> compute_new_indices(const std::vector<int>& category) {
    RowRemapping<O> output(category.size());
    for (size_t i = 0; i < category.size(); ++i) {
        auto& sofar = output.per_category[category[i]];
        output.new_indices[i] = sofar;
        ++sofar;
    }

    // Computing the permutation, while we're here.
    O offset[3];
    offset[0] = 0;
    offset[1] = output.per_category[0];
    offset[2] = offset[1] + output.per_category[1];

    for (size_t i = 0; i < category.size(); ++i){
        output.permutation[i] = offset[category[i]] + output.new_indices[i];
    }

    return output;
}

template<typename O>
RowRemapping<O> create_empty(size_t n) {
    RowRemapping<O> output(n);
    std::iota(output.new_indices.begin(), output.new_indices.end(), 0);
    std::iota(output.permutation.begin(), output.permutation.end(), 0);
    return output;
}

template<typename C>
void counts_to_offsets(std::vector<C>& counts) {
    for (size_t i = 1; i < counts.size(); ++i) {
        counts[i] += counts[i-1];
    }
    return;
}

template<typename DataOut, typename IndexOut>
std::shared_ptr<Matrix<DataOut, IndexOut> > consolidate_submatrices(std::vector<std::shared_ptr<Matrix<DataOut, IndexOut> > > collated, size_t ncol) {
    if (collated.size() == 0) {
        typedef CompressedSparseColumnMatrix<DataOut, IndexOut> CSCMatrix;
        return std::shared_ptr<Matrix<DataOut, IndexOut> >(new CSCMatrix(0, ncol, std::vector<DataOut>(), std::vector<IndexOut>(), std::vector<size_t>()));
    } else if (collated.size() == 1) {
        return collated[0];
    } else {
        return make_DelayedBind<0>(std::move(collated));
    }
}

template<typename DataOut, typename IndexOut, typename RowCount, typename RowIndex>
std::shared_ptr<Matrix<DataOut, IndexOut> > consolidate_submatrices(
    const std::vector<RowCount>& per_category,
    std::vector<RowIndex> row8,
    std::vector<uint8_t> dat8,
    std::vector<size_t> ptr8,
    std::vector<RowIndex> row16,
    std::vector<uint16_t> dat16,
    std::vector<size_t> ptr16,
    std::vector<RowIndex> row32,
    std::vector<uint32_t> dat32,
    std::vector<size_t> ptr32)
{
    std::vector<std::shared_ptr<Matrix<DataOut, IndexOut> > > collated;
    size_t NC = ptr8.size() - 1;

    if (per_category[0]) {
        typedef CompressedSparseColumnMatrix<DataOut, IndexOut, decltype(dat8), decltype(row8), decltype(ptr8)> CSCMatrix8;
        collated.emplace_back(new CSCMatrix8(per_category[0], NC, std::move(dat8), std::move(row8), std::move(ptr8)));
    }

    if (per_category[1]) {
        typedef CompressedSparseColumnMatrix<DataOut, IndexOut, decltype(dat16), decltype(row16), decltype(ptr16)> CSCMatrix16;
        collated.emplace_back(new CSCMatrix16(per_category[1], NC, std::move(dat16), std::move(row16), std::move(ptr16)));
    }

    if (per_category[2]) {
        typedef CompressedSparseColumnMatrix<DataOut, IndexOut, decltype(dat32), decltype(row32), decltype(ptr32)> CSCMatrix32;
        collated.emplace_back(new CSCMatrix32(per_category[2], NC, std::move(dat32), std::move(row32), std::move(ptr32)));
    }

    return consolidate_submatrices(std::move(collated), NC);
}

}
/**
 * @endcond
 */

}

#endif
