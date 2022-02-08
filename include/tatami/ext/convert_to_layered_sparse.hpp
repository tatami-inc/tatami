#ifndef TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP
#define TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP

#include "LayeredMatrixData.hpp"
#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

#include <cstdint>
#include <vector>
#include <memory>
#include <limits>

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

template<typename RowIndex, typename DataOut = double, typename IndexOut = int, typename DataIn, typename IndexIn>
RowRemapping<RowIndex> convert_by_row(
    const Matrix<DataIn, IndexIn>* incoming,
    std::vector<RowIndex>& row8,
    std::vector<uint8_t>& dat8,
    std::vector<size_t>& ptr8,
    std::vector<RowIndex>& row16,
    std::vector<uint32_t>& dat16,
    std::vector<size_t>& ptr16,
    std::vector<RowIndex>& row32,
    std::vector<uint32_t>& dat32,
    std::vector<size_t>& ptr32) 
{
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    if (!NC) {
        return create_empty<RowIndex>(NR);
    }

    auto wrk = incoming->new_workspace(true);
    std::vector<int> category(NR);
    std::vector<DataIn> dbuffer(NC);
    std::vector<IndexIn> ibuffer;

    auto ptr8p1 = ptr8.data() + 1; // making sure additions are offset by 1, for easier cumsumming.
    auto ptr16p1 = ptr16.data() + 1;
    auto ptr32p1 = ptr32.data() + 1;

    // First pass: categorizing each row and making a running tally of the entries per column.
    auto update_pointer = [&](size_t column, int chosen) -> void {
        switch (chosen) {
            case 0:
                ++(ptr8p1[column]);
                break;
            case 1:
                ++(ptr16p1[column]);
                break;
            case 2:
                ++(ptr32p1[column]);
                break;
            default:
                break;
        }
    };

    if (incoming->sparse()) {
        ibuffer.resize(NC);

        for (size_t i = 0; i < NR; ++i) {
            auto out = incoming->sparse_row(i, dbuffer.data(), ibuffer.data(), wrk.get());
            if (*std::min_element(out.value, out.value + out.number) < 0) {
                throw std::runtime_error("all values in the input matrix must be non-negative");
            }

            auto maxs = *std::max_element(out.value, out.value + out.number);
            int chosen = categorize_row(maxs);
            category[i] = chosen;

            for (size_t j = 0; j < out.number; ++j) {
                if (out.value[j]) {
                    update_pointer(out.index[j], chosen);
                }
            }
        }

    } else {
        for (size_t i = 0; i < NR; ++i) {
            auto out = incoming->row(i, dbuffer.data(), wrk.get());
            if (*std::min_element(out, out + NC) < 0) {
                throw std::runtime_error("all values in the input matrix must be non-negative");
            }

            auto maxs = *std::max_element(out, out + NC);
            int chosen = categorize_row(maxs);
            category[i] = chosen;

            for (size_t j = 0; j < NC; ++j) {
                if (out[j]) {
                    update_pointer(j, chosen);
                }
            }
        }
    }

    // Computing offsets.
    counts_to_offsets(ptr8);
    row8.resize(ptr8.back());
    dat8.resize(ptr8.back());

    counts_to_offsets(ptr16);
    row16.resize(ptr16.back());
    dat16.resize(ptr16.back());

    counts_to_offsets(ptr32);
    row32.resize(ptr32.back());
    dat32.resize(ptr32.back());

    auto remapping = compute_new_indices<RowIndex>(category);
    const auto& new_indices = remapping.new_indices;

    auto sofar8 = ptr8;
    auto sofar16 = ptr16;
    auto sofar32 = ptr32;

    // Second pass: categorizing each row and making a running tally of the entries per column.
    auto store_values = [&](DataIn value, RowIndex row, IndexIn column, int chosen) -> void {
        switch (chosen) {
            case 0:
                {
                    auto& sofar = sofar8[column];
                    dat8[sofar] = value;
                    row8[sofar] = row;
                    ++sofar;
                }
                break;
            case 1:
                {
                    auto& sofar = sofar16[column];
                    dat16[sofar] = value;
                    row16[sofar] = row;
                    ++sofar;
                }
                break;
            case 2:
                {
                    auto& sofar = sofar32[column];
                    dat32[sofar] = value;
                    row32[sofar] = row;
                    ++sofar;
                }
                break;
            default:
                break;
        }
    };

    if (incoming->sparse()) {
        for (size_t i = 0; i < NR; ++i) {
            auto out = incoming->sparse_row(i, dbuffer.data(), ibuffer.data(), wrk.get());
            auto r = new_indices[i];
            auto chosen = category[i];
            for (size_t j = 0; j < out.number; ++j) {
                if (out.value[j]) {
                    store_values(out.value[j], r, out.index[j], chosen);
                }
            }
        }

    } else {
        for (size_t i = 0; i < NR; ++i) {
            auto out = incoming->row(i, dbuffer.data(), wrk.get());
            auto r = new_indices[i];
            auto chosen = category[i];
            for (size_t j = 0; j < NC; ++j) {
                if (out[j]) {
                    store_values(out[j], r, j, chosen);
                }
            }
        }
    }

    return remapping;
}

template<typename RowIndex, typename DataOut = double, typename IndexOut = int, typename DataIn, typename IndexIn>
RowRemapping<RowIndex> convert_by_column(
    const Matrix<DataIn, IndexIn>* incoming,
    std::vector<RowIndex>& row8,
    std::vector<uint8_t>& dat8,
    std::vector<size_t>& ptr8,
    std::vector<RowIndex>& row16,
    std::vector<uint32_t>& dat16,
    std::vector<size_t>& ptr16,
    std::vector<RowIndex>& row32,
    std::vector<uint32_t>& dat32,
    std::vector<size_t>& ptr32) 
{
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();
    if (!NC) {
        return create_empty<RowIndex>(NR);
    }

    auto wrk = incoming->new_workspace(false);
    std::vector<DataIn> dbuffer(NR);
    std::vector<IndexIn> ibuffer;
    std::vector<DataIn> mins(NR), maxs(NR);
    std::vector<size_t> count(NR);

    // First pass: extracting the range per row and the number of non-zeros in each column.
    if (incoming->sparse()) {
        ibuffer.resize(NR);
        for (size_t i = 0; i < NC; ++i) {
            auto out = incoming->sparse_column(i, dbuffer.data(), ibuffer.data(), wrk.get());
            if (i) {
                for (size_t j = 0; j < out.number; ++j) {
                    auto val = out.value[j];
                    auto idx = out.index[j];
                    if (mins[idx] > val) {
                        mins[idx] = val;
                    } else if (maxs[idx] < val) {
                        maxs[idx] = val;
                    }
                    count[idx] += (val != 0);
                }
            } else {
                for (size_t j = 0; j < out.number; ++j) {
                    auto val = out.value[j];
                    auto idx = out.index[j];
                    mins[idx] = val;
                    maxs[idx] = val;
                    count[idx] += (val != 0);
                }
            }
        }
        
    } else {
        for (size_t i = 0; i < NC; ++i) {
            auto out = incoming->column(i, dbuffer.data(), wrk.get());
            if (i) {
                for (size_t j = 0; j < NR; ++j) {
                    auto val = out[j];
                    if (mins[j] > val) {
                        mins[j] = val;
                    } else if (maxs[j] < val) {
                        maxs[j] = val;
                    }
                    count[j] += (val != 0);
                }
            } else {
                for (size_t j = 0; j < NR; ++j) {
                    auto val = out[j];
                    mins[j] = val;
                    maxs[j] = val;
                    count[j] += (val != 0);
                }
            }
        }
    }
    
    // Categorizing each row.
    std::vector<int> category(NR);
    std::vector<size_t> totals(3);
    for (size_t i = 0; i < NR; ++i) {
        if (mins[i] < 0) {
            throw std::runtime_error("all values in the input matrix must be non-negative");
        }

        auto chosen = categorize_row(maxs[i]);
        category[i] = chosen;
        totals[chosen] += count[i];
    }

    row8.reserve(totals[0]);
    dat8.reserve(totals[0]);
    row16.reserve(totals[1]);
    dat16.reserve(totals[1]);
    row32.reserve(totals[2]);
    dat32.reserve(totals[2]);

    auto remapping = compute_new_indices<RowIndex>(category);
    const auto& new_indices = remapping.new_indices;

    auto ptr8p1 = ptr8.data() + 1; // making sure additions are offset by 1, for easier cumsumming.
    auto ptr16p1 = ptr16.data() + 1;
    auto ptr32p1 = ptr32.data() + 1;

    auto store_values = [&](DataIn value, IndexIn row, size_t column, int chosen) -> void {
        auto r = new_indices[row];
        switch (chosen) {
            case 0:
                dat8.push_back(value);
                row8.push_back(r);
                ++(ptr8p1[column]);
                break;
            case 1:
                dat16.push_back(value);
                row16.push_back(r);
                ++(ptr16p1[column]);
                break;
            case 2:
                dat32.push_back(value);
                row32.push_back(r);
                ++(ptr32p1[column]);
                break;
            default:
                break;
        }
    };

    if (incoming->sparse()) {
        for (size_t i = 0; i < NC; ++i) {
            auto out = incoming->sparse_column(i, dbuffer.data(), ibuffer.data(), wrk.get());
            for (size_t j = 0; j < out.number; ++j) {
                if (out.value[j]) {
                    store_values(out.value[j], out.index[j], i, category[out.index[j]]);
                }
            }
        }

    } else {
        for (size_t i = 0; i < NC; ++i) {
            auto out = incoming->column(i, dbuffer.data(), wrk.get());
            for (size_t j = 0; j < NR; ++j) {
                if (out[j]) {
                    store_values(out[j], j, i, category[j]);
                }
            }
        }
    }

    counts_to_offsets(ptr8);
    counts_to_offsets(ptr16);
    counts_to_offsets(ptr32);

    return remapping;
}

template<typename RowIndex, typename DataOut = double, typename IndexOut = int, typename DataIn, typename IndexIn>
LayeredMatrixData<DataOut, IndexOut> convert_internal(const Matrix<DataIn, IndexIn>* incoming) {
    size_t NR = incoming->nrow();
    size_t NC = incoming->ncol();

    std::vector<RowIndex> row8;
    std::vector<uint8_t>  dat8;
    std::vector<size_t>   ptr8(NC + 1);

    std::vector<RowIndex> row16;
    std::vector<uint32_t> dat16;
    std::vector<size_t>   ptr16(NC + 1);

    std::vector<RowIndex> row32;
    std::vector<uint32_t> dat32;
    std::vector<size_t>   ptr32(NC + 1);

    RowRemapping<RowIndex> remap;
    if (incoming->prefer_rows()) {
        remap = convert_by_row(
            incoming,
            row8, dat8, ptr8,
            row16, dat16, ptr16,
            row32, dat32, ptr32);
    } else {
         remap = convert_by_column(
            incoming,
            row8, dat8, ptr8,
            row16, dat16, ptr16,
            row32, dat32, ptr32);
    }

    // Creating the matrix.
    std::vector<std::shared_ptr<Matrix<DataOut, IndexOut> > > collated;
    const auto& per_category = remap.per_category;

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

    LayeredMatrixData<DataOut, IndexOut> output;
    output.permutation.swap(remap.permutation);

    if (collated.size() == 0) {
        typedef CompressedSparseColumnMatrix<DataOut, IndexOut, decltype(dat8), decltype(row8), decltype(ptr8)> CSCMatrix;
        output.matrix.reset(new CSCMatrix(0, NC, std::move(dat8), std::move(row8), std::move(ptr8)));
    } else if (collated.size() == 1) { 
        output.matrix = collated[0];
    } else {
        output.matrix = make_DelayedBind<0>(std::move(collated));
    }

    return output;
}
/**
 * @endcond
 */

}

/**
 * @param incoming A `tatami::Matrix` object containing non-negative integers.
 *
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of data value for the output `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 * @tparam Matrix Realized `tatami::Matrix` class to be converted.
 *
 * This function converts an existing sparse integer matrix into a layered sparse matrix.
 * The aim is to reduce memory usage by storing each gene's counts in the smallest unsigned integer type that can hold them.
 * See the documentation for `LayeredMatrixData` for more details.
 *
 * On a related note, this function will automatically use `uint16_t` values to store the internal row indices if the number of rows is less than 65536.
 * This aims to further reduce memory usage in most cases, e.g., gene count matrices usually have fewer than 50000 rows.
 * Note that the internal storage is orthogonal to the choice of `IDX` in the `tatami::Matrix` interface.
 */
template<typename T = double, typename IDX = int, class Matrix>
LayeredMatrixData<T, IDX> convert_to_layered_sparse(const Matrix* incoming) {
    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (incoming->nrow() <= max16) {
        return layered_utils::convert_internal<uint16_t, T, IDX>(incoming);
    } else {
        return layered_utils::convert_internal<IDX, T, IDX>(incoming);
    }
}

}

#endif
