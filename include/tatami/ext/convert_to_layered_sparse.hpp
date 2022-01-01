#ifndef TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP
#define TATAMI_CONVERT_TO_LAYERED_SPARSE_HPP

#include "LayeredMatrixData.hpp"
#include "../stats/ranges.hpp"
#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

#include <cstdint>
#include <vector>
#include <memory>

namespace tatami {

/**
 * @cond
 */
template<typename ROW, typename T = double, typename IDX = int>
LayeredMatrixData<T, IDX> convert_to_layered_sparse_internal(const Matrix<T, IDX>* incoming) {
    auto ranges = row_ranges(incoming);
    const auto& mins = ranges.first;
    const auto& maxs = ranges.second;

    // Using the permutation vector to store the reindexing vector for the time being.
    LayeredMatrixData<T, IDX> output;
    auto& reindexed = output.permutation;
    reindexed.resize(maxs.size());

    // Choosing a category for each row.
    std::vector<uint8_t> category(maxs.size());
    constexpr double max8 = std::numeric_limits<uint8_t>::max();
    constexpr double max16 = std::numeric_limits<uint16_t>::max();
    std::vector<size_t> per_category(3);
    
    for (size_t i = 0; i < maxs.size(); ++i) {
        if (mins[i] < 0) {
            throw std::runtime_error("all values must be non-negative");
        }
        if (maxs[i] <= max8) {
            category[i] = 0;
        } else if (maxs[i] <= max16) {
            category[i] = 1;
        } else {
            category[i] = 2;
        }

        auto& sofar = per_category[category[i]];
        reindexed[i] = sofar;
        ++sofar;
    }

    // Filling the categories.
    std::vector<ROW>      row8;
    std::vector<uint32_t> col8;
    std::vector<uint8_t>  dat8;

    std::vector<ROW>      row16;
    std::vector<uint32_t> col16;
    std::vector<uint16_t> dat16;

    std::vector<ROW>      row32;
    std::vector<uint32_t> col32;
    std::vector<uint32_t> dat32;

    auto add_nonzero_element = [&](size_t r, size_t c, T val) -> void {
        auto cat = category[r];
        auto r2 = reindexed[r];
        switch (cat) {
        case 0:
            row8.push_back(r2);
            col8.push_back(c);
            dat8.push_back(val);
            break;
        case 1:
            row16.push_back(r2);
            col16.push_back(c);
            dat16.push_back(val);
            break;
        case 2:
            row32.push_back(r2);
            col32.push_back(c);
            dat32.push_back(val);
            break;
        }
    };

    if (incoming->sparse()) {
        if (incoming->prefer_rows()) {
            auto wrk = incoming->new_workspace(true);
            std::vector<IDX> ibuffer(incoming->ncol());
            std::vector<T> vbuffer(incoming->ncol());

            for (size_t r = 0; r < incoming->nrow(); ++r) {
                auto range = incoming->sparse_row(r, vbuffer.data(), ibuffer.data(), wrk.get());
                for (size_t i = 0; i < range.number; ++i) {
                    if (range.value[i]) {
                        add_nonzero_element(r, range.index[i], range.value[i]);
                    }
                }
            }

        } else {
            auto wrk = incoming->new_workspace(false);
            std::vector<IDX> ibuffer(incoming->nrow());
            std::vector<T> vbuffer(incoming->nrow());

            for (size_t c = 0; c < incoming->ncol(); ++c) {
                auto range = incoming->sparse_column(c, vbuffer.data(), ibuffer.data(), wrk.get());
                for (size_t i = 0; i < range.number; ++i) {
                    if (range.value[i]) {
                        add_nonzero_element(range.index[i], c, range.value[i]);
                    }
                }
            }
        }

    } else {
        if (incoming->prefer_rows()) {
            auto wrk = incoming->new_workspace(true);
            std::vector<T> buffer(incoming->ncol());

            for (size_t r = 0; r < incoming->nrow(); ++r) {
                auto ptr = incoming->row(r, buffer.data(), wrk.get());
                for (size_t c = 0; c < incoming->ncol(); ++c) {
                    if (ptr[c] != 0) {
                        add_nonzero_element(r, c, ptr[c]);
                    }
                }
            }

        } else {
            auto wrk = incoming->new_workspace(false);
            std::vector<T> buffer(incoming->nrow());

            for (size_t c = 0; c < incoming->ncol(); ++c) {
                auto ptr = incoming->column(c, buffer.data(),  wrk.get());
                for (size_t r = 0; r < incoming->nrow(); ++r) {
                    if (ptr[r] != 0) {
                        add_nonzero_element(r, c, ptr[r]);
                    }
                }
            }
        }
    }

    // Creating the matrix.
    std::vector<std::shared_ptr<Matrix<T, IDX> > > collated;
    size_t NC = incoming->ncol();

    auto create_sparse_matrix = [](size_t nr, size_t nc, auto values, auto rows, auto cols) -> auto { 
        auto indptrs = compress_sparse_triplets<false>(nr, nc, values, rows, cols);
        typedef CompressedSparseColumnMatrix<T, IDX, decltype(values), decltype(rows), decltype(indptrs)> CSCMatrix;
        return std::shared_ptr<Matrix<T, IDX> >(new CSCMatrix(nr, nc, std::move(values), std::move(rows), std::move(indptrs)));
    };

    if (per_category[0]) {
        collated.push_back(create_sparse_matrix(per_category[0], NC, std::move(dat8), std::move(row8), std::move(col8)));
    }
    if (per_category[1]) {
        collated.push_back(create_sparse_matrix(per_category[1], NC, std::move(dat16), std::move(row16), std::move(col16)));
    }
    if (per_category[2]) {
        collated.push_back(create_sparse_matrix(per_category[2], NC, std::move(dat32), std::move(row32), std::move(col32)));
    }

    if (collated.size() == 0) {
        output.matrix = create_sparse_matrix(0, NC, std::move(dat8), std::move(row8), std::move(col8));
    } else if (collated.size() == 1) { 
        output.matrix = collated[0];
    } else {
        output.matrix = make_DelayedBind<0>(std::move(collated));
    }

    // Converting the reindexing vector into the permutation vector.
    {
        std::vector<size_t> offset(3);
        offset[1] = per_category[0];
        offset[2] = per_category[0] + per_category[1];
        for (size_t i = 0 ; i < category.size(); ++i) {
            output.permutation[i] += offset[category[i]];
        }
    }

    return output;
}
/**
 * @endcond
 */

/**
 * @param incoming A `tatami::Matrix` object containing non-negative integers.
 *
 * @return A `LayeredMatrixData` object.
 *
 * @tparam T Type of value in the `tatami::Matrix` interface.
 * @tparam IDX Integer type for the index.
 *
 * This function converts an existing sparse integer matrix into a layered sparse matrix.
 * The aim is to reduce memory usage by storing each gene's counts in the smallest unsigned integer type that can hold them.
 * See the documentation for `LayeredMatrixData` for more details.
 *
 * On a related note, this function will automatically use `uint16_t` values to store the internal row indices if the number of rows is less than 65536.
 * This aims to further reduce memory usage in most cases, e.g., gene count matrices usually have fewer than 50000 rows.
 * Note that the internal storage is orthogonal to the choice of `IDX` in the `tatami::Matrix` interface.
 */
template<typename T, typename IDX>
LayeredMatrixData<T, IDX> convert_to_layered_sparse(const Matrix<T, IDX>* incoming) {
    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (incoming->nrow() <= max16) {
        return convert_to_layered_sparse_internal<uint16_t, T, IDX>(incoming);
    } else {
        return convert_to_layered_sparse_internal<IDX, T, IDX>(incoming);
    }
}

}

#endif
