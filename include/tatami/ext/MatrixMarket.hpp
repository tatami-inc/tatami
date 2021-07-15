#ifndef SCRAN_LOAD_MATRIX_MARKET_HPP
#define SCRAN_LOAD_MATRIX_MARKET_HPP

#include <limits>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <fstream>

#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

namespace tatami {

namespace MatrixMarket {

template<typename T>
size_t process_triplet_line(const char* buffer, T& arg1, T& arg2, T& arg3) {
    const char * original = buffer;
    char *next;
    arg1 = std::strtol(buffer, &next, 10);
    buffer = next;

    arg2 = std::strtol(buffer, &next, 10);
    buffer = next;

    arg3 = std::strtol(buffer, &next, 10);
    buffer = next;

    // TODO: check for positive integer values.

    return buffer - original;
}

struct LineAssignments {
    LineAssignments() : rows_per_category(3), lines_per_category(3) {}

    std::vector<uint8_t> category;
    std::vector<size_t> index;
    std::vector<size_t> rows_per_category;
    std::vector<size_t> lines_per_category;
    std::vector<size_t> lines_per_row;
    std::vector<size_t> permutation;

    size_t nrows, ncols, nlines;
    bool passed_preamble = false;

    bool preamble(const char* buffer) {
        if (buffer[0] == '%') {
            // TODO: should probably check for 'coordinate integer'.
            return true;
        }

        process_triplet_line(buffer, nrows, ncols, nlines);
        category.resize(nrows);
        index.resize(nrows);
        permutation.resize(nrows);
        lines_per_row.resize(nrows);
        return false;
    }

    void add(const char* buffer) {
        // Assigning each line to a block based on its integer size.
        constexpr long max8 = std::numeric_limits<uint8_t>::max();
        constexpr long max16 = std::numeric_limits<uint16_t>::max();

        long row, col, data;
        process_triplet_line(buffer, row, col, data);
        
        if (data > max16) {
            category[row] = std::max(category[row], static_cast<uint8_t>(2));
        } else if (data > max8) {
            category[row] = std::max(category[row], static_cast<uint8_t>(1));
        }
        ++lines_per_row[row];

        return;
    }

    void finish() {
        // Computing the number of features and lines in each block.
        auto iIt = index.begin();
        auto cIt = lines_per_row.begin();
        for (auto f : category) {
            auto& current = rows_per_category[f];

            (*iIt) = current;
            ++current;
            ++iIt;

            lines_per_category[f] += *cIt;
            ++cIt;
        }

        // Computing the permutation.
        auto cumsum = rows_per_category;
        size_t last = 0;
        for (auto& x : cumsum) {
            std::swap(x, last);
            last += x;
        }
        for (size_t i = 0; i < nrows; ++i) {
            permutation[i] = cumsum[category[i]] + index[i];
        }

        return;
    }
};

template<typename ROW>
struct LayeredBuilder {
private:
    LineAssignments assign;

    std::vector<ROW>      row8;
    std::vector<uint32_t> col8;
    std::vector<uint8_t>  dat8;
    size_t counter8 = 0;

    std::vector<ROW>      row16;
    std::vector<uint32_t> col16;
    std::vector<uint16_t> dat16;
    size_t counter16 = 0;

    std::vector<ROW>      row32;
    std::vector<uint32_t> col32;
    std::vector<uint32_t> dat32;
    size_t counter32 = 0;
public:
    LayeredBuilder(LineAssignments ass) : 
        assign(std::move(ass)),
        row8(assign.lines_per_category[0]),
        col8(assign.lines_per_category[0]),
        dat8(assign.lines_per_category[0]),
        row16(assign.lines_per_category[1]),
        col16(assign.lines_per_category[1]),
        dat16(assign.lines_per_category[1]),
        row32(assign.lines_per_category[2]),
        col32(assign.lines_per_category[2]),
        dat32(assign.lines_per_category[2]) {}

public:
    bool preamble(const char* buffer) {
        if (buffer[0] == '%') {
            return true;
        }
        return false;
    }

    void add(const char* buffer) {
        long row, col, data;
        process_triplet_line(buffer, row, col, data);

        auto idx = assign.index[row];
        switch (assign.category[row]) {
            case 0:
                row8[counter8] = idx;
                col8[counter8] = col;
                dat8[counter8] = data;
                ++counter8;
                break;
            case 1:
                row16[counter16] = idx;
                col16[counter16] = col;
                dat16[counter16] = data;
                ++counter16;
                break;
            case 2:
                row32[counter32] = idx;
                col32[counter32] = col;
                dat32[counter32] = data;
                ++counter32;
                break;
        }
    }

private:
    template<typename T, typename IDX, typename U, typename V, typename W>
    std::shared_ptr<Matrix<T, IDX> > create_sparse_matrix(size_t nr, size_t nc, U& values, V& rows, W& cols) {
        auto indptrs = compress_sparse_triplets<false>(nr, nc, values, rows, cols);
        return std::shared_ptr<Matrix<T, IDX> >(new CompressedSparseColumnMatrix<T, IDX, U, V, decltype(indptrs)>(nr, nc, std::move(values), std::move(rows), std::move(indptrs)));
    }

public:
    template<typename T, typename IDX>
    std::shared_ptr<Matrix<T, IDX> > finish() {
        std::vector<std::shared_ptr<Matrix<T, IDX> > > collated;
        if (assign.rows_per_category[0]) {
            collated.push_back(create_sparse_matrix<T, IDX>(assign.rows_per_category[0], assign.ncols, dat8, row8, col8));
        }
        if (assign.rows_per_category[1]) {
            collated.push_back(create_sparse_matrix<T, IDX>(assign.rows_per_category[1], assign.ncols, dat16, row16, col16));
        }
        if (assign.rows_per_category[2]) {
            collated.push_back(create_sparse_matrix<T, IDX>(assign.rows_per_category[2], assign.ncols, dat32, row32, col32));
        }

        if (collated.size() == 0) {
            return create_sparse_matrix<T, IDX>(0, assign.ncols, dat8, row8, col8);
        } else if (collated.size() == 1) { 
            return collated[0];
        } else {
            return make_DelayedBind<0>(std::move(collated));
        }
    }
};

template<typename T = double, typename IDX = int>
struct LayeredMatrixData {
    std::shared_ptr<Matrix<T, IDX> > matrix;
    std::vector<size_t> permutation;
};

template<typename T = double, typename IDX = int>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix(const char * filepath) {
    auto process = [](auto& obj, auto& in) -> void {
        std::string line;
        while (std::getline(in, line) && obj.preamble(line.c_str())) {}
        while (std::getline(in, line)) {
            obj.add(line.c_str());
        }
        return;
    };

    std::ifstream in(filepath);
    LineAssignments ass;
    process(ass, in);
    ass.finish();

    LayeredMatrixData<T, IDX> output;
    output.permutation = ass.permutation;

    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (ass.nrows > max16) {
        std::ifstream in(filepath);
        LayeredBuilder<uint16_t> builder(std::move(ass));
        process(builder, in);
        output.matrix = builder.finish<T, IDX>();
    } else {
        std::ifstream in(filepath);
        LayeredBuilder<uint16_t> builder(std::move(ass));
        process(builder, in);
        output.matrix = builder.finish<T, IDX>();
    }

    return output;
}

}

}

#endif
