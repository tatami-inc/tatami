#ifndef SCRAN_LOAD_MATRIX_MARKET_HPP
#define SCRAN_LOAD_MATRIX_MARKET_HPP

#include <iostream>
#include <fstream>
#include <limits>
#include <cstdint>
#include <algorithm>
#include <vector>

#include "../base/CompressedSparseMatrix.hpp"
#include "../base/DelayedBind.hpp"
#include "../utils/compress_sparse_triplets.hpp"

namespace tatami {

namespace MatrixMarket {

template<class STREAM>
std::tuple<size_t, size_t, size_t> process_headers(STREAM& input) {
    // Ignore comment headers. TODO: check it's actually coordinate integer, or
    // if it's coordinate real.
    while (input.peek() == '%') {
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Read number of rows and columns
    size_t num_row, num_col, num_lines;
    input >> num_row>> num_col >> num_lines;
    
    return std::make_tuple(num_row, num_col, num_lines);
}

struct LineAssignments {
    LineAssignments(size_t n) : category(n), index(n), nfeatures(3), nlines(3), permutation(n) {}
    std::vector<uint8_t> category;
    std::vector<size_t> index;
    std::vector<size_t> nfeatures;
    std::vector<size_t> nlines;
    std::vector<size_t> permutation;
};

template<class STREAM>
LineAssignments assign_lines(STREAM& input) {
    auto details = process_headers(input);
    auto num_rows = std::get<0>(details);
    auto num_lines = std::get<2>(details);

    LineAssignments output(num_rows);
    auto& flags = output.category;
    std::vector<size_t> counter(num_rows);
    constexpr size_t max8 = std::numeric_limits<uint8_t>::max();
    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    
    // Assigning each link to a block based on its integer size.
    for (size_t l = 0; l < num_lines; ++l) {
        double data;
        int row, col;
        input >> row >> col >> data;
        
        if (data > max16) {
            flags[row] = std::max(flags[row], static_cast<uint8_t>(2));
        } else if (data > max8) {
            flags[row] = std::max(flags[row], static_cast<uint8_t>(1));
        }
        ++counter[row];
    }

    // Computing the number of features and lines in each block.
    auto iIt = output.index.begin();
    auto cIt = counter.begin();
    for (auto f : flags) {
        auto& current = output.nfeatures[f];

        (*iIt) = current;
        ++current;
        ++iIt;

        output.nlines[f] += *cIt;
        ++cIt;
    }

    // Computing the permutation.
    auto cumsum = output.nfeatures;
    size_t last = 0;
    for (auto& x : cumsum) {
        std::swap(x, last);
        last += x;
    }
    for (size_t i = 0; i < num_rows; ++i) {
        output.permutation[i] = cumsum[flags[i]] + output.index[i];
    }

    return output;
}

template<typename T, typename IDX, typename U, typename V, typename W>
std::shared_ptr<Matrix<T, IDX> > create_sparse_matrix(size_t nr, size_t nc, U& values, V& rows, W& cols) {
    auto indptrs = compress_sparse_triplets<false>(nr, nc, values, rows, cols);
    return std::shared_ptr<Matrix<T, IDX> >(new CompressedSparseColumnMatrix<T, IDX, U, V, decltype(indptrs)>(nr, nc, std::move(values), std::move(rows), std::move(indptrs)));
}

template<typename T = double, typename IDX = int, typename ROW, class STREAM>
std::shared_ptr<Matrix<T, IDX> > create_layered_sparse_matrix_internal(STREAM& input, size_t num_cols, size_t num_lines, const LineAssignments& assign) {
    std::vector<ROW>      row8(assign.nlines[0]);
    std::vector<uint32_t> col8(assign.nlines[0]);
    std::vector<uint8_t>  dat8(assign.nlines[0]);
    size_t counter8 = 0;

    std::vector<ROW>      row16(assign.nlines[1]);
    std::vector<uint32_t> col16(assign.nlines[1]);
    std::vector<uint16_t> dat16(assign.nlines[1]);
    size_t counter16 = 0;

    std::vector<ROW>      row32(assign.nlines[2]);
    std::vector<uint32_t> col32(assign.nlines[2]);
    std::vector<uint32_t> dat32(assign.nlines[2]);
    size_t counter32 = 0;

    for (size_t l = 0; l < num_lines; ++l) {
        double data;
        int row, col;
        input >> row >> col >> data;

        auto cat = assign.category[row];
        auto idx = assign.index[row];

        if (cat == 0) {
            row8[counter8] = idx;
            col8[counter8] = col;
            dat8[counter8] = data;
            ++counter8;
        } else if (cat == 1) {
            row16[counter16] = idx;
            col16[counter16] = col;
            dat16[counter16] = data;
            ++counter16;
        } else if (cat == 2) {
            row32[counter32] = idx;
            col32[counter32] = col;
            dat32[counter32] = data;
            ++counter32;
        }
    }

    std::vector<std::shared_ptr<Matrix<T, IDX> > > collated;
    if (assign.nfeatures[0]) {
        collated.push_back(create_sparse_matrix<T, IDX>(assign.nfeatures[0], num_cols, dat8, row8, col8));
    }
    if (assign.nfeatures[1]) {
        collated.push_back(create_sparse_matrix<T, IDX>(assign.nfeatures[1], num_cols, dat16, row16, col16));
    }
    if (assign.nfeatures[2]) {
        collated.push_back(create_sparse_matrix<T, IDX>(assign.nfeatures[2], num_cols, dat32, row32, col32));
    }

    if (collated.size() == 0) {
        return create_sparse_matrix<T, IDX>(0, num_cols, dat8, row8, col8);
    } else if (collated.size() == 1) { 
        return collated[0];
    } else {
        return make_DelayedBind<0>(std::move(collated));
    }
}

template<typename T = double, typename IDX = int, class STREAM>
std::shared_ptr<Matrix<T, IDX> > create_layered_sparse_matrix(STREAM& input, const LineAssignments& assign) {
    auto output = process_headers(input);
    auto num_rows = std::get<0>(output);
    auto num_cols = std::get<1>(output);
    auto num_lines = std::get<2>(output);

    constexpr size_t max16 = std::numeric_limits<uint16_t>::max();
    if (num_rows > max16) {
        return create_layered_sparse_matrix_internal<T, IDX, uint32_t>(input, num_cols, num_lines, assign);
    } else {
        return create_layered_sparse_matrix_internal<T, IDX, uint16_t>(input, num_cols, num_lines, assign);
    }
}

template<typename T = double, typename IDX = int>
struct LayeredMatrixData {
    std::shared_ptr<Matrix<T, IDX> > matrix;
    std::vector<size_t> permutation;
};

template<typename T = double, typename IDX = int>
LayeredMatrixData<T, IDX> load_layered_sparse_matrix(const char * filepath) {
    std::ifstream in(filepath);
    auto ass = assign_lines(in);
    in.close();

    LayeredMatrixData<T, IDX> output;
    output.permutation = ass.permutation;

    in.open(filepath);
    output.matrix = create_layered_sparse_matrix<T, IDX>(in, ass); 

    return output;
}


}

}

#endif
