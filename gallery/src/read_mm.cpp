#include <exception>
#include <iostream>

/**
 * INTRODUCTION:
 *
 * This is little utility that demonstrates how to call the  
 * matrix market readers. In this case, we're creating a layered
 * sparse matrix, but you can just create a regular sparse
 * matrix with 'load_sparse_matrix()`. Note that these functions
 * assume that we're dealing with non-negative count data in the
 * input file; floats and negative values are not supported.
 */

#include "tatami/ext/MatrixMarket.hpp"

int main(int argc, char* argv[]) {
    if (argc!=2) {
        std::cerr << "need one argument" << std::endl;
        return 1;
    }

    auto x = tatami::MatrixMarket::load_layered_sparse_matrix_from_file(argv[1]);
    std::cout << "Dimensions: " << x.matrix->nrow() << " x " << x.matrix->ncol() << std::endl;
    return 0;
}
