#include "tatami/tatami.hpp"
#include <vector>
#include <iostream>
#include <numeric>
#include <iomanip>

/* INTRODUCTION: 
 *
 * In this example, we will demonstrate how the interface type does not
 * necessarily need to be the same as the storage type. This allows developers
 * to use a more memory-efficient representation (e.g., by using a smaller
 * integer type or by reducing the precision of floats) while still retaining
 * compatibility with downstream functions that expect a certain type.
 */

template<class IT>
void print_vector(IT start, IT end) {
    bool first = true;
    std::cout << "[ "; 
    for (IT it = start; it != end; ++it) {
        if (!first) {
            std::cout << ", ";
        }
        std::cout << std::setw(6) << std::fixed << std::setprecision(2) << *it;
        first = false;
    }
    std::cout << " ]" << std::endl;
}

int main(int argc, char** argv) {
    std::vector<uint8_t> rows = { 3, 5, 0, 1, 8, 4, 7, 5, 6, 9, 0, 1, 2, 3, 4 };
    std::vector<uint8_t> cols = { 0, 1, 3, 2, 0, 4, 1, 2, 4, 0, 1, 3, 0, 3, 2 };
    std::vector<char> vals = { -4, 1, -7, 12, 12, -1, -4, 2, 3, 4, -1, 5, -8, 1, 2 };
    auto indptrs = tatami::compress_sparse_triplets<false>(10, 5, vals, rows, cols);

    /* Here, we create a matrix where the values are chars and the indices are
     * 8-bit unsigned integers. This enables us to be extremely space-efficient
     * if we know that the values will not exceed the limits of the type.
     */
    std::shared_ptr<tatami::NumericMatrix> mat(new tatami::CompressedSparseColumnMatrix<double, int, decltype(vals), decltype(rows)>(10, 5, vals, rows, indptrs));

    std::cout << "Matrix preview: " << std::endl;
    auto wrk = mat->dense_row();
    std::vector<double> buffer(mat->ncol());
    for (size_t i = 0; i < mat->nrow(); ++i) {
        auto ptr = wrk->fetch(i, buffer.data());
        print_vector(ptr, ptr + mat->ncol());
    }
    std::cout << std::endl;

    /* However, this improved memory efficiency does come at the cost of mandatory
     * type conversions (and the associated implicit copy). For example, it is no
     * longer possible to return a direct pointer to the underlying store. (In 
     * contrast, if vals was a vector of doubles, the code below would print 'false'.)
     */
    std::vector<double> vbuffer(mat->nrow());
    std::vector<int> ibuffer(mat->nrow());
    auto swrk = mat->sparse_column();
    auto range = swrk->fetch(0, vbuffer.data(), ibuffer.data());
    std::cout << "Using buffer instead of underlying store: " << (range.value==vbuffer.data() ? "true" : "false") << std::endl;

    return 0;
}
