#include <random>
#include <vector>
#include <chrono>
#include <iostream>

#include "tatami/base/CompressedSparseMatrix.hpp"
#include "tatami/utils/compress_sparse_triplets.hpp"

/* INTRODUCTION:
 *
 * Here, we compare the performance with and without a workspace. The workspace
 * allows a matrix representation to re-use information across multiple calls.
 * For example, for row access from compressed sparse column matrices, we can
 * cache the last position of the index pointers to expedite the access to the
 * next (consecutive) row. File-backed representations can also store handles and
 * caches in the workspace for later retrieval.
 *
 * This example demonstrates the kind of speed-ups that are possible by using
 * workspaces.  We create a large compressed sparse column sparse matrix and
 * then record the time required to iterate across all rows and count the
 * number of non-zero elements.
 */

std::shared_ptr<tatami::NumericMatrix> generate_sparse_matrix(size_t nr, size_t nc, double density) {
    std::vector<int> i, j;
    std::vector<double> x;

    std::mt19937_64 generator(1234567);
    std::uniform_real_distribution<double> distu;
    std::normal_distribution<double> distn;

    for (size_t c = 0; c < nc; ++c) {
        for (size_t r = 0; r < nr; ++r) {
            if (distu(generator) <= density) {
                i.push_back(r);
                j.push_back(c);
                x.push_back(distn(generator));
            }
        }
    }

    // Get column-major representation.
    auto indptrs = tatami::compress_sparse_triplets<false>(nr, nc, x, i, j);
    return std::shared_ptr<tatami::NumericMatrix>(new tatami::CompressedSparseColumnMatrix<double, int>(nr, nc, std::move(x), std::move(i), std::move(indptrs)));
}

int main() {
    auto mat = generate_sparse_matrix(10000, 10000, 0.1);

    std::vector<double> xbuffer(mat->ncol());
    std::vector<int> ibuffer(mat->ncol());

    auto start = std::chrono::high_resolution_clock::now();
    auto wrk = mat->sparse_row_workspace();
    int sum = 0;
    for (size_t r = 0; r < mat->nrow(); ++r) {
        auto range = mat->row(r, xbuffer.data(), ibuffer.data(), wrk.get());
        sum += range.number;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Workspace access time: " << duration.count() << " for " << sum << " non-zero elements" << std::endl;

    // The latest version of tatami mandates a workspace for all calls, so the
    // times are somewhat inflated by repeated allocations for the workspace.
    start = std::chrono::high_resolution_clock::now();
    sum = 0;
    for (size_t r = 0; r < mat->nrow(); ++r) {
        auto wrk = mat->sparse_row_workspace(); 
        auto range = mat->row(r, xbuffer.data(), ibuffer.data(), wrk.get());
        sum += range.number;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "No workspace access time: " << duration.count() << " for " << sum << " non-zero elements" << std::endl;

    // Testing the inflation in time.
    start = std::chrono::high_resolution_clock::now();
    tatami::SparseRowWorkspace* ptr;
    for (size_t r = 0; r < mat->nrow(); ++r) {
        auto wrk = mat->sparse_row_workspace(); 
        ptr = wrk.get(); // avoid optimizing out the entire loop.
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "  (Inflation from workspace creation time:" << " " << duration.count() << ")" << std::endl; 

    return 0;
}
