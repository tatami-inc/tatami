#include "tatami/tatami.h"
#include "print_vector.h"
#include <vector>

/* INTRODUCTION: 
 *
 * In this example, we will be computing the row sums of a tatami numeric
 * matrix, with a focus on seeing how to parallelize this. In this case, we are
 * using the OpenMP parallelization framework, but the same principles apply to
 * other mechanisms, e.g., with std::thread or Intel TBB.
 */

/* ATTEMPT 1:
 *
 * This is the simplest implementation whereby we just loop across all rows,
 * extract each row and compute the sum. We parallelize using the simple openMP
 * for pragma; nothing too dramatic happening here. This is thread-safe as all
 * tatami methods are read-only; there are no problems due to state changes
 * inside 'p'. Note that the buffer needs to be created per thread inside the
 * parallel region.
 */
std::vector<double> rowsums_simple(std::shared_ptr<tatami::numeric_matrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NR);

    # pragma omp parallel
    {
        std::vector<double> buffer(NC);

        #pragma omp for 
        for (size_t i = 0; i < NR; ++i) {
            auto ptr = p->row(i, buffer.data());
            output[i] = std::accumulate(ptr, ptr + NR, 0.0);
        }
    }

    return output;
}

/* ATTEMPT 2:
 *
 * A more complex approach is to use workspaces to share information between
 * calls. We create one workspace per thread within an OpenMP parallel region,
 * and then we distribute work acros the loop as usual. 
 *
 * Notice how we use a static schedule with no chunk size arguments. This is
 * because workspaces are generally most beneficial when dealing with
 * extraction of a block of consecutive rows within each thread. Our current
 * schedule ensures that each thread gets one chunk of consecutive jobs. Of
 * course, if extraction is cheap compared to the actual calculation, other
 * schedules may be preferable.
 */
std::vector<double> rowsums_work(std::shared_ptr<tatami::numeric_matrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NR);

    #pragma omp parallel
    {
        std::vector<double> buffer(NC);
        auto wrk = p->new_workspace(true);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < NR; ++i) {
            auto ptr = p->row(i, buffer.data(), wrk.get());
            output[i] = std::accumulate(ptr, ptr + NR, 0.0);
        }
    }

    return output;
}

int main(int argc, char** argv) {
    std::vector<int> rows = { 3, 5, 0, 1, 8, 4, 7, 5, 6, 9, 0, 1, 2, 3, 4 };
    std::vector<int> cols = { 0, 1, 3, 2, 0, 4, 1, 2, 4, 0, 1, 3, 0, 3, 2 };
    std::vector<double> vals = { -0.40, 0.14, -0.17, 1.20, 1.20, -1.10, -0.42, 2.10, 0.38, 0.40, -1.10, 0.57, -0.89, 1.60, 0.27 };

    auto indptrs = tatami::compress_sparse_triplets<false>(10, 5, vals, rows, cols);
    std::shared_ptr<tatami::numeric_matrix> mat(new tatami::CompressedSparseColumnMatrix<double, int>(10, 5, vals, rows, indptrs));

    std::cout << "Matrix preview: " << std::endl;
    std::vector<double> buffer(mat->ncol());
    for (size_t i = 0; i < mat->nrow(); ++i) {
        auto ptr = mat->row(i, buffer.data());
        print_vector(ptr, ptr + mat->ncol());
    }
    std::cout << std::endl;

    {
        std::cout << "Simple:" << std::endl;
        auto cs = rowsums_simple(mat); 
        print_vector(cs.begin(), cs.end());
        std::cout << std::endl;
    }

    {
        std::cout << "Workspaced:" << std::endl;
        auto cs = rowsums_work(mat); 
        print_vector(cs.begin(), cs.end());
        std::cout << std::endl;
    }

    return 0;
}
