#include "tatami/tatami.h"
#include <vector>
#include <iostream>
#include <numeric>
#include <iomanip>

/* INTRODUCTION: 
 *
 * In this example, we will be computing the column sums of a tatami numeric
 * matrix. The task itself is pretty simple - iterate across all columns,
 * extract the columnar values, compute the sum, and job done. However, there
 * are several ways that this can actually be implemented, depending on the
 * time we are willing to spend and the efficiency that we want.
 */


/* ATTEMPT 1:
 *
 * This is the simplest implementation whereby we just loop across all
 * columns, extract each column and compute the sum. 
 */
std::vector<double> colsums_simple(std::shared_ptr<tatami::numeric_matrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);

    /* Extraction requires a buffer that can be filled up with the column's
     * values. This may or may not be used, depending on the actual
     * implementation of the column() method for a particular tatami
     * subclass.
     */
    std::vector<double> buffer(NR);

    for (size_t i = 0; i < NC; ++i) {
        auto ptr = p->column(i, buffer.data());
        output[i] = std::accumulate(ptr, ptr + NR, 0.0);
    }

    return output;
}

/* ATTEMPT 2:
 *
 * For extraction of consecutive columns, we can supply a "workspace" that
 * shares data/information between columns for more efficient extraction. The
 * exact nature of the sharing (if any) depends on the underlying matrix
 * implementation, e.g., open file handles, pointers to memory locations. This
 * mechanism also makes extraction easy to parallelize as each thread/process
 * gets it own workspace and the underlying matrix is const on read access.
 */
std::vector<double> colsums_work(std::shared_ptr<tatami::numeric_matrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);
    std::vector<double> buffer(NR);

    /* We set 'false' here to indicate that we want our workspace to operate
     * for column extraction, not for row extraction. For example, no workspace
     * is necessary when performing column-wise extraction on a column-major
     * matrix, in which case the returned 'wrk' will be a null pointer.
     */
    auto wrk = p->new_workspace(false);

    for (size_t i = 0; i < NC; ++i) {
        auto ptr = p->column(i, buffer.data(), wrk);
        output[i] = std::accumulate(ptr, ptr + NR, 0.0);
    }

    return output;
}

/* ATTEMPT 3:
 *
 * If our matrix is sparse, we can speed up the entire process by only extracting
 * and computing column sums on the non-zero elements. This requires a separate
 * chunk of code to extract only the non-zero values. Of course, it only makes 
 * sense to write a high-efficiency extraction path if our downstream calculations
 * can take advantage of the sparsity.
 */
std::vector<double> colsums_sparse(std::shared_ptr<tatami::numeric_matrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);
    std::vector<double> buffer(NR);
    auto wrk = p->new_workspace(false);

    /* When performing sparse extractions, we need an additional buffer for the
     * indices - in this case, the row indices of the non-zero elements.
     */
    std::vector<int> ibuffer(NR);

    for (size_t i = 0; i < NC; ++i) {
        const double* ptr = NULL;
        size_t n = NR;

        /* Sparse extractions return a struct containing 'number', the number
         * of non-zero elements; 'value', a pointer to the non-zero values; and
         * 'index', the (row) indices of the non-zero values.
         */
        if (p->sparse()) {
            auto range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk);
            n = range.number;
            ptr = range.value;
        } else {
            ptr = p->column(i, buffer.data(), wrk);
        }

        output[i] = std::accumulate(ptr, ptr + n, 0.0);
    }

    return output;
}

/* ATTEMPT 4:
 *
 * tatami matrices have the capability to report which dimension they prefer to
 * iterate over. For example, row-major matrices will prefer to perform
 * iterations over the rows rather than the columns, as the former is more
 * cache-friendly. In such cases, we can write yet another code path to
 * support row-based extraction and compute column sums by adding each row onto
 * a running sum. Again, this is only possible for certain calculations; if we
 * really need an entire column's values at once, there's really no way around
 * calling column().
 */
std::vector<double> colsums_preferred(std::shared_ptr<tatami::numeric_matrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);

    // Deciding whether or not to perform row-wise extraction.
    if (p->prefer_rows()) {
        std::vector<double> buffer(NC);
        auto wrk = p->new_workspace(true);
        std::vector<int> ibuffer(NC);

        for (size_t i = 0; i < NR; ++i) {
            if (p->sparse()) {
                auto range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk);
                for (size_t j = 0; j < range.number; ++j) {
                    output[range.index[j]] += range.value[j];
                }
            } else {
                auto ptr = p->column(i, buffer.data(), wrk);
                for (size_t j = 0; j < NC; ++j) {
                    output[j] += ptr[j];
                }
            }
        }
    } else {
        std::vector<double> buffer(NR);
        auto wrk = p->new_workspace(false);
        std::vector<int> ibuffer(NR);

        for (size_t i = 0; i < NC; ++i) {
            const double* ptr = NULL;
            size_t n = NR;

            if (p->sparse()) {
                auto range = p->sparse_column(i, buffer.data(), ibuffer.data(), wrk);
                n = range.number;
                ptr = range.value;
            } else {
                ptr = p->column(i, buffer.data(), wrk);
            }

            output[i] = std::accumulate(ptr, ptr + n, 0.0);
        }
    }

    return output;
}

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
    std::vector<int> rows = { 3, 5, 0, 1, 8, 4, 7, 5, 6, 9, 0, 1, 2, 3, 4 };
    std::vector<int> cols = { 0, 1, 3, 2, 0, 4, 1, 2, 4, 0, 1, 3, 0, 3, 2 };
    std::vector<double> vals = { -0.40, 0.14, -0.17, 1.20, 1.20, -1.10, -0.42, 2.10, 0.38, 0.40, -1.10, 0.57, -0.89, 1.60, 0.27 };

    auto indptrs = tatami::compress_sparse_triplets(10, 5, vals, rows, cols);
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
        auto cs = colsums_simple(mat); 
        print_vector(cs.begin(), cs.end());
        std::cout << std::endl;
    }
    
    {
        std::cout << "Workspace:" << std::endl;
        auto cs = colsums_work(mat); 
        print_vector(cs.begin(), cs.end());
        std::cout << std::endl;
    }

    {
        std::cout << "Sparse:" << std::endl;
        auto cs = colsums_sparse(mat); 
        print_vector(cs.begin(), cs.end());
        std::cout << std::endl;
    }

    {
        std::cout << "Preferred:" << std::endl;
        auto cs = colsums_preferred(mat); 
        print_vector(cs.begin(), cs.end());
        std::cout << std::endl;
    }

    return 0;
}
