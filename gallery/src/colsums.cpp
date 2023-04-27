#include "tatami/tatami.hpp"
#include "print_vector.h"
#include <vector>

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
 * This is the simplest implementation whereby we just loop across all columns,
 * extract each column and compute the sum. This creates a "extractor" that
 * shares data/information between columns for more efficient extraction. The
 * exact nature of the sharing (if any) depends on the underlying matrix
 * implementation, e.g., open file handles, pointers to memory locations. This
 * mechanism also makes extraction easy to parallelize as each thread/process
 * gets it own extractor and the underlying matrix is const on read access.
 */
std::vector<double> colsums_simple(std::shared_ptr<tatami::NumericMatrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);

    /* Extraction requires a buffer that can be filled up with the column's
     * values. This may or may not be used, depending on the actual
     * implementation of the column() method for a particular tatami
     * subclass.
     */
    std::vector<double> buffer(NR);

    /* Constructing the extractor for column-wise extraction. Other extractors
     * can be constructed that only extract a slice from each column.
     */
    auto wrk = p->dense_column();

    for (size_t i = 0; i < NC; ++i) {
        auto ptr = wrk->fetch(i, buffer.data());
        output[i] = std::accumulate(ptr, ptr + NR, 0.0);
    }

    return output;
}

/* ATTEMPT 2:
 *
 * If our matrix is sparse, we can speed up the entire process by only extracting
 * and computing column sums on the non-zero elements. This requires a separate
 * chunk of code to extract only the non-zero values. Of course, it only makes 
 * sense to write a high-efficiency extraction path if our downstream calculations
 * can take advantage of the sparsity.
 */
std::vector<double> colsums_sparse(std::shared_ptr<tatami::NumericMatrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);
    std::vector<double> buffer(NR);

    if (p->sparse()) {
        /* When performing sparse extractions, we usually need an additional
         * buffer for the indices of the non-zero elements. However, if we don't
         * care about the indices, we can skip their extraction for efficiency.
         */
        tatami::Options opt;
        opt.sparse_extract_index = false;
        auto wrk = p->sparse_column(opt);

        for (size_t i = 0; i < NC; ++i) {
            /* Sparse extractions return a struct containing 'number', the number
             * of non-zero elements; 'value', a pointer to the non-zero values; and
             * 'index', the (row) indices of the non-zero values. In this case,
             * 'index' will be set to NULL as we don't need the indices for summation.
             */
            auto range = wrk->fetch(i, buffer.data(), NULL);
            output[i] = std::accumulate(range.value, range.value + range.number, 0.0);
        }
    } else {
        // Copied from colsums_simple.
        auto wrk = p->dense_column();
        for (size_t i = 0; i < NC; ++i) {
            auto ptr = wrk->fetch(i, buffer.data());
            output[i] = std::accumulate(ptr, ptr + NR, 0.0);
        }
    }

    return output;
}

/* ATTEMPT 4:
 *
 * tatami matrices can report which dimension they prefer to iterate over. For
 * example, row-major matrices will prefer to perform iterations over the rows
 * rather than the columns, as the former is more cache-friendly. In such
 * cases, we can write yet another code path to support row-based extraction
 * and compute column sums by adding each row onto a running sum. Again, this
 * is only possible for certain calculations; if we really need an entire
 * column's values at once, there's really no way around calling column().
 *
 * Incidentally, this is the approach that tatami::column_sums() uses.
 */
std::vector<double> colsums_preferred(std::shared_ptr<tatami::NumericMatrix> p) {
    size_t NR = p->nrow(), NC = p->ncol();
    std::vector<double> output(NC);

    // Deciding whether or not to perform row-wise extraction.
    if (p->prefer_rows()) {
        std::vector<double> buffer(NC);
        if (p->sparse()) {
            // This time, we actually do need the sparse indices.
            auto wrk = p->sparse_row();
            std::vector<int> ibuffer(NC);

            for (size_t i = 0; i < NR; ++i) {
                auto range = wrk->fetch(i, buffer.data(), ibuffer.data());
                for (size_t j = 0; j < range.number; ++j) {
                    output[range.index[j]] += range.value[j];
                }
            }
        } else {
            auto wrk = p->dense_row();
            for (size_t i = 0; i < NR; ++i) {
                auto ptr = wrk->fetch(i, buffer.data());
                for (size_t j = 0; j < NC; ++j) {
                    output[j] += ptr[j];
                }
            }
        }
    } else {
        // Copied from colsums_sparse.
        std::vector<double> buffer(NR);
        if (p->sparse()) {
            std::vector<int> ibuffer(NR);
            tatami::Options opt;
            opt.sparse_extract_index = false;
            auto wrk = p->sparse_column(opt);

            for (size_t i = 0; i < NC; ++i) {
                auto range = wrk->fetch(i, buffer.data(), ibuffer.data());
                output[i] = std::accumulate(range.value, range.value + range.number, 0.0);
            }
        } else {
            auto wrk = p->dense_column();
            for (size_t i = 0; i < NC; ++i) {
                auto ptr = wrk->fetch(i, buffer.data());
                output[i] = std::accumulate(ptr, ptr + NR, 0.0);
            }
        }
    }

    return output;
}

int main(int argc, char** argv) {
    std::vector<int> rows = { 3, 5, 0, 1, 8, 4, 7, 5, 6, 9, 0, 1, 2, 3, 4 };
    std::vector<int> cols = { 0, 1, 3, 2, 0, 4, 1, 2, 4, 0, 1, 3, 0, 3, 2 };
    std::vector<double> vals = { -0.40, 0.14, -0.17, 1.20, 1.20, -1.10, -0.42, 2.10, 0.38, 0.40, -1.10, 0.57, -0.89, 1.60, 0.27 };

    auto indptrs = tatami::compress_sparse_triplets<false>(10, 5, vals, rows, cols);
    std::shared_ptr<tatami::NumericMatrix> mat(new tatami::CompressedSparseColumnMatrix<double, int>(10, 5, vals, rows, indptrs));

    std::cout << "Matrix preview: " << std::endl;
    std::vector<double> buffer(mat->ncol());
    auto wrk = mat->dense_row();
    for (size_t i = 0; i < mat->nrow(); ++i) {
        auto ptr = wrk->fetch(i, buffer.data());
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
