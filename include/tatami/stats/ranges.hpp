#ifndef TATAMI_STATS_RANGES_HPP
#define TATAMI_STATS_RANGES_HPP

#include "../base/Matrix.hpp"
#include "apply.hpp"
#include <vector>
#include <algorithm>

/**
 * @file ranges.hpp
 *
 * Compute row and column ranges from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<typename O, bool compute_max>
struct ExtremeFactory {
public:
    ExtremeFactory(O* o, size_t d1, size_t d2) : output(o), dim(d1), otherdim(d2) {}

private:
    O* output;
    size_t dim, otherdim;

public:
    struct DenseDirect {
        DenseDirect(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename V>
        void compute(size_t i, const V* ptr) {
            if (otherdim) {
                if constexpr(compute_max) {
                    output[i] = *std::max_element(ptr, ptr + otherdim);
                } else {
                    output[i] = *std::min_element(ptr, ptr + otherdim);
                }
            }
        }
    private:
        O* output;
        size_t otherdim;
    };

    DenseDirect dense_direct() {
        return DenseDirect(output, otherdim);
    }

public:
    struct SparseDirect {
        SparseDirect(O* o, size_t d2) : output(o), otherdim(d2) {}

        template<typename T, typename IDX>
        void compute(size_t i, const SparseRange<T, IDX>& range) {
            if (range.number) {
                if constexpr(compute_max) {
                    output[i] = *std::max_element(range.value, range.value + range.number);
                } else {
                    output[i] = *std::min_element(range.value, range.value + range.number);
                }

                if (range.number != otherdim) {
                    if constexpr(compute_max) {
                        if (output[i] < 0) {
                            output[i] = 0;
                        }
                    } else {
                        if (output[i] > 0) {
                            output[i] = 0;
                        }
                    }
                }
            } else if (otherdim) {
                output[i] = 0;
            }
        }
    private:
        O* output;
        size_t otherdim;
    };

    SparseDirect sparse_direct() {
        return SparseDirect(output, otherdim);
    }

public:
    struct DenseRunning {
        DenseRunning(O* o, size_t d1) : output(o), dim(d1) {}

        template<typename V>
        void add(const V* ptr) {
            if (first) {
                std::copy(ptr, ptr + dim, output);
                first = false;
            } else {
                for (size_t d = 0; d < dim; ++d) {
                    if constexpr(compute_max) {
                        if (output[d] < ptr[d]) {
                            output[d] = ptr[d];
                        }
                    } else {
                        if (output[d] > ptr[d]) {
                            output[d] = ptr[d];
                        }
                    }
                }
            }
        }

        void finish() {}
    private:
        O* output;
        size_t dim;
        bool first = true;
    };

    DenseRunning dense_running() {
        return DenseRunning(output, dim);
    }

    DenseRunning dense_running(size_t start, size_t end) {
        return DenseRunning(output + start, end - start);
    }

public:
    struct SparseRunning {
        SparseRunning(O* o, size_t d1, size_t d2, size_t s, size_t e) : output(o), collected(d1), otherdim(d2), start(s), end(e) {}

        template<typename T = double, typename IDX = int>
        void add(const SparseRange<T, IDX>& range) {
            if (first) {
                // Assume output is zero-initialized.
                for (size_t j = 0; j < range.number; ++j) {
                    ++collected[range.index[j]];
                    output[range.index[j]] = range.value[j];
                }
                first = false;
            } else {
                for (size_t j = 0; j < range.number; ++j) {
                    ++collected[range.index[j]];
                    auto& existing = output[range.index[j]];
                    if constexpr(compute_max) {
                        if (existing < range.value[j]) {
                            existing = range.value[j];
                        }
                    } else {
                        if (existing > range.value[j]) {
                            existing = range.value[j];
                        }
                    }
                }
            }
        }

        void finish() {
            for (size_t i = start; i < end; ++i) {
                if (collected[i] < otherdim) {
                    if constexpr(compute_max) {
                        if (output[i] < 0) {
                            output[i] = 0;
                        }
                    } else {
                        if (output[i] > 0) {
                            output[i] = 0;
                        }
                    }
                }
            }
        }
    private:
        O* output;
        bool first = true;
        std::vector<size_t> collected;
        size_t otherdim;
        size_t start, end;
    };

    SparseRunning sparse_running() {
        return SparseRunning(output, dim, otherdim, 0, dim);
    }

    SparseRunning sparse_running(size_t start, size_t end) {
        return SparseRunning(output, dim, otherdim, start, end);
    }
};

template<typename O>
using MaxFactory = ExtremeFactory<O, true>;

template<typename O>
using MinFactory = ExtremeFactory<O, false>;

/**
 * @endcond
 */

}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the maximum value in each column.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> column_maxs(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->ncol());
    stats::MaxFactory<Output> factory(output.data(), p->ncol(), p->nrow());
    apply<1>(p, factory);
    return output;
}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the maximum value in each row.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> row_maxs(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->nrow());
    stats::MaxFactory<Output> factory(output.data(), p->nrow(), p->ncol());
    apply<0>(p, factory);
    return output;
}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the minimum value in each column.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> column_mins(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->ncol());
    stats::MinFactory<Output> factory(output.data(), p->ncol(), p->nrow());
    apply<1>(p, factory);
    return output;
}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the minimum value in each row.
 */
template<typename Output = double, typename T, typename IDX>
std::vector<Output> row_mins(const Matrix<T, IDX>* p) {
    std::vector<Output> output(p->nrow());
    stats::MinFactory<Output> factory(output.data(), p->nrow(), p->ncol());
    apply<0>(p, factory);
    return output;
}

namespace stats {

/**
 * @cond
 */
template<typename O>
struct RangeFactory {
public:
    RangeFactory(O* min, O* max, size_t d1, size_t d2) : mins(min, d1, d2), maxs(max, d1, d2) {}

private:
    MinFactory<O> mins;
    MaxFactory<O> maxs;

public:
    struct DenseDirect {
        DenseDirect(typename MinFactory<O>::DenseDirect mn, typename MaxFactory<O>::DenseDirect mx) : mins(std::move(mn)), maxs(std::move(mx)) {}

        template<typename V>
        void compute(size_t i, const V* ptr) {
            mins.compute(i, ptr);
            maxs.compute(i, ptr);
            return;
        }
    private:
        typename MinFactory<O>::DenseDirect mins;
        typename MaxFactory<O>::DenseDirect maxs;
    };

    DenseDirect dense_direct() {
        return DenseDirect(mins.dense_direct(), maxs.dense_direct());
    }

public:
    struct SparseDirect {
        SparseDirect(typename MinFactory<O>::SparseDirect mn, typename MaxFactory<O>::SparseDirect mx) : mins(std::move(mn)), maxs(std::move(mx)) {}

        template<typename T, typename IDX>
        void compute(size_t i, const SparseRange<T, IDX>& range) {
            mins.compute(i, range);
            maxs.compute(i, range);
            return;
        }
    private:
        typename MinFactory<O>::SparseDirect mins;
        typename MaxFactory<O>::SparseDirect maxs;
    };

    SparseDirect sparse_direct() {
        return SparseDirect(mins.sparse_direct(), maxs.sparse_direct());
    }

public:
    struct DenseRunning {
        DenseRunning(typename MinFactory<O>::DenseRunning mn, typename MaxFactory<O>::DenseRunning mx) : mins(std::move(mn)), maxs(std::move(mx)) {}

        template<typename V>
        void add(const V* ptr) {
            mins.add(ptr);
            maxs.add(ptr);
            return;
        }

        void finish() {};
    private:
        typename MinFactory<O>::DenseRunning mins;
        typename MaxFactory<O>::DenseRunning maxs;
    };

    DenseRunning dense_running() {
        return DenseRunning(mins.dense_running(), maxs.dense_running());
    }

    DenseRunning dense_running(size_t start, size_t end) {
        return DenseRunning(mins.dense_running(start, end), maxs.dense_running(start, end));
    }

public:
    struct SparseRunning {
        SparseRunning(typename MinFactory<O>::SparseRunning mn, typename MaxFactory<O>::SparseRunning mx) : mins(std::move(mn)), maxs(std::move(mx)) {}

        template<typename T = double, typename IDX = int>
        void add(const SparseRange<T, IDX>& range) {
            mins.add(range);
            maxs.add(range);
            return;
        }

        void finish() {
            mins.finish();
            maxs.finish();
            return;
        };
    private:
        typename MinFactory<O>::SparseRunning mins;
        typename MaxFactory<O>::SparseRunning maxs;
    };

    SparseRunning sparse_running() {
        return SparseRunning(mins.sparse_running(), maxs.sparse_running());
    }

    SparseRunning sparse_running(size_t start, size_t end) {
        return SparseRunning(mins.sparse_running(start, end), maxs.sparse_running(start, end));
    }
};
/**
 * @endcond
 */

}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A pair of vectors, each of length equal to the number of rows.
 * The first and second vector contains the minimum and maximum value per row, respectively.
 */
template<typename Output = double, typename T, typename IDX>
std::pair<std::vector<Output>, std::vector<Output> > column_ranges(const Matrix<T, IDX>* p) {
    std::vector<Output> mins(p->ncol()), maxs(p->ncol());
    stats::RangeFactory factory(mins.data(), maxs.data(), p->ncol(), p->nrow());
    apply<1>(p, factory);
    return std::make_pair(std::move(mins), std::move(maxs));
}

/**
 * @tparam Output Type of the output value.
 * @tparam T Type of the matrix value.
 * @tparam IDX Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A pair of vectors, each of length equal to the number of rows.
 * The first and second vector contains the minimum and maximum value per row, respectively.
 */
template<typename Output = double, typename T, typename IDX>
std::pair<std::vector<Output>, std::vector<Output> > row_ranges(const Matrix<T, IDX>* p) {
    std::vector<Output> mins(p->nrow()), maxs(p->nrow());
    stats::RangeFactory factory(mins.data(), maxs.data(), p->nrow(), p->ncol());
    apply<0>(p, factory);
    return std::make_pair(std::move(mins), std::move(maxs));
}

}

#endif
