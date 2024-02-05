#ifndef TATAMI_STATS_RANGES_HPP
#define TATAMI_STATS_RANGES_HPP

#include "../base/Matrix.hpp"
#include "utils.hpp"
#include <vector>
#include <algorithm>

/**
 * @file ranges.hpp
 *
 * @brief Compute row and column ranges from a `tatami::Matrix`.
 */

namespace tatami {

namespace stats {

/**
 * @cond
 */
template<bool row_, typename Output_, typename Value_, typename Index_, class StoreMinimum_, class StoreMaximum_>
void dimension_extremes(const Matrix<Value_, Index_>* p, int threads, StoreMinimum_& min_out, StoreMaximum_& max_out) {
    auto dim = (row_ ? p->nrow() : p->ncol());
    auto otherdim = (row_ ? p->ncol() : p->nrow());
    const bool direct = p->prefer_rows() == row_;

    constexpr bool store_min = !std::is_same<StoreMinimum_, bool>::value;
    constexpr bool store_max = !std::is_same<StoreMaximum_, bool>::value;

    if (!otherdim) {
        if constexpr(store_min) {
            std::fill(min_out, min_out + dim, 0);
        }
        if constexpr(store_max) {
            std::fill(max_out, max_out + dim, 0);
        }
        return;
    }

    if (p->sparse()) {
        Options opt;
        opt.sparse_ordered_index = false;

        if (direct) {
            opt.sparse_extract_index = false;
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<row_, true>(p, s, l, opt);
                std::vector<Value_> vbuffer(otherdim);

                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), NULL);

                    if (out.number) {
                        if constexpr(store_min) {
                            auto minned = *std::min_element(out.value, out.value + out.number);
                            if (minned > 0 && out.number != otherdim) {
                                minned = 0;
                            }
                            min_out[i] = minned;
                        }
                        if constexpr(store_max) {
                            auto maxed = *std::max_element(out.value, out.value + out.number);
                            if (maxed < 0 && out.number != otherdim) {
                                maxed = 0;
                            }
                            max_out[i] = maxed;
                        }
                    } else {
                        if constexpr(store_min) {
                            min_out[i] = 0;
                        }
                        if constexpr(store_max) {
                            max_out[i] = 0;
                        }
                    }
                }
            }, dim, threads);

        } else {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<!row_, true>(p, 0, otherdim, s, l, opt);
                auto len = ext->block_length;
                std::vector<Value_> vbuffer(len);
                std::vector<Index_> ibuffer(len);
                std::vector<Index_> counter(len);

                for (Index_ i = 0; i < otherdim; ++i) {
                    auto out = ext->fetch(i, vbuffer.data(), ibuffer.data());
                    for (Index_ j = 0; j < out.number; ++j) {
                        auto idx = out.index[j];
                        auto& c = counter[idx - s];
                        auto val = static_cast<Output_>(out.value[j]);
                        if constexpr(store_min) {
                            auto& last = min_out[idx];
                            if (c == 0 || last > val) {
                                last = val;
                            }
                        } 
                        if constexpr(store_max) {
                            auto& last = max_out[idx];
                            if (c == 0 || last < val) {
                                last = val;
                            }
                        }
                        ++c;
                    }
                }

                // Handling the zeros.
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto c = counter[i - s];
                    if (c == 0) {
                        if constexpr(store_min) {
                            min_out[i] = 0;
                        }
                        if constexpr(store_max) {
                            max_out[i] = 0;
                        }
                    } else if (c < otherdim) {
                        if constexpr(store_min) {
                            if (min_out[i] > 0) {
                                min_out[i] = 0;
                            }
                        }
                        if constexpr(store_max) {
                            if (max_out[i] < 0) {
                                max_out[i] = 0;
                            }
                        }
                    }
                }
            }, dim, threads);
        }

    } else {
        if (direct) {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<row_, false>(p, s, l);
                std::vector<Value_> buffer(otherdim);
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    auto ptr = ext->fetch(i, buffer.data());
                    if constexpr(store_min) {
                        min_out[i] = *std::min_element(ptr, ptr + otherdim);
                    }
                    if constexpr(store_max) {
                        max_out[i] = *std::max_element(ptr, ptr + otherdim);
                    } 
                }
            }, dim, threads);

        } else {
            parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = consecutive_extractor<!row_, false>(p, 0, otherdim, s, l);
                auto len = ext->block_length;
                std::vector<Value_> buffer(len);

                // We already have a otherdim > 0 check above.
                auto ptr = ext->fetch(0, buffer.data());
                if constexpr(store_min) {
                    std::copy(ptr, ptr + len, min_out + s);
                }
                if constexpr(store_max) {
                    std::copy(ptr, ptr + len, max_out + s);
                }

                for (Index_ i = 1; i < otherdim; ++i) {
                    auto ptr = ext->fetch(i, buffer.data());
                    for (Index_ d = 0; d < len; ++d) {
                        auto idx = d + s;
                        auto val = static_cast<Output_>(ptr[d]);
                        if constexpr(store_min) {
                            auto& last = min_out[idx];
                            if (last > val) {
                                last = val;
                            }
                        }
                        if constexpr(store_max) {
                            auto& last = max_out[idx];
                            if (last < val) {
                                last = val;
                            }
                        }
                    }
                }
            }, dim, threads);
        }
    }

    return;
}
/**
 * @endcond
 */

}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of columns.
 * On output, this contains the maximum value in each column.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_maxs(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    bool temp = false;
    stats::dimension_extremes<false, Output_>(p, threads, temp, output);
    return;
}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the maximum value in each column.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> column_maxs(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    column_maxs(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows.
 * On output, this contains the maximum value in each row.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void row_maxs(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    bool temp = false;
    stats::dimension_extremes<true, Output_>(p, threads, temp, output);
    return;
}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the maximum value in each row.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> row_maxs(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->nrow());
    row_maxs(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of columns.
 * On output, this contains the minimum value in each column.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_mins(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    bool temp = false;
    stats::dimension_extremes<false, Output_>(p, threads, output, temp);
    return; 
}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the minimum value in each column.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> column_mins(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    column_mins(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows.
 * On output, this is filled with the minimum value in each row.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void row_mins(const Matrix<Value_, Index_>* p, Output_* output, int threads = 1) {
    bool temp = false;
    stats::dimension_extremes<true, Output_>(p, threads, output, temp);
    return;
}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of rows, containing the minimum value in each row.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> row_mins(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->nrow());
    row_mins(p, output.data(), threads);
    return output;
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] min_output Pointer to an array of length equal to the number of rows.
 * On output, this contains the minimum value per row.
 * @param[out] max_output Pointer to an array of length equal to the number of rows.
 * On output, this contains the maximum value per row.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void column_ranges(const Matrix<Value_, Index_>* p, Output_* min_output, Output_* max_output, int threads = 1) {
    stats::dimension_extremes<false, Output_>(p, threads, min_output, max_output);
    return;
}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A pair of vectors, each of length equal to the number of columns.
 * The first and second vector contains the minimum and maximum value per column, respectively.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<std::vector<Output_>, std::vector<Output_> > column_ranges(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> mins(p->ncol()), maxs(p->ncol());
    column_ranges(p, mins.data(), maxs.data(), threads);
    return std::make_pair(std::move(mins), std::move(maxs));
}

/**
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] min_output Pointer to an array of length equal to the number of rows.
 * On output, this contains the minimum value per row.
 * @param[out] max_output Pointer to an array of length equal to the number of rows.
 * On output, this contains the maximum value per row.
 * @param threads Number of threads to use.
 */
template<typename Value_, typename Index_, typename Output_>
void row_ranges(const Matrix<Value_, Index_>* p, Output_* min_output, Output_* max_output, int threads = 1) {
    stats::dimension_extremes<true, Output_>(p, threads, min_output, max_output);
    return;
}

/**
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param threads Number of threads to use.
 *
 * @return A pair of vectors, each of length equal to the number of rows.
 * The first and second vector contains the minimum and maximum value per row, respectively.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<std::vector<Output_>, std::vector<Output_> > row_ranges(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> mins(p->nrow()), maxs(p->nrow());
    row_ranges(p, mins.data(), maxs.data(), threads);
    return std::make_pair(std::move(mins), std::move(maxs));
}

}

#endif
