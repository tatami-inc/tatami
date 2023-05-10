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

                    // If there's no non-zero values, we just leave the min/max as zero.
                    if (out.number) {
                        if constexpr(store_min) {
                            auto minned = *std::min_element(out.value, out.value + out.number);
                            if (minned < 0 || out.number == otherdim) {
                                min_out[i] = minned;
                            }
                        }
                        if constexpr(store_max) {
                            auto maxed = *std::max_element(out.value, out.value + out.number);
                            if (maxed > 0 || out.number == otherdim) {
                                max_out[i] = maxed;
                            }
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
                        auto val = static_cast<Output_>(out.value[j]);
                        if constexpr(store_min) {
                            auto& last = min_out[idx];
                            if (i == 0 || last > val) {
                                last = val;
                            }
                        } 
                        if constexpr(store_max) {
                            auto& last = max_out[idx];
                            if (i == 0 || last < val) {
                                last = val;
                            }
                        }
                        ++counter[idx - s];
                    }
                }

                // Handling the zeros.
                for (Index_ i = s, e = s + l; i < e; ++i) {
                    if (counter[i - s] < otherdim) {
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

                for (Index_ i = 0; i < otherdim; ++i) {
                    auto ptr = ext->fetch(i, buffer.data());
                    if (i) {
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
                    } else {
                        if constexpr(store_min) {
                            std::copy(ptr, ptr + len, min_out.data() + s);
                        }
                        if constexpr(store_max) {
                            std::copy(ptr, ptr + len, max_out.data() + s);
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
 * @param threads Number of threads to use.
 *
 * @return A vector of length equal to the number of columns, containing the maximum value in each column.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> column_maxs(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> output(p->ncol());
    bool temp = false;
    stats::dimension_extremes<false, Output_>(p, threads, temp, output);
    return output;
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
    bool temp = false;
    stats::dimension_extremes<true, Output_>(p, threads, temp, output);
    return output;
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
    bool temp = false;
    stats::dimension_extremes<false, Output_>(p, threads, output, temp);
    return output;
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
    bool temp = false;
    stats::dimension_extremes<true, Output_>(p, threads, output, temp);
    return output;
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
std::pair<std::vector<Output_>, std::vector<Output_> > column_ranges(const Matrix<Value_, Index_>* p, int threads = 1) {
    std::vector<Output_> mins(p->ncol()), maxs(p->ncol());
    stats::dimension_extremes<false, Output_>(p, threads, mins, maxs);
    return std::make_pair(std::move(mins), std::move(maxs));
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
    stats::dimension_extremes<true, Output_>(p, threads, mins, maxs);
    return std::make_pair(std::move(mins), std::move(maxs));
}

}

#endif
