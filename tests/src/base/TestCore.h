#ifndef TEST_CORE_H 
#define TEST_CORE_H 

#include <gtest/gtest.h>

#include <vector>
#include <algorithm>

template<class BASE>
class TestCore : public BASE {
private:
    static void copy(const double* ptr, std::vector<double>& output) {
        if (ptr != output.data()) {
            std::copy(ptr, ptr + output.size(), output.begin());
        }
        return;
    }

    static void copy(const tatami::sparse_range<double, int>& range, std::vector<double>& output, size_t shift = 0) {
        std::fill(output.begin(), output.end(), 0);
        for (size_t i = 0; i < range.number; ++i) {
            if (i) {
                ASSERT_TRUE(range.index[i] > range.index[i-1]);
            }
            output[range.index[i] - shift] = range.value[i];
        }
        return;
    }

    template<bool ROW, class M> 
    static size_t otherdim (const M* ptr) {
        if constexpr(ROW) {
            return ptr->ncol();
        } else {
            return ptr->nrow();
        }
    }

protected:
    template<bool ROW, class M>
    std::vector<double> extract_dense(const M* ptr, size_t i) {
        std::vector<double> output(otherdim<ROW>(ptr));

        const double* out = NULL;
        if constexpr(ROW) {
            out = ptr->row(i, output.data());
        } else {
            out = ptr->column(i, output.data());
        }

        copy(out, output);
        return output;
    }

    template<bool ROW, class M>
    std::vector<double> extract_dense(const M* ptr, size_t i, tatami::workspace* work) {
        std::vector<double> output(otherdim<ROW>(ptr));

        const double* out = NULL;
        if constexpr(ROW) {
            out = ptr->row(i, output.data(), work);
        } else {
            out = ptr->column(i, output.data(), work);
        }

        copy(out, output);
        return output;
    }

protected:
    template<bool ROW, class M>
    std::vector<double> extract_dense(const M* ptr, size_t i, size_t first, size_t last) {
        std::vector<double> output(last - first);

        const double* out = NULL;
        if constexpr(ROW) {
            out = ptr->row(i, output.data(), first, last);
        } else {
            out = ptr->column(i, output.data(), first, last);
        }

        copy(out, output);
        return output;
    }

    template<bool ROW, class M>
    std::vector<double> extract_dense(const M* ptr, size_t i, size_t first, size_t last, tatami::workspace* work) {
        std::vector<double> output(last - first);

        const double* out = NULL;
        if constexpr(ROW) {
            out = ptr->row(i, output.data(), first, last, work);
        } else {
            out = ptr->column(i, output.data(), first, last, work);
        }

        copy(out, output);
        return output;
    }

protected:
    template<bool ROW, class M>
    std::vector<double> extract_sparse(const M* ptr, size_t i) {
        std::vector<double> output(otherdim<ROW>(ptr));
        std::vector<double> vbuffer(output.size());
        std::vector<int> ibuffer(output.size());

        tatami::sparse_range<double, int> out;
        if constexpr(ROW) {
            out = ptr->sparse_row(i, vbuffer.data(), ibuffer.data());
        } else {
            out = ptr->sparse_column(i, vbuffer.data(), ibuffer.data());
        }

        copy(out, output);
        return output;
    }

    template<bool ROW, class M>
    std::vector<double> extract_sparse(const M* ptr, size_t i, tatami::workspace* work) {
        std::vector<double> output(otherdim<ROW>(ptr));
        std::vector<double> vbuffer(output.size());
        std::vector<int> ibuffer(output.size());

        tatami::sparse_range<double, int> out;
        if constexpr(ROW) {
            out = ptr->sparse_row(i, vbuffer.data(), ibuffer.data(), work);
        } else {
            out = ptr->sparse_column(i, vbuffer.data(), ibuffer.data(), work);
        }

        copy(out, output);
        return output;
    }

protected:
    template<bool ROW, class M>
    std::vector<double> extract_sparse(const M* ptr, size_t i, size_t first, size_t last) {
        std::vector<double> output(last - first);
        std::vector<double> vbuffer(output.size());
        std::vector<int> ibuffer(output.size());

        tatami::sparse_range<double, int> out;
        if constexpr(ROW) {
            out = ptr->sparse_row(i, vbuffer.data(), ibuffer.data(), first, last);
        } else {
            out = ptr->sparse_column(i, vbuffer.data(), ibuffer.data(), first, last);
        }

        copy(out, output, first);
        return output;
    }

    template<bool ROW, class M>
    std::vector<double> extract_sparse(const M* ptr, size_t i, size_t first, size_t last, tatami::workspace* work) {
        std::vector<double> output(last - first);
        std::vector<double> vbuffer(output.size());
        std::vector<int> ibuffer(output.size());

        tatami::sparse_range<double, int> out;
        if constexpr(ROW) {
            out = ptr->sparse_row(i, vbuffer.data(), ibuffer.data(), first, last, work);
        } else {
            out = ptr->sparse_column(i, vbuffer.data(), ibuffer.data(), first, last, work);
        }

        copy(out, output, first);
        return output;
    }

protected:
    static std::pair<size_t, size_t> wrap_intervals(size_t first, size_t last, size_t max) {
        size_t diff = last - first;
        first %= max;
        last = std::min(max, first + diff);
        return std::make_pair(first, last);
    }

};

#endif
