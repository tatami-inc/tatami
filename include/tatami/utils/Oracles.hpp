#ifndef TATAMI_ORACLES_HPP
#define TATAMI_ORACLES_HPP

#include "../base/Options.hpp"
#include <numeric>

/**
 * @file Oracles.hpp
 *
 * @brief Predict future accesses during iteration.
 */

namespace tatami {

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses from a known sequence.
 */
template<typename Index_>
struct FixedOracle : public Oracle<Index_> {
    /**
     * @param r Pointer to a constant array of indices.
     * The underlying array should be valid for the lifetime of this `FixedOracle` instance.
     * @param n Length of the array at `r`.
     */
    FixedOracle(const Index_* r, size_t n) : reference(r), length(n) {}

    size_t predict(Index_* predicted, size_t number) {
        auto upto = counter + number;
        auto out = reference + counter;
        if (upto >= length) {
            number = length - counter;
            upto = length;
        }
        counter = upto;
        std::copy(out, out + number, predicted);
        return number;
    }
private:
    const Index_* reference;

    size_t length;

    size_t counter = 0;
};

/**
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @brief Predict future accesses of a consecutive sequence.
 */
template<typename Index_>
struct ConsecutiveOracle : public Oracle<Index_> {
    /**
     * @param s Start index of the consecutive sequence.
     * @param l Length of the sequence.
     */
    ConsecutiveOracle(Index_ s, Index_ l) : end(s + l), counter(s) {}

    size_t predict(Index_* buffer, size_t number) {
        auto upto = counter + number;
        if (upto >= end) {
            number = end - counter;
            upto = end;
        }
        std::iota(buffer, buffer + number, counter);
        counter = upto;
        return number;
    }
private:
    size_t end;

    size_t counter = 0;
};

/**
 * @brief Stream predictions from the oracle.
 *
 * This provides a friendlier interface to the `Oracle` classes that manages the prediction buffer.
 * Consumers of `Oracle` predictions can receive and process one prediction at a time via the `next()` method,
 * without worrying about the allocation and consumption of the buffer filled by `Oracle::predict()`.
 *
 * The `OracleStream` is primarily intended for `Matrix` developers to use within subclass implementations that can exploit oracular predictions.
 * Users should not have to construct `OracleStream` instances.
 * 
 * @tparam Index_ Integer type of the row/column indices.
 */
template<typename Index_>
struct OracleStream {
    /**
     * Default constructor.
     */
    OracleStream() = default;

    /**
     * @param o Pointer to an `Oracle` instance.
     */
    OracleStream(std::unique_ptr<Oracle<Index_> > o) : oracle(std::move(o)) {}

    bool active() const {
        return oracle.get() != nullptr;
    }

    /**
     * Replace the existing `Oracle`.
     * This assumes that existing `Oracle`'s predictions have been completely consumed by `next()`;
     * any unused predictions will be discarded.
     *
     * @param o Pointer to an `Oracle` instance.
     */
    void set(std::unique_ptr<Oracle<Index_> > o) {
        predictions.clear();
        used = 0;
        oracle = std::move(o);
    }

    /**
     * @param[out] prediction The next prediction, set on output.
     *
     * @return Whether the stream of predictions has been exhausted.
     * If `false`, the value of `prediction` may be used on output;
     * otherwise, it should be ignored.
     */
    bool next(Index_& prediction) {
        if (used == predictions.size()) {
            predictions.resize(100); // take up to 100 predictions at a time.
            used = 0;

            auto filled = oracle->predict(predictions.data(), predictions.size());
            predictions.resize(filled);
            if (filled == 0) {
                return false;
            }
        }

        prediction = predictions[used];
        ++used;
        return true;
    }

    /**
     * Move the state backwards by one prediction.
     * This is guaranteed to work at least once after calling `next()`;
     * multiple calls may or may not work, depending on whether the previous state has already been discarded.
     * It is most useful for loops that need to inspect the next prediction before deciding whether to use it.
     *
     * @return Whether or not the state was moved back successfully.
     * This is guaranteed to be `true` once after a call to `next()`.
     */
    bool back() {
        // Easy to see that, in next(), any successful call leaves us with used > 0.
        // So at least one invocation of unwind() will be successful.
        if (!used) {
            return false;
        }
        --used;
        return true; 
    }

private:
    std::unique_ptr<Oracle<Index_> > oracle;
    std::vector<Index_> predictions;
    size_t used = 0;
};

}

#endif
