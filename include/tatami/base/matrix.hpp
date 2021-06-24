#ifndef TATAMI_MATRIX_H
#define TATAMI_MATRIX_H

#include "content_type.hpp"
#include "workspace.hpp"

/**
 * @file matrix.hpp
 *
 * Virtual class for all matrices.
 */

namespace tatami {

/**
 * @brief Virtual class for all matrices.
 */
class matrix {
public:
    virtual ~matrix() {}

    /**
     * @return Number of rows.
     */
    virtual size_t nrow() const = 0;

    /**
     * @return Number of columns.
     */
    virtual size_t ncol() const = 0;

    /**
     * @param row Should a workspace for row extraction be returned?
     *
     * @return A shared pointer to a `workspace` for row or column extraction, or a null pointer if no workspace is required.
     * Defaults to returning a null pointer if no specialized method is provided in derived classes.
     */
    virtual std::shared_ptr<workspace> new_workspace(bool row) const { return nullptr; }

    /**
     * @return A `content_type` specifying the type of the values in the matrix.
     * Defaults to `_unknown` if no specialized method is provided in derived classes.
     */
    virtual content_type type() const { return _unknown; }

    /**
     * @return Is this matrix sparse?
     * Defaults to `false` if no specialized method is provided in derived classes.
     */
    virtual bool sparse() const { return false; }

    /**
     * @return The preferred dimension for extracting values.
     * If `true`, row-wise extraction is preferred; if `false`, column-wise extraction is preferred.
     * Defaults to `false` if no specialized method is provided in derived classes.
     */
    virtual bool prefer_rows() const { return false; }

    /**
     * @return A `pair` containing the number of matrix elements that prefer row-level access (`first`) or column-level access (`second`).
     *
     * This method is useful for determining the return value of `prefer_rows()` in combined matrices consisting of both row- and column-preferred submatrices.
     * In such cases, the net preference can be determined based on the combined size of the submatrices for each preference.
     *
     * For simpler matrices, the return value contains the total size of the matrix in one of the `double`s and zero in the other.
     */
    virtual std::pair<double, double> dimension_preference () const {
        double size = static_cast<double>(nrow()) * static_cast<double>(ncol());
        if (prefer_rows()) {
            return std::make_pair(size, 0.0);
        } else {
            return std::make_pair(0.0, size);
        }
    }
};

}

#endif
