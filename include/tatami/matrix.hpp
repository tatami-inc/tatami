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
    virtual workspace_ptr new_workspace(bool row) const { return nullptr; }

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
};

}

#endif
