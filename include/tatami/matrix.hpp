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
     * @return A pointer to a workspace for row or column extraction, or a null pointer if no workspace is required.
     * Defaults to returning a null pointer if no specialized method is provided in derived classes.
     */
    virtual workspace* create_workspace(bool row) const { return NULL; }

    /**
     * @return A `content_type` specifying the type of the values in the matrix.
     * Defaults to `_unknown` if no specialized method is provided in derived classes.
     */
    virtual content_type type() const { return _unknown; }

    /**
     * @return Is this matrix sparse?
     * Defaults to `false` if no specialized method is provided in derived classes.
     */
    virtual bool is_sparse() const { return false; }

    /**
     * @return The preferred dimension for extracting values.
     * If 0, row-wise extraction is preferred; if 1, column-wise extraction is preferred.
     * Defaults to 1 if no specialized method is provided in derived classes.
     */
    virtual int preferred_dimension() const { return 1; }
};

}

#endif
