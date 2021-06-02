#ifndef TATAMI_WORKSPACE_H
#define TATAMI_WORKSPACE_H

#include <memory>

/**
 * @file workspace.hpp
 *
 * Defines the virtual base `workspace` class.
 */

namespace tatami {

/**
 * @brief Virtual workspace class.
 *
 * Instances of this class cannot be constructed; it is only provided as a base class for `matrix::new_workspace()` methods to return appropriate objects.
 */
class workspace {
public:
    virtual ~workspace() {}
protected:
    workspace() {}
};

/**
 * Shared pointer to a `workspace` instance.
 */
typedef std::shared_ptr<workspace> workspace_ptr; 

}

#endif
