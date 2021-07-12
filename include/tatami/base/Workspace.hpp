#ifndef TATAMI_WORKSPACE_H
#define TATAMI_WORKSPACE_H

#include <memory>

/**
 * @file Workspace.hpp
 *
 * Defines the virtual base `Workspace` class.
 */

namespace tatami {

/**
 * @brief Virtual workspace class.
 *
 * Instances of this class cannot be constructed; it is only provided as a base class for `matrix::new_workspace()` methods to return appropriate objects.
 */
class Workspace {
public:
    virtual ~Workspace() {}
protected:
    Workspace() {}
};

}

#endif
