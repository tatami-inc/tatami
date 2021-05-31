#ifndef TATAMI_WORKSPACE_H
#define TATAMI_WORKSPACE_H

/**
 * @file workspace.hpp
 *
 * Defines the virtual base `workspace` class.
 */

namespace tatami {

/**
 * @brief Virtual workspace class.
 *
 * Instances of this class cannot be constructed; it is only provided as a base class for `create_workspace()` methods to return appropriate objects.
 */
class workspace {
public:
    virtual ~workspace() {}
protected:
    workspace() {}
};

}

#endif
