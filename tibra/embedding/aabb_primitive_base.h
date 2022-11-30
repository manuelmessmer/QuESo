// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef AABB_PRIMITIVE_BASE_INCLUDE_H
#define AABB_PRIMITIVE_BASE_INCLUDE_H

// External includes
#include <cstddef>
#include <array>

#include "containers/point_types.h"

// Global variables
constexpr double kEpsilon = 1e-14;

///@name TIBRA Classes
///@{

/// Forward declaration
class AABB_primitive;

/**
 * @class  AABB_primitive
 * @author Manuel Messmer
 * @brief  Base class for aabb primitives. Derived classes must override intersect().
*/
class AABB_primitive_base
{
public:

    ///@}
    ///@name Operations
    ///@{
    virtual bool intersect(const AABB_primitive &aabb) const = 0;

    ///@}
};
///@} // end classes

#endif // AABB_PRIMITIVE_BASE_INCLUDE_H