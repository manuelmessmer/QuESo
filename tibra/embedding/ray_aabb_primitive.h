
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef RAY_AABB_PRIMITIVE_INCLUDE_H
#define RAY_AABB_PRIMITIVE_INCLUDE_H

/// External includes

/// Project includes
#include "embedding/aabb_primitive.h"

///@name TIBRA Classes
///@{

/**
 * @class  Ray
 * @author Manuel Messmer
 * @brief  Ray to be used in AABB_tree. Direction of ray must be positive oriented (x>0, y>0, z>0).
 *         Check if orientation is actually positive is omitted for better performance.
*/
class Ray_AABB_primitive : public AABB_primitive_base {

public:
    ///@name Type Definitions
    ///@{
    typedef AABB_primitive_base::IndexType IndexType;
    typedef AABB_primitive_base::Vector3d Vector3d;
    typedef AABB_primitive_base::Vector3i Vector3i;

    ///@}
    ///@name Life cycle
    ///@{

    /// Constructor.
    Ray_AABB_primitive(const Vector3d& Origin, const Vector3d& Direction) :
        mOrigin(Origin), mDirection(Direction)
    {
        mInvDirection[0] = 1.0 / Direction[0];
        mInvDirection[1] = 1.0 / Direction[1];
        mInvDirection[2] = 1.0 / Direction[2];
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if origin of ray is inside aabb.
    ///@param aabb
    ///@return bool
    bool inside(const AABB_primitive &aabb) const;

    ///@brief Returns true if ray intersects aabb.
    ///@param aabb
    ///@return bool
    bool intersect(const AABB_primitive &aabb) const override;

    ///@brief Returns true if ray intersects triangle (only checks intersections in positive direction).
    ///@param v0 Triangle Vertex 1
    ///@param v1 Triangle Vertex 2
    ///@param v2 Triangle Vertex 3
    ///@param t Distance to intersection.
    ///@param u Parametric coordinate 1.
    ///@param v Parametric coordinate 2.
    ///@return bool
    bool intersect( const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                    double &t, double &u, double &v) const;

    ///@}
private:

    ///@name Private Members
    ///@{
    Vector3d mOrigin;
    Vector3d mDirection;
    Vector3d mInvDirection;
    Vector3i mSign;
    ///@}
};
///@} // End class Ray

#endif