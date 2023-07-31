// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef GEOMETRY_QUERY_INCLUDE_H
#define GEOMETRY_QUERY_INCLUDE_H

/// Project includes
#include "containers/triangle_mesh.hpp"
#include "embedding/ray_aabb_primitive.h"
#include "embedding/aabb_primitive.h"
#include "embedding/aabb_tree.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  Polygon
 * @author Manuel Messmer
 * @brief Provides functions to query geometry tests. Uses AABB tree to accelerate this process.
*/
class GeometryQuery {

public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    GeometryQuery(const TriangleMesh& rTriangleMesh, bool MeshIsClosed = true)
        : mTriangleMesh(rTriangleMesh), mTree(rTriangleMesh), mMeshIsClosed(MeshIsClosed)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns true, if rPoint is within the bounding box of the given triangle mesh.
    /// @param rPoint
    /// @return bool.
    bool IsWithinBoundingBox(const PointType& rPoint) const;

    /// @brief Ray tracing to check, if a Point is inside or outside of the given triangle mesh.
    /// @details Calls: IsInsideOpen or IsInsideClosed depending on mMeshIsClosed.
    /// @param rRay Ray.
    /// @return std::pair<bool, bool> first-is_inside second-test_successful.
    std::pair<bool, bool> IsInside( const Ray_AABB_primitive& rRay ) const;

    /// @brief Returns distance to closed triangle.
    /// @param rRay
    /// @return double
    double DistanceToClosestTriangle( const Ray_AABB_primitive& rRay ) const;

    /// @brief Returns true, if the AABB intersects with the given triangle mesh.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param Tolerance Reduces size of AABB.
    /// @return bool
    bool DoIntersect(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const;

    /// @brief Returns a vector of ids of all triangles that intersect with the AABB.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param Tolerance Reduces size of AABB.
    /// @return Unique<std::vector<IndexType>>
    Unique<std::vector<IndexType>> GetIntersectedTriangleIds(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const;

private:
    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Perfoms classical ray tracing. IsInside=true, if the number of intersected triangles is odd.
    /// @param rRay Ray.
    /// @return std::pair<bool, bool> first-is_inside second-test_successful
    std::pair<bool, bool> IsInsideClosed( const Ray_AABB_primitive& rRay ) const;

    /// @brief IsInside=true, if the nearest triangle is back facing.
    /// @param rRay Ray.
    /// @return std::pair<bool, bool> first-is_inside second-test_successful
    std::pair<bool, bool> IsInsideOpen( const Ray_AABB_primitive& rRay ) const;

    ///@}
    ///@name Private Members
    ///@{

    const TriangleMesh& mTriangleMesh;
    AABB_tree mTree;
    bool mMeshIsClosed;
    ///@}

}; // End GeometryQuery class
///@} End TIBRA classes

} // End namespace tibra

#endif // GEOMETRY_QUERY_INCLUDE_H