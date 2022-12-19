// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_INCLUDE_H
#define BREP_OPERATOR_INCLUDE_H

//// STL includes
#include <memory>
//// Project includes
#include "containers/triangle_mesh.h"
#include "containers/element.h"
#include "containers/boundary_integration_point.h"
#include "embedding/trimmed_domain.h"
#include "embedding/aabb_tree.h"
#include "embedding/clipper.h"
#include "io/io_utilities.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  BRepOperator
 * @author Manuel Messmer
 * @brief  Provides geometrical operations for Brep models.
 * @details Uses AABB Tree for fast search.
*/
class BRepOperator {

public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
    typedef std::unique_ptr<BoundaryIPVectorType> BoundaryIPVectorPtrType;
    typedef std::unique_ptr<TrimmedDomainBase> TrimmedDomainBasePtrType;

    enum IntersectionStatus {Inside, Outside, Trimmed};
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param rTriangleMesh
    BRepOperator(const TriangleMesh& rTriangleMesh, const Parameters& rParameters)
        : mTriangleMesh(rTriangleMesh), mTree(rTriangleMesh), mParameters(rParameters)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    bool IsInside(const PointType& rPoint) const;

    ///@brief Returns intersections state of element.
    ///@param rElement
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const Element& rElement) const;

    ///@brief Returns intersections state of element.
    ///@param rLowerBound
    ///@param rUpperBound
    ///@param Tolerance Tolerance reduces element slightly. If Tolerance=0 touch is detected as intersection.
    ///                 If Tolerance>0, touch is not detected as intersection.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance=EPS1) const;

    /// @brief Returns ptr to trimmed domain. Trimmed domain holds cipped mesh. (not closed).
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param rParam Parameters
    /// @return TrimmedDomainBasePtrType (std::unique_ptr)
    TrimmedDomainBasePtrType GetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound ) const;

    ///@brief Return ids of triangles that intersect AABB.
    ///@param rLowerBound of AABB.
    ///@param rUpperBound of AABB.
    ///@return std::unique_ptr<std::vector<IndexType>> containing ids.
    std::unique_ptr<std::vector<IndexType>> GetIntersectedTriangleIds( const PointType& rLowerBound, const PointType& rUpperBound ) const;

    ///@brief Clips triangle mesh by AABB.
    ///@param rLowerBound of AABB.
    ///@param rUpperBound of AABB.
    ///@return std::unique_ptr<TriangleMesh>. Clipped mesh.
    std::unique_ptr<TriangleMesh> ClipTriangleMesh(const PointType& rLowerBound, const PointType& rUpperBound) const;

    ///@}

private:
    ///@name Private Operations
    ///@{

    ///@brief Returns true if point is inside AABB.
    ///@param rPoint Query point.
    ///@param rLowerBound of AABB.
    ///@param rUpperBound of AABB.
    ///@return bool
    inline bool IsContained(const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound) const{
        if(    rPoint[0] < rLowerBound[0]
            || rPoint[0] > rUpperBound[0]
            || rPoint[1] < rLowerBound[1]
            || rPoint[1] > rUpperBound[1]
            || rPoint[2] < rLowerBound[2]
            || rPoint[2] > rUpperBound[2] )
        {
            return false;
        }

        return true;
    }

    ///@}
    ///@name Private Members
    ///@{

    AABB_tree mTree;
    const TriangleMesh& mTriangleMesh;
    const Parameters& mParameters;
    ///@}
}; // End BRepOperator class

///@} // End TIBRA classes

} // End namespace tibra

#endif // BREP_OPERATOR_INCLUDE_H