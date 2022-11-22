// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_BASE_INCLUDE_H
#define BREP_OPERATOR_BASE_INCLUDE_H

/// External includes
#include <memory>

/// Project includes
#include "geometries/element.h"
#include "geometries/triangle_mesh.h"
#include "geometries/boundary_integration_point.h"
#include "utilities/parameters.h"

///@name TIBRA Classes
///@{

/**
 * @class  BRepOperatorBase
 * @author Manuel Messmer
 * @brief  Base class for BRepOperator
*/
class BRepOperatorBase {

public:
    ///@name Type Definitions
    ///@{
    typedef TriangleMesh::Vector3d PointType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
    typedef std::unique_ptr<BoundaryIPVectorType> BoundaryIPVectorPtrType;

    enum IntersectionStatus {Inside, Outside, Trimmed};

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    virtual bool IsInside(const PointType& rPoint) const = 0;

    ///@brief Returns intersections state of element.
    ///@param rElement
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    virtual IntersectionStatus GetIntersectionState(const Element& rElement){
        const auto& lower_bound = rElement.GetGlobalLowerPoint();
        const auto& upper_bound = rElement.GetGlobalUpperPoint();
        return GetIntersectionState(lower_bound, upper_bound);
    }

    ///@brief Returns intersections state of element.
    ///@param rLowerBound of AABB.
    ///@param rUpperBound of AABB.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    virtual IntersectionStatus GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance=1e-8) const = 0;

    ///@brief Returns boundary integration points of element.
    ///@param rElement
    ///@param [out] rpBoundaryIps
    ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
    virtual bool ComputeBoundaryIps(Element& rElement, BoundaryIPVectorPtrType& rpBoundaryIps, const Parameters& rParam) const = 0;

};

#endif