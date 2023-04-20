// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BREP_OPERATOR_INCLUDE_H
#define BREP_OPERATOR_INCLUDE_H

//// STL includes
#include <memory>
//// Project includes
#include "containers/triangle_mesh.hpp"
#include "containers/element.hpp"
#include "containers/boundary_integration_point.hpp"
#include "embedding/trimmed_domain.h"
#include "embedding/geometry_query.h"
#include "embedding/clipper.h"
#include "embedding/brep_operator_base.h"
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
class BRepOperator : public BRepOperatorBase {

public:
    ///@name Type Definitions
    ///@{

    typedef BRepOperatorBase BaseType;
    typedef BaseType::TrimmedDomainBasePtrType TrimmedDomainBasePtrType;

    // Declare BaseType functions.
    using BaseType::GetIntersectionState;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param rTriangleMesh
    ///@param rParameters TIBRA parameters.
    BRepOperator(const TriangleMesh& rTriangleMesh, const Parameters& rParameters)
        : BaseType(rParameters), mTriangleMesh(rTriangleMesh), mGeometryQuery(rTriangleMesh, true)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    bool IsInside(const PointType& rPoint) const override;

    ///@brief Returns intersections state of element.
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@param Tolerance Tolerance reduces size of element/AABB slightly. Default: SNAPTOL. If Tolerance=0 touch is detected as intersection.
    ///                 If Tolerance>0, touch is not detected as intersection.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance = SNAPTOL) const override;

    /// @brief Returns ptr to trimmed domain. Trimmed domain contains intersection mesh.(see: GetTriangleMesh())
    /// @param rLowerBound Lower bound of AABB.
    /// @param rUpperBound Upper bound of AABB.
    /// @return TrimmedDomainBasePtrType (Unique)
    TrimmedDomainBasePtrType pGetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound ) const override;

    ///@brief Clips triangle mesh by AABB.
    ///       Will NOT keep triangles that are categorized to be on one of the six planes of AABB.
    ///       This is a requirement for the intersection algorithm (see: TrimemdDomain and TrimmedDomainOnPlane).
    ///@see pClipTriangleMeshUnique().
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@return Unique<TriangleMesh>. Clipped mesh.
    Unique<TriangleMesh> pClipTriangleMesh(const PointType& rLowerBound, const PointType& rUpperBound ) const;

    /// @brief Returns true, if AABB is intersected by at least one triangle.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param Tolerance Reduces size of AABB.
    /// @return bool.
    bool IsTrimmed(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance = SNAPTOL) const override;

    /// @brief Returns true if rPoint lies on bounded side of clipped mesh (clipped by AABB).
    ///        Ray tracing trough the center of at least 10 triangles (or maximum number of triangles, if n_max < 10) is performed.
    ///        The majority decides about the classification of rPoint. Note that this function is much more efficient than IsInside.
    ///        However, rPoint must be close to AABB. This is e.g. used to classify an aabb next to a trimmed aabb.
    /// @param rPoint Query Point.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @return bool
    bool OnBoundedSideOfClippedSection( const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound ) const override;

    ///@brief ProtoType: Clips triangle mesh by AABB. This function keeps triangles that are categorized on the planes of AABB.
    ///       However, to avoid that triangles are assigned twice to both adjacent AABB's, they are only assigned to the positive planes (+x, +y, +z).
    ///       This is a requirement for the application of boundary conditions.
    ///@todo This needs improvement. Probably global function that cuts every plane only once, to guarantee that triangles on the planes are not assigned twice.
    ///@see pClipTriangleMesh()
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@return Unique<TriangleMesh>. Clipped mesh.
    Unique<TriangleMesh> pClipTriangleMeshUnique(const PointType& rLowerBound, const PointType& rUpperBound ) const;
    ///@}

private:
    ///@name Private Operations
    ///@{

    ///@brief Return ids of triangles that intersect AABB.
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@param Tolerance Positve tolerance reduces extent of AABB. This can be used to neglect touching triangles.
    ///@return Unique<std::vector<IndexType>> containing ids.
    Unique<std::vector<IndexType>> GetIntersectedTriangleIds( const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance ) const;

    ///@}
    ///@name Private Members
    ///@{

    const TriangleMesh& mTriangleMesh;
    GeometryQuery mGeometryQuery;

    ///@}
}; // End BRepOperator class

///@} // End TIBRA classes

} // End namespace tibra

#endif // BREP_OPERATOR_INCLUDE_H