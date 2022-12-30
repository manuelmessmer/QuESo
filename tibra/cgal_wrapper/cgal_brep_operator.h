// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CGAL_BREP_OPERATOR_INCLUDE_H
#define CGAL_BREP_OPERATOR_INCLUDE_H

//// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

//// Project includes
#include "cgal_wrapper/cgal_utilities.h"
#include "containers/triangle_mesh.h"
#include "containers/element.h"
#include "containers/integration_point.h"
#include "embedding/brep_operator_base.h"
#include "utilities/parameters.h"
#include "io/io_utilities.h"

namespace tibra {
namespace cgal {

///@name TIBRA Classes
///@{

/**
 * @class  CGAL brep operator
 * @author Manuel Messmer
 * @brief Provides geometrical operations for Brep models using CGAL (optional dependency of TIBRA).
*/
class CGALBRepOperator : public BRepOperatorBase {

public:
    ///@name Type Definitions
    ///@{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernalType;
    typedef CGALKernalType::Point_3 CGALPointType;
    typedef CGAL::Surface_mesh<CGALPointType> CGALMeshType;
    typedef CGAL::Side_of_triangle_mesh<CGALMeshType, CGALKernalType> CGALInsideTestType;
    typedef CGAL::AABB_face_graph_triangle_primitive<CGALMeshType> CGALAABBFaceGraphsPrimitivesType;
    typedef CGAL::AABB_traits<CGALKernalType, CGALAABBFaceGraphsPrimitivesType> CGALAABBFaceGraphsTraits;
    typedef CGAL::AABB_tree<CGALAABBFaceGraphsTraits> CGALAABBTreeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@param rTriangleMesh
    CGALBRepOperator(const TriangleMesh& rTriangleMesh, const Parameters& rParameters ) : BRepOperatorBase(rParameters)
    {
        // Copy triangle mesh to CGAL mesh.
        cgal::CGALUtilities::CopyMesh(rTriangleMesh, mCGALMesh);
        // Initialize inside test and AABB tree.
        mpCGALInsideTest = std::make_unique<CGALInsideTestType>(mCGALMesh);
        CGAL::Polygon_mesh_processing::build_AABB_tree(mCGALMesh, mCGALAABBTree);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns intersections state of AABB.
    ///@param rLowerBound
    ///@param rUpperBound
    ///@param Tolerance
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance) const override;

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    bool IsInside(const PointType& rPoint) const override;

    /// @brief Returns ptr to trimmed domain.
    /// @param rLowerBound of AABB.
    /// @param rUpperBound of AABB.
    /// @param rParam Parameters
    /// @return TrimmedDomainBasePtrType (std::unique_ptr)
    TrimmedDomainBasePtrType GetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound) const override;

    ///@}
private:

    ///@name Private member variables
    ///@{
    CGALMeshType mCGALMesh;
    /// CGALInsideTestType must be ptr, since it does not have default constructor.
    std::unique_ptr<CGALInsideTestType> mpCGALInsideTest;
    CGALAABBTreeType mCGALAABBTree;
    ///@}

}; // End Class CGALBRepOperator
///@} // End TIBRA classes

} // End namespace cgal
} // End namespace tibra

#endif // CGAL_BREP_OPERATOR_INCLUDE_H