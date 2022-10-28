// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CGAL_BREP_OPERATOR_INCLUDE_H
#define CGAL_BREP_OPERATOR_INCLUDE_H

// External includes
#include <map>

/// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

// Project includes
#include "geometries/triangle_mesh.h"
#include "embedding/brep_operator_base.h"
#include "geometries/element.h"
#include "geometries/integration_point.h"
#include "utilities/parameters.h"
#include "io/io_utilities.h"

namespace cgal {

///@name TIBRA Classes
///@{

/**
 * @class  CGAL brep operator
 * @author Manuel Messmer
 * @brief Provides geometrical operations for Brep models using cgal.
*/
class BRepOperator : public BRepOperatorBase {

public:
    ///@name Type Definitions
    ///@{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernalType;
    typedef CGALKernalType::Point_3 CGALPointType;
    typedef CGAL::Surface_mesh<CGALPointType> CGALMeshType;
    typedef std::unique_ptr<CGALMeshType> CGALMeshPtrType;
    typedef CGAL::Side_of_triangle_mesh<CGALMeshType, CGALKernalType> CGALInsideTestType;
    typedef std::array<double,3> PointType;

    typedef CGAL::AABB_face_graph_triangle_primitive<CGALMeshType> CGALAABBFaceGraphsPrimitivesType;
    typedef CGAL::AABB_traits<CGALKernalType, CGALAABBFaceGraphsPrimitivesType> CGALAABBFaceGraphsTraits;
    typedef CGAL::AABB_tree<CGALAABBFaceGraphsTraits> CGALAABBTreeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param rTriangleMesh
    BRepOperator(const TriangleMesh& rTriangleMesh )
    {
        // Copy triangle mesh to CGAL mesh.
        std::map<IndexType, CGALMeshType::Vertex_index> index_map{};
        const auto& r_vertices = rTriangleMesh.GetVertices();
        const auto v_it_begin = r_vertices.begin();
        for( IndexType i = 0; i < rTriangleMesh.NumOfVertices(); ++i){
            auto v = *(v_it_begin + i);
            auto index1 = mCGALMesh.add_vertex(CGALPointType(v[0], v[1], v[2]));
            index_map.insert( std::pair<IndexType, CGALMeshType::Vertex_index>( i, index1) );
        }

        for( IndexType i = 0; i < rTriangleMesh.NumOfTriangles(); ++i){
            const auto& ids = rTriangleMesh.VertexIds(i);
            mCGALMesh.add_face( index_map[ids[0]], index_map[ids[1]], index_map[ids[2]]  );
        }
        mpCGALInsideTest = std::make_unique<CGALInsideTestType>(mCGALMesh);
        CGAL::Polygon_mesh_processing::build_AABB_tree(mCGALMesh, mCGALAABBTree);
    }

    ///@}
    ///@name Operatios
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

    ///@brief Returns boundary integration points of element.
    ///@param rElement
    ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
    bool ComputeBoundaryIps(Element& rElement, BoundaryIPVectorPtrType& rpBoundaryIps, const Parameters& rParam) const override;

    ///@brief Computes intersection mesh between rGeometry and rCube.
    ///@param rGeometry
    ///@param rCube
    ///@param rElement
    ///@param rParam
    ///@return bool Success?
    bool ComputeIntersectionMesh(const CGALMeshType& rGeometry, CGALMeshType& rCube,
                                 Element& rElement, const Parameters& rParam){

                                 }

private:

    CGALMeshType mCGALMesh;
    std::unique_ptr<CGALInsideTestType> mpCGALInsideTest;
    CGALAABBTreeType mCGALAABBTree;

}; // End Class BRepOperator

} // End namespace cgal

#endif // CGAL_BREP_OPERATOR_INCLUDE_H