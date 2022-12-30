// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MESH_UTILITIES_INCLUDE_H
#define MESH_UTILITIES_INCLUDE_H

//// Project includes
#include "containers/triangle_mesh.h"

namespace tibra {

///@name TIBRA Classes
///@{

///
/**
 * @class  MeshUtilities
 * @author Manuel Messmer
 * @brief  Provides operators for the TriangleMesh.
*/
class MeshUtilities {
public:
    ///@name Type Definitions
    ///@{
    typedef std::unique_ptr<TriangleMesh> TriangleMeshPtrType;

    ///@}
    ///@name Public Operations
    ///@{

    /// @brief Refines triangle mesh. Always conducts one refinement loop, such that area < 0.5*area_max.
    /// @param rTriangleMesh Triangle mesh to refine.
    /// @param MinNumberOfTriangles Minimum number of triangles in final mesh.
    /// @todo May needs improvement, more parameters etc.
    static void Refine(TriangleMesh& rTriangleMesh, IndexType MinNumberOfTriangles);

    /// @brief Appends rTriangleMesh by rNewMesh.
    /// @param rTriangleMesh Mesh to be appended.
    /// @param rNewMesh New mesh to be inserted in rTriangleMesh.
    static void Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh);

    /// @brief Append rTriangleMesh with some triangles in rNewMesh (given by Indices).
    /// @param rTriangleMesh Mesh to be appended.
    /// @param rNewMesh New mesh to be inserted in rTriangleMesh.
    /// @param rIndices Indices of triangles to be copied.
    static void Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh, const std::vector<IndexType>& rIndices);

    ///@brief Return meshed cuboid.
    ///@param rLowerPoint
    ///@param rUpperPoint
    ///@return std::unique_ptr<TriangleMesh>
    static TriangleMeshPtrType pGetCuboid(const PointType& rLowerPoint, const PointType& rUpperPoint);

    ///@}
}; // End class MeshUtilities
///@} // End TIBRA Classes
} // End namespace tibra

#endif // MESH_UTILITIES_INCLUDE_H