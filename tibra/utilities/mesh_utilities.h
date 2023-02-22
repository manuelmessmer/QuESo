// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MESH_UTILITIES_INCLUDE_H
#define MESH_UTILITIES_INCLUDE_H

//// Project includes
#include "containers/triangle_mesh.hpp"

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
    typedef Unique<TriangleMesh> TriangleMeshPtrType;

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
    ///@return Unique<TriangleMesh>
    static TriangleMeshPtrType pGetCuboid(const PointType& rLowerPoint, const PointType& rUpperPoint);

    ///@brief Returns enclosed volume by triangle mesh.
    ///@param rTriangleMesh
    ///@return double
    static double Volume(const TriangleMesh& rTriangleMesh);

    /// @brief Returns enclosed volume by triangle mesh. Uses divergence theorem to compute volume.
    /// Only aplpies divergence theorem in direction Dir: 0-x, 1-y, 2-z.
    /// @param rTriangleMesh
    /// @param Dir
    /// @return double
    static double Volume(const TriangleMesh& rTriangleMesh, IndexType Dir);

    ///@brief Returns enclosed volume by triangle mesh (OMP-version).
    ///@param rTriangleMesh
    ///@return double
    static double VolumeOMP(const TriangleMesh& rTriangleMesh);

    ///@brief Returns true if rTriangleMesh represents a closed volume.
    ///@param rTriangleMesh
    ///@param Tolerance Default: 1e-5.
    ///@return bool
    static bool IsClosed(const TriangleMesh& rTriangleMesh, double Tolerance = 1e-5);

    /// @brief Prototype so far.
    /// @param rTriangleMesh
    /// @return
    static double MaxAspectRatio(const TriangleMesh& rTriangleMesh);
    ///@}
}; // End class MeshUtilities
///@} // End TIBRA Classes
} // End namespace tibra

#endif // MESH_UTILITIES_INCLUDE_H