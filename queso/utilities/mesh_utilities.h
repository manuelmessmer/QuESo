//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef MESH_UTILITIES_INCLUDE_H
#define MESH_UTILITIES_INCLUDE_H

//// STL includes
#include <utility>

//// Project includes
#include "queso/containers/triangle_mesh_view.hpp"

namespace queso {

///@name QuESo Classes
///@{

///
/**
 * @class  MeshUtilities
 * @author Manuel Messmer
 * @brief  Provides operators for the TriangleMesh.
*/
namespace MeshUtilities {

using TriangleMeshPtrType = Unique<TriangleMesh>;

/// @brief Refines a mesh until a minimum triangle count is reached.
/// @param rTriangleMesh Mesh that is refined in-place.
/// @param MinNumberOfTriangles Requested minimum number of triangles.
/// @details Refinement keeps the original surface geometry while increasing triangle resolution.
void Refine(TriangleMesh& rTriangleMesh, IndexType MinNumberOfTriangles);

/// @brief Appends one mesh to another.
/// @param rTriangleMesh Destination mesh (modified in-place).
/// @param rNewMesh Source mesh that is appended to the destination.
/// @details Vertex and triangle connectivity are remapped consistently during append.
void Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh);

/// @brief Creates a cuboid surface mesh from axis-aligned bounds.
/// @param rLowerPoint Lower corner of the cuboid in global coordinates.
/// @param rUpperPoint Upper corner of the cuboid in global coordinates.
/// @return Unique pointer to the generated cuboid triangle mesh.
TriangleMeshPtrType pGetCuboid(const PointType& rLowerPoint, const PointType& rUpperPoint);

/// @brief Computes total surface area of a triangle mesh view.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Total mesh area.
double Area(const TriangleMeshView& rTriangleMeshView);

/// @brief Computes total surface area in parallel.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Total mesh area.
/// @note Uses OpenMP parallelization when available.
double AreaOMP(const TriangleMeshView& rTriangleMeshView);

/// @brief Computes enclosed volume of a closed triangle mesh.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Signed volume enclosed by the mesh.
double Volume(const TriangleMeshView& rTriangleMeshView);

/// @brief Computes directional contribution to enclosed volume.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @param Dir Direction index used for directional volume evaluation.
/// @return Signed directional volume contribution.
double Volume(const TriangleMeshView& rTriangleMeshView, IndexType Dir);

/// @brief Computes enclosed volume in parallel.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Signed volume enclosed by the mesh.
/// @note Uses OpenMP parallelization when available.
double VolumeOMP(const TriangleMeshView& rTriangleMeshView);

/// @brief Estimates overall mesh quality.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Scalar quality estimate based on triangle shape metrics.
double EstimateQuality(const TriangleMeshView& rTriangleMeshView);

/// @brief Returns the maximum triangle aspect ratio in the mesh.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Maximum aspect ratio over all triangles.
double MaxAspectRatio(const TriangleMeshView& rTriangleMeshView);

/// @brief Returns the average triangle aspect ratio in the mesh.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Average aspect ratio over all triangles.
double AverageAspectRatio(const TriangleMeshView& rTriangleMeshView);

/// @brief Computes axis-aligned bounding box of a triangle mesh.
/// @param rTriangleMeshView Non-owning view of the mesh.
/// @return Pair of points `{lower_bound, upper_bound}` in global coordinates.
std::pair<PointType, PointType> BoundingBox(const TriangleMeshView& rTriangleMeshView);

} // namespace MeshUtilities
///@} // End QuESo Classes
} // End namespace queso

#endif // MESH_UTILITIES_INCLUDE_H

