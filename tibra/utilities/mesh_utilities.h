// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MESH_UTILITIES_INCLUDE_H
#define MESH_UTILITIES_INCLUDE_H

//// Project includes
#include "containers/triangle_mesh.h"

namespace tibra {

class MeshUtilities {

public:
    static void Refine(TriangleMesh& rTriangleMesh, IndexType MinNumberOfTriangles);

    static void Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh);

    static void Append(TriangleMesh& rTriangleMesh, const TriangleMesh& rNewMesh, const std::vector<IndexType>& rIndices);
}; // End class MeshUtilities

} // End namespace tibra

#endif // MESH_UTILITIES_INCLUDE_H