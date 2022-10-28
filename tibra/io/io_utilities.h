// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

// External includes
#include <fstream>      // std::ofstream

// Project includes
#include "geometries/element_container.h"
#include "geometries/triangle_mesh.h"
#include "geometries/boundary_integration_point.h"

class IO{

public:

  static bool WriteMeshToVTK(const TriangleMesh& rTriangleMesh,
                             const char* Filename,
                             const bool Binary);

  static bool WriteMeshToSTL(const TriangleMesh& rTriangleMesh,
                             const char* Filename,
                             const bool Binary);

  static bool ReadMeshFromSTL(TriangleMesh& rTriangleMesh,
                              const char* Filename);

  static bool WriteDisplacementToVTK(const std::vector<std::array<double,3>>& rDisplacement,
                                     const char* Filename,
                                     const bool Binary);

  static bool WriteElementsToVTK(const ElementContainer& rElementContainer,
                                 const char* Filename,
                                 const bool Binary);

  static bool WritePointsToVTK(const ElementContainer& rElementContainer,
                               const char* Type,
                               const char* Filename,
                               const bool Binary);

  /// TODO: templetize
  static bool WritePointsToVTK(const std::vector<BoundaryIntegrationPoint>& pPoints,
                               const char* Filename,
                               const bool Binary);

private:

  template<typename T>
  static void SwapEnd(T& var)
  {
    char* varArray = reinterpret_cast<char*>(&var);
    for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
      std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
  }

  template<typename T>
  static void WriteBinary(std::ofstream& stream, T& var){
    SwapEnd(var);
    stream.write(reinterpret_cast<char*>(&var), sizeof(T));
  }

};

#endif // IO_UTILTIES_H