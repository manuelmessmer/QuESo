// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

//// STL includes
#include <fstream>
//// Project includes
#include "containers/element_container.h"
#include "containers/triangle_mesh.h"
#include "containers/boundary_integration_point.h"

namespace tibra {

class IO{

public:

  typedef std::size_t IndexType;
  typedef std::size_t SizeType;

  static bool WriteMeshToVTK(const TriangleMesh& rTriangleMesh,
                             const char* Filename,
                             const bool Binary);

  static bool WriteMeshToSTL(const TriangleMesh& rTriangleMesh,
                             const char* Filename,
                             const bool Binary);

  static bool ReadMeshFromSTL(TriangleMesh& rTriangleMesh,
                              const char* Filename);

  static bool WriteDisplacementToVTK(const std::vector<Vector3d>& rDisplacement,
                                     const char* Filename,
                                     const bool Binary);

  static bool WriteElementsToVTK(const ElementContainer& rElementContainer,
                                 const char* Filename,
                                 const bool Binary);

  static bool WritePointsToVTK(const ElementContainer& rElementContainer,
                               const char* Type,
                               const char* Filename,
                               const bool Binary);

  template<typename Type>
  static bool WritePointsToVTK(const std::vector<Type>& pPoints,
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

} // End namespace tibra

#endif // IO_UTILTIES_H