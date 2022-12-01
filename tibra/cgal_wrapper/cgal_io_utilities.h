// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CGAL_IO_UTILTIES_H
#define CGAL_IO_UTILTIES_H

//// STL includes
#include <fstream>

namespace tibra {
namespace cgal {

class IO{

public:
  template<typename SM>
  static bool WriteMeshToVTK(const SM& rSurfaceMesh,
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

} // End namespace cgal
} // End namespace tibra

#endif // CGAL_IO_UTILTIES_H