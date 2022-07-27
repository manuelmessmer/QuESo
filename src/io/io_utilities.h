#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

#include "geometries/element_container.h"

namespace IO{

template<typename T>
void SwapEnd(T& var)
{
  char* varArray = reinterpret_cast<char*>(&var);
  for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
    std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

template<typename PM>
void polygon_mesh_to_vtk(const PM& pmesh,//PolygonMesh
                                      const char* filename,
                                      const bool binary);


void WriteElementsToVTK(ElementContainer& rElementContainer,
                        const char* filename, const bool binary);

void WritePointsToVTK(ElementContainer& rElementContainer, const char* type,
                        const char* filename,
                        const bool binary);

}

#endif // IO_UTILTIES_H