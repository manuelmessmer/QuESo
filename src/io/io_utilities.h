#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

#include "geometries/element_container.h"
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace IO{

template<typename PM>
void polygon_mesh_to_vtk(const PM& pmesh,//PolygonMesh
                                      const char* filename);

void WriteElementsToVTK(ElementContainer& rElementContainer,
                        const char* filename);

void WritePointsToVTK(ElementContainer& rElementContainer, const char* type,
                        const char* filename);

void WriteVTKPolyDataToVTK(const vtkSmartPointer<vtkPolyData> pPolyMesh,
                        const char* filename);
}

#endif // IO_UTILTIES_H