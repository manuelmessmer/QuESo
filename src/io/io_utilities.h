// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// Project includes
#include "geometries/element_container.h"

namespace IO_Utilities {

void WriteVTK(const vtkSmartPointer<vtkPolyData> pPolyMesh, //PolygonMesh
                const char* filename);

void WriteElementsToVTK(ElementContainer& rElementContainer, //PolygonMesh
                        const char* filename);

} // end namespace

#endif // IO_UTILTIES