
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MESH_UTILITIES_INCLUDE_H
#define MESH_UTILITIES_INCLUDE_H

// Vtk inlude
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

// Project include
#include "utilities/parameters.h"
#include <chrono>

namespace MeshUtilities {
    bool DoIntersect(vtkSmartPointer<vtkPolyData> pInput1, vtkSmartPointer<vtkPolyData> pInput2);

    vtkSmartPointer<vtkPolyData> GetIntersection(vtkSmartPointer<vtkPolyData> pInput1, vtkSmartPointer<vtkPolyData> pInput2, const Parameters& rParam);

    vtkSmartPointer<vtkPolyData> IsotropicRemeshing(vtkSmartPointer<vtkPolyData> pInput, const Parameters& rParam, int id);
}

#endif //MESH_UTILITIES_INCLUDE_H