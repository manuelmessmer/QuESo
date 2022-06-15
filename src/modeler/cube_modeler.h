
#ifndef CUBE_MODELER_INCLUDE_H
#define CUBE_MODELER_INCLUDE_H

//VTK includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkHexahedron.h>


namespace CubeModeler {

typedef std::array<double,3> PointType;

vtkSmartPointer<vtkPolyData> GetPolyData( std::array<double,3> lower_point, std::array<double,3> upper_point);

vtkSmartPointer<vtkHexahedron> GetHexahedron( std::array<double,3> lower_point, std::array<double,3> upper_point);

}// End Namespace TestHelper

#endif // CUBE_MODELER_INCLUDE_H