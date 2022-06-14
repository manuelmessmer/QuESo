
#ifndef CUBE_MODELER_INCLUDE_H
#define CUBE_MODELER_INCLUDE_H

//VTK includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataNormals.h>

typedef std::array<double,3> PointType;

namespace CubeModeler {

vtkSmartPointer<vtkPolyData> make_cube_3( std::array<double,3> lower_point, std::array<double,3> upper_point) {

    vtkSmartPointer<vtkCubeSource> new_cube = vtkSmartPointer<vtkCubeSource>::New();
    std::array<double,6> bounds = {lower_point[0], upper_point[0], lower_point[1], upper_point[1], lower_point[2], upper_point[2]};
    new_cube->SetBounds( bounds.data() );
    new_cube->Update();

    vtkSmartPointer<vtkPolyDataNormals> normal_generator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_generator->SetInputData(new_cube->GetOutput());
    normal_generator->ComputePointNormalsOff();
    normal_generator->ComputeCellNormalsOn();
    normal_generator->Update();

    return normal_generator->GetOutput();
}

}// End Namespace TestHelper

#endif // CUBE_MODELER_INCLUDE_H