
//VTK includes
#include <vtkCubeSource.h>
#include <vtkPolyDataNormals.h>

// Project includes
#include "modeler/cube_modeler.h"

vtkSmartPointer<vtkPolyData> CubeModeler::GetPolyData( std::array<double,3> lower_point, std::array<double,3> upper_point) {

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

vtkSmartPointer<vtkHexahedron> CubeModeler::GetHexahedron( std::array<double,3> lower_point, std::array<double,3> upper_point){

    auto hexahedron = vtkSmartPointer<vtkHexahedron>::New();
    for (int i = 0; i < hexahedron->GetNumberOfPoints(); ++i)
    {
        hexahedron->GetPointIds()->SetId(i, i);
    }
    hexahedron->GetPoints()->SetPoint(0,lower_point[0], lower_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(1,upper_point[0], lower_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(2,upper_point[0], upper_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(3,lower_point[0], upper_point[1], lower_point[2]);
    hexahedron->GetPoints()->SetPoint(4,lower_point[0], lower_point[1], upper_point[2]);
    hexahedron->GetPoints()->SetPoint(5,upper_point[0], lower_point[1], upper_point[2]);
    hexahedron->GetPoints()->SetPoint(6,upper_point[0], upper_point[1], upper_point[2]);
    hexahedron->GetPoints()->SetPoint(7,lower_point[0], upper_point[1], upper_point[2]);

    return hexahedron;
}

