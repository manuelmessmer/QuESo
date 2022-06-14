// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// VTK includes
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>

// Project includes
#include "io/io_utilities.h"
#include "utilities/mapping_utilities.h"
#include "modeler/cube_modeler.h"

void IO_Utilities::WriteVTK(const vtkSmartPointer<vtkPolyData> pPolyMesh, //PolygonMesh
                const char* filename)
  {

    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(pPolyMesh);
    writer->Write();
  }

void IO_Utilities::WriteElementsToVTK(ElementContainer& rElementContainer, //PolygonMesh
                        const char* filename){

    vtkSmartPointer<vtkPolyData> tmp_poly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkAppendPolyData> append_filter = vtkSmartPointer<vtkAppendPolyData>::New();
    // Todo: makehexahedron
    const auto begin_el_itr = rElementContainer.begin();
    for( int i = 0; i < rElementContainer.size(); ++i){
      auto el_itr = *(begin_el_itr + i);
      auto lower_point = el_itr->GetGlobalLowerPoint();
      auto upper_point = el_itr->GetGlobalUpperPoint();
      auto cube_tmp = CubeModeler::make_cube_3(lower_point, upper_point);


      append_filter->AddInputData(cube_tmp);
      append_filter->AddInputData(tmp_poly);
      append_filter->Update();
      tmp_poly = append_filter->GetOutput();
    }

    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(tmp_poly);
    writer->Write();
  }

