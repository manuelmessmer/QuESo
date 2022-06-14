// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// VTK includes
#include <vtkPolyDataWriter.h>

// Project includes
#include "io/io_utilities.h"

void IO_Utilities::WriteVTK(const vtkSmartPointer<vtkPolyData> pPolyMesh, //PolygonMesh
                const char* filename)
  {

    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(pPolyMesh);
    writer->Write();
  }

