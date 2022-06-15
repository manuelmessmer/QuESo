// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// VTK includes
#include <vtkPolyDataWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkMergeCells.h>
#include <vtkIdTypeArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataObject.h>
#include <vtkDataWriter.h>
#include <vtkVertex.h>

// Project includes
#include "io/io_utilities.h"
#include "utilities/mapping_utilities.h"
#include "modeler/cube_modeler.h"

void IO_Utilities::WritePolyDataToVTK(const vtkSmartPointer<vtkPolyData> pPolyMesh, //PolygonMesh
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

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // Todo: makehexahedron
    auto merge = vtkSmartPointer<vtkMergeCells>::New();

    merge->SetTotalNumberOfPoints(8*rElementContainer.size());
    merge->SetTotalNumberOfCells(rElementContainer.size());
    merge->SetTotalNumberOfDataSets(rElementContainer.size());
    merge->SetPointMergeTolerance(1e-5);
    merge->SetUseGlobalCellIds(1);
    merge->SetUseGlobalIds(0);


    merge->SetUnstructuredGrid(grid);
    const auto begin_el_itr = rElementContainer.begin();
    for( int i = 0; i < rElementContainer.size(); ++i){
      auto el_itr = *(begin_el_itr + i);
      auto lower_point = el_itr->GetGlobalLowerPoint();
      auto upper_point = el_itr->GetGlobalUpperPoint();

      auto tmp_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
      auto hexahedron = CubeModeler::GetHexahedron(lower_point, upper_point);
      tmp_grid->SetPoints(hexahedron->GetPoints());
      tmp_grid->InsertNextCell(hexahedron->GetCellType(), hexahedron->GetPointIds() );
      auto ids1 = vtkSmartPointer<vtkIdTypeArray>::New();
      ids1->SetName("GlobalCellIds");
      ids1->SetNumberOfValues(1);
      tmp_grid->GetCellData()->SetGlobalIds(ids1);
      tmp_grid->GetCellData()->GetGlobalIds()->SetTuple1(0,i);

      merge->MergeDataSet(tmp_grid);

    }
    merge->Finish();

    auto grid_writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    grid_writer->SetFileName(filename);
    grid_writer->SetInputData(grid);
    grid_writer->Write();
  }

void IO_Utilities::WritePointsToVTK(ElementContainer& rElementContainer, const char* type,  //PolygonMesh
                                    const char* filename){
    auto p_points = rElementContainer.pGetPoints(type);
    auto vtk_points = vtkSmartPointer<vtkPoints>::New();
    auto vtk_verts = vtkSmartPointer<vtkCellArray>::New();


    auto weights = vtkSmartPointer<vtkDoubleArray>::New();
    auto polys = vtkSmartPointer<vtkCellArray>::New();
    polys->InitTraversal();

    const auto begin_points_it_ptr = p_points->begin();
    const Parameters& param = (*rElementContainer.begin())->GetParameters();
    for(int i = 0; i < p_points->size(); ++i){
      auto points_it = *(begin_points_it_ptr + i);
      auto point_global = MappingUtilities::FromLocalToGlobalSpace(*points_it, param.PointA(), param.PointB() );
      auto id =  vtk_points->InsertNextPoint(point_global.data());
      weights->InsertNextTuple1(points_it->GetWeight());

      vtkSmartPointer<vtkVertex> tmp_vertex = vtkSmartPointer<vtkVertex>::New();
      tmp_vertex->GetPointIds()->SetId(0,i);
      polys->InsertNextCell(tmp_vertex);
    }
    auto vtk_poly_data = vtkSmartPointer<vtkPolyData>::New();

    vtk_poly_data->SetPoints(vtk_points);
    vtk_poly_data->SetPolys(polys);
    auto point_data = vtk_poly_data->GetPointData();
    weights->SetName("Weights");
    point_data->SetScalars(weights);

    WritePolyDataToVTK(vtk_poly_data, filename);
}
