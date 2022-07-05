
// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/number_utils.h>
#include <CGAL/IO/Complex_3_in_triangulation_3_to_vtk.h>
// Boost
//#include <CGAL/boost/graph/graph_traits.h>
#include <CGAL/boost/graph/properties.h>
// VTK includes
#include <vtkUnstructuredGrid.h>
#include <vtkMergeCells.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkPolyDataWriter.h>

// Project includes
#include "io/io_utilities.h"
#include "modeler/modeler.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Mesh; // Todo: Rename to Polyhedron
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;


template<typename PM>
void IO::polygon_mesh_to_vtk(const PM& pmesh,//PolygonMesh
                                      const char* filename)
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<PM, CGAL::vertex_point_t>::const_type VPMap;
  typedef typename boost::property_map_value<PM, CGAL::vertex_point_t>::type Point_3;
  VPMap vpmap = get(CGAL::vertex_point, pmesh);

  vtkPoints* const vtk_points = vtkPoints::New();
  vtkCellArray* const vtk_cells = vtkCellArray::New();

  vtk_points->Allocate(num_vertices(pmesh));
  vtk_cells->Allocate(num_faces(pmesh));

  std::map<vertex_descriptor, vtkIdType> Vids;
  vtkIdType inum = 0;

  for(vertex_descriptor v : vertices(pmesh))
  {
    const Point_3& p = get(vpmap, v);
    vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
                                CGAL::to_double(p.y()),
                                CGAL::to_double(p.z()));
    Vids[v] = inum++;
  }
  for(face_descriptor f : faces(pmesh))
  {
    vtkIdList* cell = vtkIdList::New();
    for(halfedge_descriptor h :
                  halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      cell->InsertNextId(Vids[target(h, pmesh)]);
    }
    vtk_cells->InsertNextCell(cell);
    cell->Delete();
  }

  vtkSmartPointer<vtkUnstructuredGrid> usg =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  usg->SetPoints(vtk_points);
  vtk_points->Delete();

  usg->SetCells(5,vtk_cells);
  vtk_cells->Delete();

  // Write the unstructured grid
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename);
  writer->SetInputData(usg);
  writer->Write();
}

template void IO::polygon_mesh_to_vtk<Mesh>(const Mesh& pmesh,//PolygonMesh
                                                const char* filename);
template void IO::polygon_mesh_to_vtk<SurfaceMesh>(const SurfaceMesh& pmesh,//PolygonMesh
                                                const char* filename);


void IO::WriteElementsToVTK(ElementContainer& rElementContainer, //PolygonMesh
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
    auto hexahedron = Modeler::GetVTKHexahedron(lower_point, upper_point);
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


void IO::WritePointsToVTK(ElementContainer& rElementContainer, const char* type,  //PolygonMesh
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
    auto points_it = (begin_points_it_ptr + i);
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

  WriteVTKPolyDataToVTK(vtk_poly_data, filename);
}

void IO::WriteVTKPolyDataToVTK(const vtkSmartPointer<vtkPolyData> pPolyMesh, //PolygonMesh
                               const char* filename)
{
  vtkSmartPointer<vtkPolyDataWriter> writer =
    vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(filename);
  writer->SetInputData(pPolyMesh);
  writer->Write();
 }


