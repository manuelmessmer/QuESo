
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
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>

#include "io/vtk_utilities.h"

namespace CGAL{

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Mesh; // Todo: Rename to Polyhedron
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;

template<typename PM>
void polygon_mesh_to_vtkUnstructured_(const PM& pmesh,//PolygonMesh
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

template void polygon_mesh_to_vtkUnstructured_<Mesh>(const Mesh& pmesh,//PolygonMesh
                                                const char* filename);
template void polygon_mesh_to_vtkUnstructured_<SurfaceMesh>(const SurfaceMesh& pmesh,//PolygonMesh
                                                const char* filename);
}
