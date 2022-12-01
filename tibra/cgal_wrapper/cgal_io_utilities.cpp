// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
// Boost
#include <CGAL/boost/graph/properties.h>
//// Project includes
#include "cgal_wrapper/cgal_io_utilities.h"

namespace tibra {
namespace cgal {

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernalType;
typedef CGAL::Surface_mesh<CGALKernalType::Point_3> CGALMeshType;
typedef std::size_t IndexType;
typedef std::size_t SizeType;

template<typename SM>
bool IO::WriteMeshToVTK(const SM& rSurfaceMesh,
                        const char* Filename,
                        const bool Binary)
{
  typedef typename boost::graph_traits<SM>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<SM>::face_descriptor     face_descriptor;
  typedef typename boost::graph_traits<SM>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<SM, CGAL::vertex_point_t>::const_type VPMap;
  typedef typename boost::property_map_value<SM, CGAL::vertex_point_t>::type Point_3;

  VPMap vpmap = CGAL::get(CGAL::vertex_point, rSurfaceMesh);

  const SizeType num_elements = rSurfaceMesh.number_of_faces();
  const SizeType num_points = rSurfaceMesh.number_of_vertices();



  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::out | std::ios::binary);
  else
    file.open(Filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(Binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;

  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_points << " double" << std::endl;

  std::map<vertex_descriptor, IndexType> vids;
  IndexType inum = 0;
  auto raw_points = rSurfaceMesh.points().data()->cartesian_begin();
  for(vertex_descriptor v : vertices(rSurfaceMesh))
  {
    const Point_3& p = get(vpmap, v);
    if( Binary ){
      double rx = CGAL::to_double(p.x());
      double ry = CGAL::to_double(p.y());
      double rz = CGAL::to_double(p.z());

      WriteBinary(file, rx);
      WriteBinary(file, ry);
      WriteBinary(file, rz);
    }
    else {
      file << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << ' ' << CGAL::to_double(p.z()) << std::endl;
    }
    vids[v] = inum++;
  }
  file << std::endl;

  // Write Cells
  file << "Cells " << num_elements << " " << num_elements*4 << std::endl;

  for(face_descriptor f : faces(rSurfaceMesh))
  {
    if( Binary ){
      int k = 3;
      WriteBinary(file, k);
      for(halfedge_descriptor h :
                    halfedges_around_face(halfedge(f, rSurfaceMesh), rSurfaceMesh))
      {
        k = vids[target(h, rSurfaceMesh)];
        WriteBinary(file, k);
      }
    }
    else {
      file << 3;
      for(halfedge_descriptor h :
                    halfedges_around_face(halfedge(f, rSurfaceMesh), rSurfaceMesh))
      {
          file << ' ' << vids[target(h, rSurfaceMesh)];
      }
      file << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
        int k = 5;
        WriteBinary(file, k);
    }
    else {
      file << 5 << std::endl;
    }
  }
  file << std::endl;
  file.close();

  return true;
}

// Instantiation
template bool IO::WriteMeshToVTK<CGALMeshType>(const CGALMeshType& rSurfaceMesh,//PolygonMesh
                                              const char* Filename,
                                              const bool Binary);

} // End namespace cgal
} // End namespace tibra