
// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/number_utils.h>

// Boost
//#include <CGAL/boost/graph/graph_traits.h>
#include <CGAL/boost/graph/properties.h>

// Project includes
#include "io/io_utilities.h"
#include "modeler/modeler.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Mesh; // Todo: Rename to Polyhedron
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;
typedef std::size_t SizeType;

template<typename PM>
void IO::polygon_mesh_to_vtk(const PM& pmesh,//PolygonMesh
                                      const char* filename, const bool binary)
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<PM, CGAL::vertex_point_t>::const_type VPMap;
  typedef typename boost::property_map_value<PM, CGAL::vertex_point_t>::type Point_3;

  VPMap vpmap = CGAL::get(CGAL::vertex_point, pmesh);

  const SizeType num_elements = pmesh.number_of_faces();
  const SizeType num_points = pmesh.number_of_vertices();

  std::ofstream file;
  if(binary)
    file.open(filename, std::ios::out | std::ios::binary);
  else
    file.open(filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;

  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_points << " double" << std::endl;


  std::map<vertex_descriptor, IndexType> Vids;
  IndexType inum = 0;
  for(vertex_descriptor v : vertices(pmesh))
  {
    const Point_3& p = get(vpmap, v);
    if( binary ){

      double rx = CGAL::to_double(p.x());
      double ry = CGAL::to_double(p.y());
      double rz = CGAL::to_double(p.z());
      SwapEnd(rx);
      SwapEnd(ry);
      SwapEnd(rz);

      file.write(reinterpret_cast<char*>(&rx), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz), sizeof(double));

    }
    else {
      file << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << ' ' << CGAL::to_double(p.z()) << std::endl;

    }
    Vids[v] = inum++;
  }
  file << std::endl;

  // Write Cells
  file << "Cells " << num_elements << " " << num_elements*4 << std::endl;

  for(face_descriptor f : faces(pmesh))
  {
    if( binary ){
      int k = 3;
      SwapEnd(k);
      file.write(reinterpret_cast<char*>(&k), sizeof(int));
      for(halfedge_descriptor h :
                    halfedges_around_face(halfedge(f, pmesh), pmesh))
      {
        k = Vids[target(h, pmesh)];
        SwapEnd(k);
        file.write(reinterpret_cast<char*>(&k), sizeof(int));

      }
    }
    else {
      file << 3;
      for(halfedge_descriptor h :
                    halfedges_around_face(halfedge(f, pmesh), pmesh))
      {
          file << ' ' << Vids[target(h, pmesh)];
      }
      file << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( binary ){
        int k = 5;
        SwapEnd(k);
        file.write(reinterpret_cast<char*>(&k), sizeof(int));
    }
    else {
      file << 5 << std::endl;
    }
  }
  file << std::endl;

  // file << "CELL_DATA 4" << std::endl;
  // file << "GLOBAL_IDS GlobalCellIds vtkIdType" << std::endl;
  // for( int i = 0; i < num_elements; ++i){
  //     if( binary ){
  //       int k = i;
  //       SwapEnd(k);
  //       file.write(reinterpret_cast<char*>(&k), sizeof(int));
  //     }
  //     else {
  //       file << i << ' ';
  //       if( i%4 == 0 && i > 0){
  //         file << std::endl;
  //       }
  //     }
  // }
  // file << std::endl;

  file.close();

}

// template void IO::polygon_mesh_to_vtk<Mesh>(const Mesh& pmesh,//PolygonMesh
//                                                 const char* filename,
//                                                 const bool binary);
template void IO::polygon_mesh_to_vtk<SurfaceMesh>(const SurfaceMesh& pmesh,//PolygonMesh
                                                const char* filename,
                                                const bool binary);

void IO::WriteElementsToVTK(ElementContainer& rElementContainer, //PolygonMesh
                        const char* filename, const bool binary){


  const SizeType num_elements = rElementContainer.size();

  std::ofstream file;
  if(binary)
    file.open(filename, std::ios::out | std::ios::binary);
  else
    file.open(filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;

  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_elements*8 << " double" << std::endl;

  const auto begin_el_itr = rElementContainer.begin();
  for( int i = 0; i < rElementContainer.size(); ++i){
    auto el_itr = *(begin_el_itr + i);
    auto lower_point = el_itr->GetGlobalLowerPoint();
    auto upper_point = el_itr->GetGlobalUpperPoint();

    if( binary ){

      double rx0 = lower_point[0];
      SwapEnd(rx0);
      double rx1 = upper_point[0];
      SwapEnd(rx1);
      double ry0 = lower_point[1];
      SwapEnd(ry0);
      double ry1 = upper_point[1];
      SwapEnd(ry1);
      double rz0 = lower_point[2];
      SwapEnd(rz0);
      double rz1 = upper_point[2];
      SwapEnd(rz1);

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

      file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz1), sizeof(double));
    }
    else {
      file << lower_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
      file << upper_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
      file << upper_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
      file << lower_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
      file << lower_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
      file << upper_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
      file << upper_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
      file << lower_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
    }
  }
  file << std::endl;
  // Write Cells
  file << "Cells " << num_elements << " " << num_elements*9 << std::endl;
  for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
    if( binary ){
      int k = 8;
      SwapEnd(k);
      file.write(reinterpret_cast<char*>(&k), sizeof(int));
      for( int j = 0; j < 8; ++j){
        k = 8*i+j;
        SwapEnd(k);
        file.write(reinterpret_cast<char*>(&k), sizeof(int));
      }
    }
    else {
      file << 8 << ' ' << 8*i     << ' ' << 8*i + 1 << ' ' << 8*i + 2 << ' ' << 8*i + 3
                << ' ' << 8*i + 4 << ' ' << 8*i + 5 << ' ' << 8*i + 6 << ' ' << 8*i + 7 << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << rElementContainer.size() << std::endl;
  for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
    if( binary ){
        int k = 12;
        SwapEnd(k);
        file.write(reinterpret_cast<char*>(&k), sizeof(int));
    }
    else {
      file << 12 << std::endl;
    }
  }
  file << std::endl;

  // file << "CELL_DATA 4" << std::endl;
  // file << "GLOBAL_IDS GlobalCellIds vtkIdType" << std::endl;
  // for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
  //     if( binary ){
  //       int k = i;
  //       SwapEnd(k);
  //       file.write(reinterpret_cast<char*>(&k), sizeof(int));
  //     }
  //     else {
  //       file << i << ' ';
  //       if( i%4 == 0 && i > 0){
  //         file << std::endl;
  //       }
  //     }
  // }
  // file << std::endl;

  file.close();
}


void IO::WritePointsToVTK(ElementContainer& rElementContainer, const char* type,
                                    const char* filename,
                                    const bool binary){

  auto p_points = rElementContainer.pGetPoints(type);
  const auto begin_points_it_ptr = p_points->begin();
  const int num_points = p_points->size();
  const int num_elements = p_points->size();

  std::ofstream file;
  if(binary)
    file.open(filename, std::ios::out | std::ios::binary);
  else
    file.open(filename);

  file << "# vtk DataFile Version 4.1" << std::endl;
  file << "vtk output" << std::endl;
  if(binary)
    file << "BINARY"<< std::endl;
  else
    file << "ASCII"<< std::endl;


  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  file << "POINTS " << num_points << " double" << std::endl;

  const Parameters& param = (*rElementContainer.begin())->GetParameters();
  for(int i = 0; i < num_points; ++i){
    auto points_it = (begin_points_it_ptr + i);
    auto point_global = MappingUtilities::FromLocalToGlobalSpace(*points_it, param.PointA(), param.PointB() );

    if( binary ){

      double rx = point_global[0];
      double ry = point_global[1];
      double rz = point_global[2];
      SwapEnd(rx);
      SwapEnd(ry);
      SwapEnd(rz);

      file.write(reinterpret_cast<char*>(&rx), sizeof(double));
      file.write(reinterpret_cast<char*>(&ry), sizeof(double));
      file.write(reinterpret_cast<char*>(&rz), sizeof(double));

    }
    else {
      file << point_global[0] << ' ' << point_global[1] << ' ' << point_global[2] << std::endl;
    }
  }
  file << std::endl;

  //Write Cells
  file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( binary ){
      int k = 1;
      SwapEnd(k);
      file.write(reinterpret_cast<char*>(&k), sizeof(int));

      k = i;
      SwapEnd(k);
      file.write(reinterpret_cast<char*>(&k), sizeof(int));

    }
    else {
      file << 1 << ' ' << i << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( binary ){
        int k = 1;
        SwapEnd(k);
        file.write(reinterpret_cast<char*>(&k), sizeof(int));
    }
    else {
      file << 1 << std::endl;
    }
  }
  file << std::endl;

  file << "POINT_DATA " << num_points << std::endl;
  file << "SCALARS Weights double 1" << std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  for(int i = 0; i < num_points; ++i){
      auto points_it = (begin_points_it_ptr + i);

      if( binary ){
        double rw = points_it->GetWeight();
        SwapEnd(rw);
        file.write(reinterpret_cast<char*>(&rw), sizeof(double));
      }
      else {
        file << points_it->GetWeight() << std::endl;
      }
  }
  file << std::endl;

  file.close();
}



