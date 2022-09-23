// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
// Boost
#include <CGAL/boost/graph/properties.h>

// External includes
//#include <fstream>      // std::ofstream

// Project includes
#include "io/io_utilities.h"
#include "modeler/modeler.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;
typedef std::size_t SizeType;

template<typename SM>
void IO::WriteMeshToVTK(const SM& rSurfaceMesh,
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
}

// Instantiation
template void IO::WriteMeshToVTK<SurfaceMesh>(const SurfaceMesh& rSurfaceMesh,//PolygonMesh
                                              const char* Filename,
                                              const bool Binary);


void IO::WriteDisplacementToVTK(const std::vector<std::array<double,3>>& rDisplacement,
                                const char* Filename,
                                const bool Binary){

  // const SizeType num_elements = rSurfaceMesh.number_of_faces();
  // const SizeType num_points = rSurfaceMesh.number_of_vertices();
  const SizeType num_points = rDisplacement.size();

  std::ofstream file;
  if(Binary)
    file.open(Filename, std::ios::app | std::ios::binary);
  else
    file.open(Filename);

  file << "POINT_DATA " << num_points << std::endl;
  file << "VECTORS Displacement double" << std::endl;
  for(int i = 0; i < num_points; ++i){
      if( Binary ){
        double rw1 = rDisplacement[i][0];
        WriteBinary(file, rw1);
        double rw2 = rDisplacement[i][1];
        WriteBinary(file, rw2);
        double rw3 = rDisplacement[i][2];
        WriteBinary(file, rw3);
      }
      else {
        //file << points_it->GetWeight() << std::endl;
      }
  }
  file << std::endl;

  file.close();
}

void IO::WriteElementsToVTK(const ElementContainer& rElementContainer, //PolygonMesh
                            const char* Filename,
                            const bool Binary){


  const SizeType num_elements = rElementContainer.size();

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
  file << "POINTS " << num_elements*8 << " double" << std::endl;

  const auto begin_el_itr = rElementContainer.begin();
  for( int i = 0; i < rElementContainer.size(); ++i){
    auto el_itr = *(begin_el_itr + i);
    auto lower_point = el_itr->GetGlobalLowerPoint();
    auto upper_point = el_itr->GetGlobalUpperPoint();

    if( Binary ){

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
    if( Binary ){
      int k = 8;
      WriteBinary(file, k);
      for( int j = 0; j < 8; ++j){
        k = 8*i+j;
        WriteBinary(file, k);
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
    if( Binary ){
        int k = 12;
        WriteBinary(file, k);
    }
    else {
      file << 12 << std::endl;
    }
  }
  file << std::endl;

  file.close();
}


void IO::WritePointsToVTK(const ElementContainer& rElementContainer,
                          const char* type,
                          const char* Filename,
                          const bool Binary){

  auto p_points = rElementContainer.pGetPoints(type);
  const auto begin_points_it_ptr = p_points->begin();
  const int num_points = p_points->size();
  const int num_elements = p_points->size();

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

  const Parameters& param = (*rElementContainer.begin())->GetParameters();
  for(int i = 0; i < num_points; ++i){
    auto points_it = (begin_points_it_ptr + i);
    auto point_global = MappingUtilities::FromLocalToGlobalSpace(*points_it, param.PointA(), param.PointB() );

    if( Binary ){
      WriteBinary(file, point_global[0]);
      WriteBinary(file, point_global[1]);
      WriteBinary(file, point_global[2]);
    }
    else {
      file << point_global[0] << ' ' << point_global[1] << ' ' << point_global[2] << std::endl;
    }
  }
  file << std::endl;

  //Write Cells
  file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
      int k = 1;
      WriteBinary(file, k);
      k = i;
      WriteBinary(file, k);
    }
    else {
      file << 1 << ' ' << i << std::endl;
    }
  }
  file << std::endl;

  file << "CELL_TYPES " << num_elements << std::endl;
  for( int i = 0; i < num_elements; ++i){
    if( Binary ){
        int k = 1;
        WriteBinary(file, k);
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

      if( Binary ){
        double rw = points_it->GetWeight();
        WriteBinary(file, rw);
      }
      else {
        file << points_it->GetWeight() << std::endl;
      }
  }
  file << std::endl;

  file.close();
}



