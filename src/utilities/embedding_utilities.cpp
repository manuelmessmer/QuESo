// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

/// CGAL includes
// Domain
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h> //Todo: Can this be removed?
// Meshing/ Triangulation
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
// IO
#include <CGAL/IO/output_to_vtu.h>
//#include <fstream>
// Project includes
#include "utilities/embedding_utilities.h"
#include "io/io_utilities.h"

#include <chrono>
#include <stdexcept>

// Typedefs
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
// Warning: Do not use parallel tag, as global omp parallization is implemented.
typedef CGAL::Sequential_tag Concurrency_tag;
typedef K::Vector_3 Vector;

typedef boost::property_map<CGAL::Surface_mesh<Point_3>, CGAL::edge_is_feature_t>::type EIFMap;

// Namespaces
namespace PMP = CGAL::Polygon_mesh_processing;
using namespace CGAL::parameters;

bool EmbeddingUtilities::ComputeIntersectionMesh(const SurfaceMeshType& rGeometry, SurfaceMeshType& rCube, Element& rElement, const Parameters& rParam){

  // Make sure cube is triangulated.
  CGAL::Polygon_mesh_processing::triangulate_faces(rCube);

  // Make tmp copy, because PMP::corefine_and_compute_intersection does changes on geometry.
  SurfaceMeshType tmp_polyhedron = rGeometry;
  SurfaceMeshType tmp_cube = rCube;

  // Compute intersection surface mesh
  SurfaceMeshType intersection_mesh;
  bool valid_intersection = false;
  try {
    valid_intersection = PMP::corefine_and_compute_intersection(tmp_polyhedron, tmp_cube, intersection_mesh);
  }
  catch(const std::exception& exc) {
    std::cerr << "No Valid Intersetion!" << std::endl;
    return 0;
  }
  if(!valid_intersection){
    std::cerr << "No Valid Intersetion!" << std::endl;
    return 0;
  }
  const double volume_cube = CGAL::Polygon_mesh_processing::volume(tmp_cube);
  const double volume_intersection_surface_mesh = CGAL::Polygon_mesh_processing::volume(intersection_mesh);

  //Only consider intersections, which are larger than 0.5% with respect to the original element.
  // if( volume_intersection_surface_mesh/volume_cube < 0.001){
  //   // std::cout << "Warning :: Intersection neglected! Intersection Polyhedron To Cube Volume Ratio: "
  //   //   << volume_intersection_surface_mesh/volume_cube << std::endl;

  //   return 0;
  // }

  // Construct ptr to SurfaceMesh
  std::unique_ptr<SurfaceMeshType> refinend_intersection_mesh = std::make_unique<SurfaceMeshType>();

  // Remesh intersected domain until minimum number of boundary triangles is reached.
  double edge_length = rParam.InitialTriangleEdgeLength();
  int iteration_count = 0;
  while (refinend_intersection_mesh->number_of_faces() < rParam.MinimumNumberOfTriangles() && iteration_count < 10){
    refinend_intersection_mesh->clear();
    CGAL::copy_face_graph(intersection_mesh,*refinend_intersection_mesh);

    // Element::PositionType positions = refinend_intersection_mesh->points();
    // for (auto vi = refinend_intersection_mesh->vertices_begin(); vi != refinend_intersection_mesh->vertices_end(); ++vi)
    // {
    //     double x = (positions[*vi].x() - rParam.PointA()[0])/std::abs(rParam.PointA()[0] - rParam.PointB()[0]);
    //     double y = (positions[*vi].y() - rParam.PointA()[1])/std::abs(rParam.PointA()[1] - rParam.PointB()[1]);
    //     double z = (positions[*vi].z() - rParam.PointA()[2])/std::abs(rParam.PointA()[2] - rParam.PointB()[2]);
    //     Point_3 p(x,y,z);
    //     positions[*vi] = p;
    // }

    EIFMap eif = CGAL::get(CGAL::edge_is_feature, *refinend_intersection_mesh);
    PMP::detect_sharp_edges(*refinend_intersection_mesh, 10, eif);

    unsigned int nb_iterations = 3;
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(*refinend_intersection_mesh), edge_length,
                  *refinend_intersection_mesh,  PMP::parameters::number_of_iterations(nb_iterations).use_safety_constraints(false) // authorize all moves
                                        .edge_is_constrained_map(eif));
    edge_length = edge_length * 0.95*std::sqrt((double)refinend_intersection_mesh->number_of_faces() / (double) rParam.MinimumNumberOfTriangles()); // Todo: Make this better!!!
    iteration_count++;
  }

  // Element::PositionType positions = refinend_intersection_mesh->points();
  // for (auto vi = refinend_intersection_mesh->vertices_begin(); vi != refinend_intersection_mesh->vertices_end(); ++vi)
  // {
  //     double x = (positions[*vi].x() * std::abs(rParam.PointA()[0] - rParam.PointB()[0])) + rParam.PointA()[0];
  //     double y = (positions[*vi].y() * std::abs(rParam.PointA()[1] - rParam.PointB()[1])) + rParam.PointA()[1];
  //     double z = (positions[*vi].z() * std::abs(rParam.PointA()[2] - rParam.PointB()[2])) + rParam.PointA()[2];
  //     Point_3 p(x,y,z);
  //     positions[*vi] = p;
  // }

  if (refinend_intersection_mesh->number_of_faces() < rParam.MinimumNumberOfTriangles() ){
    std::cout << "Warning:: Targeted number of triangles is not reached: "
      << refinend_intersection_mesh->number_of_faces() << " / " <<  rParam.MinimumNumberOfTriangles() << std::endl;
  }

  rElement.pSetSurfaceMesh(refinend_intersection_mesh); // Ptr is lost after this function call.

  return 1;
}

