// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

/// CGAL includes
// Meshing/ Triangulation
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
// IO
#include <CGAL/IO/output_to_vtu.h>
//#include <fstream>
// Project includes
#include "cgal_wrapper/cgal_brep_operator.h"
#include "io/io_utilities.h"
#include "modeler/modeler.h"

#include <chrono>
#include <stdexcept>

namespace cgal {
// Namespaces
namespace PMP = CGAL::Polygon_mesh_processing;
using namespace CGAL::parameters;

typedef BRepOperatorBase::IntersectionStatus IntersectionStatus;
typedef BRepOperatorBase::BoundaryIPVectorType BoundaryIPVectorType;
typedef BRepOperatorBase::BoundaryIPVectorPtrType BoundaryIPVectorPtrType;

IntersectionStatus BRepOperator::GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance) const {
    // Test if center is inside or outside.
    const PointType center = { 0.5*(rUpperBound[0]+rLowerBound[0]), 0.5*(rUpperBound[1]+rLowerBound[1]), 0.5*(rUpperBound[2]+rLowerBound[2]) };
    auto status = (IsInside(center)) ? Inside : Outside;

    const CGALPointType point1(rLowerBound[0], rLowerBound[1], rLowerBound[2]);
    const CGALPointType point2(rUpperBound[0], rUpperBound[1], rUpperBound[2]);
    const CGAL::Iso_cuboid_3<CGALKernalType> tmp_cuboid( point1, point2, 0);
    if( mCGALAABBTree.do_intersect(tmp_cuboid) ){ // Perform inexact (conservative), but fast test based on AABBTree
        //if( PMP::do_intersect(rSurfaceMesh, rCubeMesh) ){ // Perform exact check
        status = Trimmed;
        //}
    }

    return status;
}

bool BRepOperator::IsInside(const PointType& rPoint) const {

    const CGALPointType tmp_point(rPoint[0], rPoint[1], rPoint[2]);
    // TODO: Add more test points, not only vertices of cube.
    CGAL::Bounded_side res = (*mpCGALInsideTest)(tmp_point);
    if (res == CGAL::ON_BOUNDED_SIDE) { return true; }
    if (res == CGAL::ON_BOUNDARY) { return false; }

    return false;
}

bool BRepOperator::ComputeBoundaryIps(Element& rElement, BoundaryIPVectorPtrType& rpBoundaryIps, const Parameters& rParam) const {


  //auto p_boundary_ips = std::make_unique<BoundaryIPVectorType>();

  // Warning: Do not use parallel tag, as global omp parallization is implemented.
  typedef CGAL::Sequential_tag Concurrency_tag;
  typedef boost::property_map<CGALMeshType, CGAL::edge_is_feature_t>::type CGALEIFMap;

  const auto& lower_bound = rElement.GetGlobalLowerPoint();
  const auto& upper_bound = rElement.GetGlobalUpperPoint();

  // Make sure cube is triangulated.
  //CGAL::Polygon_mesh_processing::triangulate_faces(rCube);

  // Make tmp copy, because PMP::corefine_and_compute_intersection does changes on geometry.
  CGALMeshType tmp_polyhedron = mCGALMesh;
  auto p_tmp_cube = Modeler::make_cube_3(lower_bound, upper_bound);

  // Compute intersection surface mesh
  CGALMeshType intersection_mesh{};
  bool valid_intersection = false;

  try {
    valid_intersection = PMP::corefine_and_compute_intersection(tmp_polyhedron, *p_tmp_cube, intersection_mesh);
  }
  catch(const std::exception& exc) {
    std::cerr << "CGAL :: CGAL BrepOperator :: ComputeIntersectionMesh :: No Valid Intersetion!" << std::endl;
    return 0;
  }
  if(!valid_intersection){
    std::cerr << "CGAL :: CGAL BrepOperator :: ComputeIntersectionMesh :: No Valid Intersetion!" << std::endl;
    return 0;
  }

  const double volume_cube = PMP::volume(*p_tmp_cube);
  const double volume_intersection_surface_mesh = PMP::volume(intersection_mesh);

  //Only consider intersections, which are larger than 0.1% with respect to the original element.
  if( volume_intersection_surface_mesh/volume_cube < 1e-3){
    // std::cout << "Warning :: Intersection neglected! Intersection Polyhedron To Cube Volume Ratio: "
    //   << volume_intersection_surface_mesh/volume_cube << std::endl;

    return 0;
  }

  // Construct ptr to SurfaceMesh
  std::unique_ptr<CGALMeshType> refinend_intersection_mesh = std::make_unique<CGALMeshType>();

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
    //     CGALPointType p(x,y,z);
    //     positions[*vi] = p;
    // }

    CGALEIFMap eif = CGAL::get(CGAL::edge_is_feature, *refinend_intersection_mesh);
    PMP::detect_sharp_edges(*refinend_intersection_mesh, 10, eif);

    unsigned int nb_iterations = 3;
    try {
      PMP::isotropic_remeshing(faces(*refinend_intersection_mesh), edge_length,
        *refinend_intersection_mesh,  PMP::parameters::number_of_iterations(nb_iterations).use_safety_constraints(false) // authorize all moves
          .edge_is_constrained_map(eif).number_of_relaxation_steps(1));
    }
    catch(const std::exception& exc) {
        if( !CGAL::is_closed(*refinend_intersection_mesh) ){
          std::cout << "Remeshing Exception: " << exc.what() << std::endl;
          std::cout << "Mesh is not closed. Knot spans will be neglected" << std::endl;
          return 0;
        }
    }
    edge_length = edge_length * 0.95*std::sqrt((double)refinend_intersection_mesh->number_of_faces() / (double) rParam.MinimumNumberOfTriangles()); // Todo: Make this better!!!
    iteration_count++;
  }

  // Element::PositionType positions = refinend_intersection_mesh->points();
  // for (auto vi = refinend_intersection_mesh->vertices_begin(); vi != refinend_intersection_mesh->vertices_end(); ++vi)
  // {
  //     double x = (positions[*vi].x() * std::abs(rParam.PointA()[0] - rParam.PointB()[0])) + rParam.PointA()[0];
  //     double y = (positions[*vi].y() * std::abs(rParam.PointA()[1] - rParam.PointB()[1])) + rParam.PointA()[1];
  //     double z = (positions[*vi].z() * std::abs(rParam.PointA()[2] - rParam.PointB()[2])) + rParam.PointA()[2];
  //     CGALPointType p(x,y,z);
  //     positions[*vi] = p;
  // }

  if (refinend_intersection_mesh->number_of_faces() < rParam.MinimumNumberOfTriangles() ){
    std::cout << "Warning:: Targeted number of triangles is not reached: "
      << refinend_intersection_mesh->number_of_faces() << " / " <<  rParam.MinimumNumberOfTriangles() << std::endl;
  }

  rElement.pSetSurfaceMesh( refinend_intersection_mesh ); // Ptr is lost after this function call.

  return 1;
}

} // End namespace cgal

