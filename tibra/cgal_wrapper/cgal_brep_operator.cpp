// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// CGAL includes
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
//// STL includes
#include <chrono>
#include <stdexcept>
//// Project includes
#include "cgal_wrapper/cgal_brep_operator.h"
#include "cgal_wrapper/cgal_cuboid_modeler.h"
#include "cgal_wrapper/cgal_trimmed_domain.h"
#include "cgal_wrapper/cgal_utilities.h"

namespace tibra {
namespace cgal {

namespace PMP = CGAL::Polygon_mesh_processing;
using namespace CGAL::parameters;

typedef BRepOperatorBase::IntersectionStatus IntersectionStatus;
typedef BRepOperatorBase::TrimmedDomainBasePtrType TrimmedDomainBasePtrType;

IntersectionStatus CGALBRepOperator::GetIntersectionState(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance) const {
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

bool CGALBRepOperator::IsInside(const PointType& rPoint) const {

    // Construct cgal point
    const CGALPointType tmp_point(rPoint[0], rPoint[1], rPoint[2]);

    // Test point
    CGAL::Bounded_side res = (*mpCGALInsideTest)(tmp_point);
    if (res == CGAL::ON_BOUNDED_SIDE) { return true; }
    if (res == CGAL::ON_BOUNDARY) { return false; }

    return false;
}

TrimmedDomainBasePtrType CGALBRepOperator::GetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound, const Parameters& rParam) const {

    typedef boost::property_map<CGALMeshType, CGAL::edge_is_feature_t>::type CGALEIFMap;
    typedef CGALKernalType::Vector_3 CGALVectorType;

    // Make tmp copy, because PMP::corefine_and_compute_intersection does changes on geometry.
    CGALMeshType tmp_polyhedron = mCGALMesh;
    auto p_tmp_cube = CuboidModeler::MakeCuboid(rLowerBound, rUpperBound);

    // Compute intersection surface mesh
    CGALMeshType intersection_mesh{};
    bool valid_intersection = false;

    try {
      valid_intersection = PMP::corefine_and_compute_intersection(tmp_polyhedron, *p_tmp_cube, intersection_mesh);
    }
    catch(const std::exception& exc) {
      std::cerr << "CGAL :: CGAL BrepOperator :: ComputeIntersectionMesh :: No Valid Intersetion!" << std::endl;
      return nullptr;
    }
    if(!valid_intersection){
      std::cerr << "CGAL :: CGAL BrepOperator :: ComputeIntersectionMesh :: No Valid Intersetion!" << std::endl;
      return nullptr;
    }

    const double volume_cube = PMP::volume(*p_tmp_cube);
    const double volume_intersection_surface_mesh = PMP::volume(intersection_mesh);

    //Only consider intersections, which are larger than 0.1% with respect to the original element.
    if( volume_intersection_surface_mesh/volume_cube < 1e-3){
      // std::cout << "Warning :: Intersection neglected! Intersection Polyhedron To Cube Volume Ratio: "
      //   << volume_intersection_surface_mesh/volume_cube << std::endl;

      return nullptr;
    }

    // Construct ptr to SurfaceMesh
    std::unique_ptr<CGALMeshType> p_refinend_intersection_mesh = std::make_unique<CGALMeshType>();

    // Remesh intersected domain until minimum number of boundary triangles is reached.
    double edge_length = rParam.InitialTriangleEdgeLength();
    int iteration_count = 0;
    while (p_refinend_intersection_mesh->number_of_faces() < rParam.MinimumNumberOfTriangles() && iteration_count < 10){
      p_refinend_intersection_mesh->clear();
      CGAL::copy_face_graph(intersection_mesh,*p_refinend_intersection_mesh);

      // Element::CGALPositionType positions = p_refinend_intersection_mesh->points();
      // for (auto vi = p_refinend_intersection_mesh->vertices_begin(); vi != p_refinend_intersection_mesh->vertices_end(); ++vi)
      // {
      //     double x = (positions[*vi].x() - rParam.PointA()[0])/std::abs(rParam.PointA()[0] - rParam.PointB()[0]);
      //     double y = (positions[*vi].y() - rParam.PointA()[1])/std::abs(rParam.PointA()[1] - rParam.PointB()[1]);
      //     double z = (positions[*vi].z() - rParam.PointA()[2])/std::abs(rParam.PointA()[2] - rParam.PointB()[2]);
      //     CGALPointType p(x,y,z);
      //     positions[*vi] = p;
      // }

      CGALEIFMap eif = CGAL::get(CGAL::edge_is_feature, *p_refinend_intersection_mesh);
      PMP::detect_sharp_edges(*p_refinend_intersection_mesh, 10, eif);

      unsigned int nb_iterations = 3;
      try {
        PMP::isotropic_remeshing(faces(*p_refinend_intersection_mesh), edge_length,
          *p_refinend_intersection_mesh,  PMP::parameters::number_of_iterations(nb_iterations).use_safety_constraints(false) // authorize all moves
            .edge_is_constrained_map(eif).number_of_relaxation_steps(1));
      }
      catch(const std::exception& exc) {
          if( !CGAL::is_closed(*p_refinend_intersection_mesh) ){
            std::cout << "Remeshing Exception: " << exc.what() << std::endl;
            std::cout << "Mesh is not closed. Knot spans will be neglected" << std::endl;
            return nullptr;
          }
      }
      edge_length = edge_length * 0.95*std::sqrt((double)p_refinend_intersection_mesh->number_of_faces() / (double) rParam.MinimumNumberOfTriangles()); // Todo: Make this better!!!
      iteration_count++;
    }

    // Element::CGALPositionType positions = p_refinend_intersection_mesh->points();
    // for (auto vi = p_refinend_intersection_mesh->vertices_begin(); vi != p_refinend_intersection_mesh->vertices_end(); ++vi)
    // {
    //     double x = (positions[*vi].x() * std::abs(rParam.PointA()[0] - rParam.PointB()[0])) + rParam.PointA()[0];
    //     double y = (positions[*vi].y() * std::abs(rParam.PointA()[1] - rParam.PointB()[1])) + rParam.PointA()[1];
    //     double z = (positions[*vi].z() * std::abs(rParam.PointA()[2] - rParam.PointB()[2])) + rParam.PointA()[2];
    //     CGALPointType p(x,y,z);
    //     positions[*vi] = p;
    // }

    if (p_refinend_intersection_mesh->number_of_faces() < rParam.MinimumNumberOfTriangles() ){
      std::cout << "Warning:: Targeted number of triangles is not reached: "
        << p_refinend_intersection_mesh->number_of_faces() << " / " <<  rParam.MinimumNumberOfTriangles() << std::endl;
    }

    // Return trimmed domain
    return std::make_unique<CGALTrimmedDomain>(p_refinend_intersection_mesh, rLowerBound, rUpperBound);
}

} // End namespace cgal
} // End namespace tibra

