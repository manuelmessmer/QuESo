// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "TIBRA_main.hpp"
#include "quadrature/single_element.h"
#include "quadrature/moment_fitting_utilities.h"
#include "geometries/triangle_3d_3n.h"
#include "quadrature/multiple_elements.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "modeler/modeler.h"

//// External includes
#include <fstream>
#include <thread>
#include <omp.h>

//TODO: Put enums inside outside test and IntegrationMethod into paramters!

void TIBRA::Run(){
  typedef BRepOperatorBase::BoundaryIPVectorType BoundaryIPVectorType;
  typedef BRepOperatorBase::BoundaryIPVectorPtrType BoundaryIPVectorPtrType;

  // Get extreme points of bounding box
  const auto point_a = mParameters.PointA();
  const auto point_b = mParameters.PointB();

  // Obtain discretization of background mesh.
  const double delta_x = std::abs(mParameters.PointA()[0] - mParameters.PointB()[0])/(mParameters.NumberOfElements()[0]);
  const double delta_y = std::abs(mParameters.PointA()[1] - mParameters.PointB()[1])/(mParameters.NumberOfElements()[1]);
  const double delta_z = std::abs(mParameters.PointA()[2] - mParameters.PointB()[2])/(mParameters.NumberOfElements()[2]);

  const int number_elements_u = mParameters.NumberOfElements()[0];
  const int number_elements_v = mParameters.NumberOfElements()[1];
  const int number_elements_w = mParameters.NumberOfElements()[2];

  const int global_number_of_elements = number_elements_u * number_elements_v * number_elements_w;
  mpElementContainer->reserve(global_number_of_elements);
  int num = 0;

  // Time Variables
  double et_check_intersect = 0.0;
  double et_compute_intersection = 0.0;
  double et_moment_fitting = 0.0;

  #pragma omp parallel for reduction(+ : et_compute_intersection) reduction(+ : et_check_intersect) reduction(+ : et_moment_fitting) schedule(dynamic)
  for( int i = 0; i < global_number_of_elements; ++i){
    // Unroll 'for' loop to enable better parallelization.
    // First walk along rows (x), then columns (y) then into depths (z).
    const int index_in_row_column_plane = i % (number_elements_u*number_elements_v);
    const int xx = index_in_row_column_plane % number_elements_u; // row
    const int yy = index_in_row_column_plane / number_elements_u; // column
    const int zz = i / (number_elements_u*number_elements_v);     // depth

    // Construct bounding box for each element
    PointType cube_point_A{point_a[0] + (xx)*delta_x, point_a[1] + (yy)*delta_y, point_a[2] + (zz)*delta_z};
    PointType cube_point_B{point_a[0] + (xx+1)*delta_x, point_a[1] + (yy+1)*delta_y, point_a[2] + (zz+1)*delta_z };

    // Map points into parametric/local space
    const auto cube_point_A_param = MappingUtilities::FromGlobalToLocalSpace(cube_point_A, mParameters.PointA(), mParameters.PointB());
    const auto cube_point_B_param = MappingUtilities::FromGlobalToLocalSpace(cube_point_B, mParameters.PointA(), mParameters.PointB());

    // Construct element and check status
    std::shared_ptr<Element> new_element = std::make_shared<Element>(i+1, cube_point_A_param, cube_point_B_param, mParameters);

    // Get surface mesh of cube
    SurfaceMeshType cube{};
    Modeler::make_cube_3(cube, cube_point_A, cube_point_B);

    // Check intersection status
    IntersectionTest::IntersectionStatus status{};
    if( mEmbeddingFlag ){
      auto t_begin_di = std::chrono::high_resolution_clock().now();
      status = static_cast<IntersectionTest::IntersectionStatus>(mpBRepOperator->GetIntersectionState(*new_element));
      auto t_end_di = std::chrono::high_resolution_clock().now();
      std::chrono::duration<double> t_delta_di = (t_end_di - t_begin_di);
      et_check_intersect += t_delta_di.count();

      // auto t_begin_di_2 = std::chrono::high_resolution_clock().now();
      // auto status_2 = mClassifier->GetIntersectionState( *new_element );
      // auto t_end_di_2 = std::chrono::high_resolution_clock().now();

      // // if( status != status_2){
      // //   std::cout << "status: " << status << ", " << status_2 << std::endl;
      // //   throw std::runtime_error("ooho");
      // // }
      // std::chrono::duration<double> t_delta_di_2 = (t_end_di_2 - t_begin_di_2);
      // et_check_intersect_2 += t_delta_di_2.count();

    }
    else { // If flag is false, consider all knotspans/ elements as inside
      status = IntersectionTest::Inside;
    }
    bool valid_element = false;
    // Distinguish between trimmed and non-trimmed elements.
    if( status == IntersectionTest::Trimmed) {
      new_element->SetIsTrimmed(true);
      auto t_begin_ci = std::chrono::high_resolution_clock().now();

      auto p_boundary_ips = std::make_unique<BoundaryIPVectorType>();
      valid_element = mpBRepOperator->ComputeBoundaryIps(*new_element, p_boundary_ips, mParameters);
      auto t_end_ci = std::chrono::high_resolution_clock().now();
      std::chrono::duration<double> t_delta_ci = (t_end_ci - t_begin_ci);
      et_compute_intersection += t_delta_ci.count();

      // If valid solve moment fitting equation
      if( valid_element ){
        auto t_begin_mf = std::chrono::high_resolution_clock().now();
        MomentFitting::CreateIntegrationPointsTrimmed(*new_element, p_boundary_ips, mParameters);
        auto t_end_mf = std::chrono::high_resolution_clock().now();
        std::chrono::duration<double> t_delta_mf = (t_end_mf - t_begin_mf);
        et_moment_fitting += t_delta_mf.count();

        if( new_element->GetIntegrationPointsTrimmed().size() == 0 ){
          valid_element = false;
        }
      }
    }
    else if( status == IntersectionTest::Inside){
      // Get standard gauss legendre points
      if( mParameters.IntegrationMethod() <= 2 ){
        SingleElement::AssembleIPs(*new_element, mParameters);
      }

      //ExportVolumeMesh(cube, new_element->GetId());
      valid_element = true;
    }

    if( valid_element ){
      #pragma omp critical
      mpElementContainer->AddElement(new_element); // After this new_element is a null_ptr. Is std::moved to container.
    }
  }

  if( mParameters.IntegrationMethod() >= 3 ){
    #pragma omp single
    MultipleElements::AssembleIPs(*mpElementContainer, mParameters);
  }

  // Average time spend for each task
  if( mParameters.EchoLevel() > 1 ){
    const int num_procs = std::thread::hardware_concurrency();
    std::cout << "#########################################\n";
    std::cout << "Elapsed times of individual tasks: \n";
    std::cout << "Detection of Trimmed Elements: --- " << et_check_intersect / ((double) num_procs) << '\n';
    std::cout << "Compute Intersection: ------------ " << et_compute_intersection / ((double) num_procs) << "\n";
    std::cout << "Moment fitting: ------------------ " << et_moment_fitting / ((double) num_procs) << "\n";
    std::cout << "#########################################\n";
  }

}

// #include <CGAL/Mesh_triangulation_3.h>
// #include <CGAL/Mesh_complex_3_in_triangulation_3.h>
// #include <CGAL/Mesh_criteria_3.h>
// #include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
// #include <CGAL/make_mesh_3.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/IO/output_to_vtu.h>

// void STLEmbedder::ExportVolumeMesh() {
//   // Triangulation
//   // std::string file_name_1 = "output/meshes/vtu/mesh_" + std::to_string(id) + ".vtu";
//   // polygon_mesh_to_vtkUnstructured_(surface_mesh, file_name_1.c_str());

//   typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
//   typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
//   typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,CGAL::Sequential_tag>::type Tr;
//   typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
//   // Criteria
//   typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

//   using namespace CGAL::parameters;
//   // Create domain
//   SurfaceMeshType tmp_surface_mesh{};
//   CGAL::copy_face_graph(mPolyhedronForExport, tmp_surface_mesh);

//   Polyhedron poly{};
//   CGAL::copy_face_graph(mPolyhedronForExport, poly);

//   Mesh_domain domain(poly);
//   domain.detect_features();

//   // Cylinder:
//   //Mesh_criteria criteria( edge_size = 0.1, cell_size=0.1);
//   // Chaumodell:
//   //Mesh_criteria criteria( edge_size = 3.0, cell_size=3.0);
//   // Gehause:
//   Mesh_criteria criteria( edge_size = 1.5, cell_size= 1.5);
//   // Mesh generation
//   C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb(), manifold_with_boundary());
//   // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
//   // Mesh_criteria new_criteria(cell_radius_edge_ratio=0.5, cell_size=0.1);
//   // // Mesh refinement (and make the output manifold)
//   // CGAL::refine_mesh_3(c3t3, domain, new_criteria);

//   std::string file_name_1 = "output/meshes/vtu/physical_mesh.vtu";
//   std::ofstream os_file(file_name_1);

//   CGAL::output_to_vtu(os_file, c3t3, CGAL::IO::ASCII);
//   os_file.close();

//   Tr& tr = c3t3.triangulation();

//   // Mesh criteria (no cell_size set)
//   std::cout << "Number of vertices in compelx: " << c3t3.number_of_vertices_in_complex() << std::endl;
//   for (auto vi = tr.finite_vertices_begin(); vi != tr.finite_vertices_end(); ++vi)
//   {
//     double x = (vi->point().x() - mParameters.PointA()[0])/std::abs(mParameters.PointA()[0] - mParameters.PointB()[0]);
//     //std::cout << "Soll: " << x << std::endl;
//     double y = (vi->point().y() - mParameters.PointA()[1])/std::abs(mParameters.PointA()[1] - mParameters.PointB()[1]);
//     double z = (vi->point().z() - mParameters.PointA()[2])/std::abs(mParameters.PointA()[2] - mParameters.PointB()[2]);
//     CGAL::Weighted_point_3<K> WeightedPoint(x, y, z);

//     vi->point() = WeightedPoint;

//     //std::cout << "Ist: " << vi->point().x() << std::endl;

//   }

//   // Output
//   std::string file_name_2 = "output/meshes/vtu/param_mesh.vtu";
//   std::ofstream os_file2(file_name_2);
//   CGAL::output_to_vtu(os_file2, c3t3, CGAL::IO::ASCII);
//   os_file2.close();

// }

//Remove this again
// #include <CGAL/Polygon_mesh_processing/remesh.h>
// #include <CGAL/Polygon_mesh_processing/detect_features.h>
// //Remove this again
// typedef boost::property_map<CGAL::Surface_mesh<Point_3>, CGAL::edge_is_feature_t>::type EIFMap;
// namespace PMP = CGAL::Polygon_mesh_processing;
// /// Functions
// void STLEmbedder::ExportSTL(){
//     // auto element_it_begin = mNeclectedElements.begin();
//     // for( int i = 0; i < mNeclectedElements.size(); ++i){
//     //   auto element_it = element_it_begin + i;
//     //   if( (*element_it)->IsNeclected() ){
//     //     std::array<double, 3> cube_point_A = MappingUtilities::FromLocalToGlobalSpace(
//     //         (*element_it)->GetLocalLowerPoint(), mParameters.PointA(), mParameters.PointB() );
//     //     std::array<double, 3> cube_point_B = MappingUtilities::FromLocalToGlobalSpace(
//     //         (*element_it)->GetLocalUpperPoint(), mParameters.PointA(), mParameters.PointB() );
//     //     cube_point_A[0] -= 0.1;
//     //     cube_point_A[1] -= 0.1;
//     //     cube_point_A[2] -= 0.1;
//     //     cube_point_B[0] += 0.1;
//     //     cube_point_B[1] += 0.1;
//     //     cube_point_B[2] += 0.1;
//     //     SurfaceMeshType cube;
//     //     CubeModeler::make_cube_3(cube, cube_point_A, cube_point_B);
//     //     cgal::BRepOperator::SubstractElementsFromSurfaceMesh(mPolyhedronForExport, cube, **element_it, mParameters);
//     //   }
//     // }

//     auto el_itr_begin = mpElementContainer->begin();
//     for( int i = 0; i < mpElementContainer->size(); ++i){
//       auto el_itr = el_itr_begin + i;
//       const int id = (*el_itr)->GetId();
//       auto s_mesh = (*el_itr)->GetSurfaceMesh();
//       Element::PositionType positions = s_mesh.points();
//       for (auto vi = s_mesh.vertices_begin(); vi != s_mesh.vertices_end(); ++vi)
//       {
//           double x = (positions[*vi].x() - mParameters.PointA()[0])/std::abs(mParameters.PointA()[0] - mParameters.PointB()[0]);
//           double y = (positions[*vi].y() - mParameters.PointA()[1])/std::abs(mParameters.PointA()[1] - mParameters.PointB()[1]);
//           double z = (positions[*vi].z() - mParameters.PointA()[2])/std::abs(mParameters.PointA()[2] - mParameters.PointB()[2]);
//           Point_3 p(x,y,z);
//           positions[*vi] = p;
//       }
//       std::string name2 = "output/meshes/vtu/" + std::to_string(id) + ".vtu";
//       CGAL::polygon_mesh_to_vtkUnstructured_(s_mesh, name2.c_str());
//     }

    // EIFMap eif = CGAL::get(CGAL::edge_is_feature, mPolyhedronForExport);
    // PMP::detect_sharp_edges(mPolyhedronForExport, 10, eif);

    // unsigned int nb_iterations = 3;
    // CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mPolyhedronForExport), 0.15,
    //               mPolyhedronForExport,  PMP::parameters::number_of_iterations(nb_iterations).use_safety_constraints(false) // authorize all moves
    //                                     .edge_is_constrained_map(eif));

    // std::string fullname= "mesh_physical_space.off";
    // std::ofstream os(fullname.c_str() ) ;
    // os << mPolyhedronForExport ;
    // std::string name = "mesh_physical_space.vtu";
    // polygon_mesh_to_vtkUnstructured_(mPolyhedronForExport, name.c_str());
    // // Map mesh into parameter space
    // Element::PositionType positions = mPolyhedronForExport.points();
    // for (auto vi = mPolyhedronForExport.vertices_begin(); vi != mPolyhedronForExport.vertices_end(); ++vi)
    // {
    //     double x = (positions[*vi].x() - mParameters.PointA()[0])/std::abs(mParameters.PointA()[0] - mParameters.PointB()[0]);
    //     double y = (positions[*vi].y() - mParameters.PointA()[1])/std::abs(mParameters.PointA()[1] - mParameters.PointB()[1]);
    //     double z = (positions[*vi].z() - mParameters.PointA()[2])/std::abs(mParameters.PointA()[2] - mParameters.PointB()[2]);
    //     Point_3 p(x,y,z);
    //     positions[*vi] = p;
    // }
    // std::string fullname2= "mesh_param_space.off";
    // std::ofstream os2(fullname2.c_str() ) ;
    // os2 << mPolyhedronForExport ;
    // std::string name2 = "mesh_param_space.vtu";
    // polygon_mesh_to_vtkUnstructured_(mPolyhedronForExport, name2.c_str());

//}