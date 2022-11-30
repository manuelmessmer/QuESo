// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "TIBRA_main.hpp"
#include "quadrature/single_element.h"
#include "quadrature/moment_fitting_utilities.h"
#include "quadrature/multiple_elements.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

//// External includes
#include <fstream>
#include <thread>
#include <omp.h>

//TODO: Put enums inside outside test and IntegrationMethod into paramters!
void TIBRA::Run(){

  typedef BRepOperatorBase::IntersectionStatus IntersectionStatus;

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
    const PointType cube_point_A{point_a[0] + (xx)*delta_x, point_a[1] + (yy)*delta_y, point_a[2] + (zz)*delta_z};
    const PointType cube_point_B{point_a[0] + (xx+1)*delta_x, point_a[1] + (yy+1)*delta_y, point_a[2] + (zz+1)*delta_z };

    // Map points into parametric/local space
    const auto cube_point_A_param = MappingUtilities::FromGlobalToLocalSpace(cube_point_A, mParameters.PointA(), mParameters.PointB());
    const auto cube_point_B_param = MappingUtilities::FromGlobalToLocalSpace(cube_point_B, mParameters.PointA(), mParameters.PointB());

    // Construct element and check status
    std::shared_ptr<Element> new_element = std::make_shared<Element>(i+1, cube_point_A_param, cube_point_B_param, mParameters);

    // Check intersection status
    IntersectionStatus status{};
    if( mEmbeddingFlag ){
      auto t_begin_di = std::chrono::high_resolution_clock().now();
      status = static_cast<IntersectionStatus>(mpBRepOperator->GetIntersectionState(*new_element));
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
      status = IntersectionStatus::Inside;
    }
    bool valid_element = false;
    // Distinguish between trimmed and non-trimmed elements.
    if( status == IntersectionStatus::Trimmed) {
      new_element->SetIsTrimmed(true);
      auto t_begin_ci = std::chrono::high_resolution_clock().now();

      auto p_trimmed_domain = mpBRepOperator->GetTrimmedDomain(cube_point_A, cube_point_B, mParameters);
      if( p_trimmed_domain ){
        new_element->pSetTrimmedDomain(p_trimmed_domain);
        valid_element = true;
      }
      auto t_end_ci = std::chrono::high_resolution_clock().now();
      std::chrono::duration<double> t_delta_ci = (t_end_ci - t_begin_ci);
      et_compute_intersection += t_delta_ci.count();

      // If valid solve moment fitting equation
      if( valid_element ){
        auto t_begin_mf = std::chrono::high_resolution_clock().now();
        MomentFitting::CreateIntegrationPointsTrimmed(*new_element, mParameters);
        auto t_end_mf = std::chrono::high_resolution_clock().now();
        std::chrono::duration<double> t_delta_mf = (t_end_mf - t_begin_mf);
        et_moment_fitting += t_delta_mf.count();

        if( new_element->GetIntegrationPointsTrimmed().size() == 0 ){
          valid_element = false;
        }
      }
    }
    else if( status == IntersectionStatus::Inside){
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
