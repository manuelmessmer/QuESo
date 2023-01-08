// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <fstream>
#include <thread>
#include <omp.h>

//// Project includes
#include "TIBRA_main.h"
#include "quadrature/single_element.h"
#include "quadrature/moment_fitting_utilities.h"
#include "quadrature/multiple_elements.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace tibra {

void TIBRA::Run(){

  typedef BRepOperatorBase::IntersectionStatus IntersectionStatus;

  // Get extreme points of bounding box
  const auto lower_bound = mParameters.LowerBound();
  const auto upper_bound = mParameters.UpperBound();

  const IndexType number_elements_x = mParameters.NumberOfElements()[0];
  const IndexType number_elements_y = mParameters.NumberOfElements()[1];
  const IndexType number_elements_z = mParameters.NumberOfElements()[2];

  // Obtain discretization of background mesh.
  const double delta_x = std::abs(lower_bound[0] - upper_bound[0])/(number_elements_x);
  const double delta_y = std::abs(lower_bound[1] - upper_bound[1])/(number_elements_y);
  const double delta_z = std::abs(lower_bound[2] - upper_bound[2])/(number_elements_z);

  const IndexType global_number_of_elements = number_elements_x * number_elements_y * number_elements_z;
  mpElementContainer->reserve(global_number_of_elements);
  IndexType num = 0;

  // Time Variables
  double et_check_intersect = 0.0;
  double et_compute_intersection = 0.0;
  double et_moment_fitting = 0.0;

  #pragma omp parallel for reduction(+ : et_compute_intersection) reduction(+ : et_check_intersect) reduction(+ : et_moment_fitting) schedule(dynamic)
  for( IndexType i = 0; i < global_number_of_elements; ++i){
    // Unroll 'for' loop to enable better parallelization.
    // First walk along rows (x), then columns (y) then into depths (z).
    const IndexType index_in_row_column_plane = i % (number_elements_x*number_elements_y);
    const IndexType xx = index_in_row_column_plane % number_elements_x; // row
    const IndexType yy = index_in_row_column_plane / number_elements_x; // column
    const IndexType zz = i / (number_elements_x*number_elements_y);     // depth

    // Construct bounding box for each element
    const PointType el_lower_bound{lower_bound[0] + (xx)*delta_x, lower_bound[1] + (yy)*delta_y, lower_bound[2] + (zz)*delta_z};
    const PointType el_upper_bound{lower_bound[0] + (xx+1)*delta_x, lower_bound[1] + (yy+1)*delta_y, lower_bound[2] + (zz+1)*delta_z };

    // Map points into parametric/local space
    const auto el_lower_bound_param = Mapping::GlobalToParam(el_lower_bound, lower_bound, upper_bound);
    const auto el_upper_bound_param = Mapping::GlobalToParam(el_upper_bound, lower_bound, upper_bound);

    // Construct element and check status
    std::shared_ptr<Element> new_element = std::make_shared<Element>(i+1, el_lower_bound_param, el_upper_bound_param, mParameters);

    // Check intersection status
    IntersectionStatus status{};
    if( mParameters.Get<bool>("embedding_flag") ){
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

      auto p_trimmed_domain = mpBRepOperator->GetTrimmedDomain(el_lower_bound, el_upper_bound);
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

        if( new_element->GetIntegrationPoints().size() == 0 ){
          valid_element = false;
        }
      }
    }
    else if( status == IntersectionStatus::Inside){
      // Get standard gauss legendre points
      if( !mParameters.GGQRuleIsUsed() ){
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

  if( mParameters.GGQRuleIsUsed() ){
    #pragma omp single
    MultipleElements::AssembleIPs(*mpElementContainer, mParameters);
  }

  // Average time spend for each task
  if( mParameters.EchoLevel() > 1 ){
    const IndexType num_procs = std::thread::hardware_concurrency();
    std::cout << "#########################################\n";
    std::cout << "Elapsed times of individual tasks: \n";
    std::cout << "Detection of Trimmed Elements: --- " << et_check_intersect / ((double) num_procs) << '\n';
    std::cout << "Compute Intersection: ------------ " << et_compute_intersection / ((double) num_procs) << "\n";
    std::cout << "Moment fitting: ------------------ " << et_moment_fitting / ((double) num_procs) << "\n";
    std::cout << "#########################################\n";
  }

}

} // End namespace tibra
