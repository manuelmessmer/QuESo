// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <fstream>
#include <thread>
#include <omp.h>

//// Project includes
#include "TIBRA_main.h"
#include "quadrature/single_element.h"
#include "quadrature/trimmed_element.h"
#include "quadrature/multiple_elements.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace tibra {

void TIBRA::Run(){

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
      Timer timer_check_intersect{};
      status = static_cast<IntersectionStatus>(mpBRepOperator->GetIntersectionState(*new_element));
      et_check_intersect += timer_check_intersect.Measure();
    }
    else { // If flag is false, consider all knotspans/ elements as inside
      status = IntersectionStatus::Inside;
    }
    bool valid_element = false;
    // Distinguish between trimmed and non-trimmed elements.
    if( status == IntersectionStatus::Trimmed) {
      new_element->SetIsTrimmed(true);
      Timer timer_compute_intersection{};
      auto p_trimmed_domain = mpBRepOperator->GetTrimmedDomain(el_lower_bound, el_upper_bound);
      if( p_trimmed_domain ){
        new_element->pSetTrimmedDomain(p_trimmed_domain);
        valid_element = true;
      }
      et_compute_intersection += timer_compute_intersection.Measure();

      // If valid solve moment fitting equation
      if( valid_element ){
        Timer timer_moment_fitting{};
        QuadratureTrimmedElement::AssembleIPs(*new_element, mParameters);
        et_moment_fitting += timer_moment_fitting.Measure();

        if( new_element->GetIntegrationPoints().size() == 0 ){
          valid_element = false;
          // if( mParameters.EchoLevel() > 3 ){
          //   auto mesh = new_element->pGetTrimmedDomain()->GetTriangleMesh();
          //   TIBRA_INFO << "Warning :: Moment fitting of trimmed element: " << new_element->GetId() << " with volume: "
          //     << MeshUtilities::Volume(mesh) << " failed! \n";

          //   std::string name = "failed_element_" + std::to_string(new_element->GetId()) + ".stl";
          //   IO::WriteMeshToSTL(mesh, name.c_str(), true);
          // }
        }
      }
    }
    else if( status == IntersectionStatus::Inside){
      // Get standard gauss legendre points
      if( !mParameters.GGQRuleIsUsed() ){
        QuadratureSingleElement::AssembleIPs(*new_element, mParameters);
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
    QuadratureMultipleElements::AssembleIPs(*mpElementContainer, mParameters);
  }

  // Average time spend for each task
  if( mParameters.EchoLevel() > 1 ){
    const IndexType num_procs = std::thread::hardware_concurrency();
    TIBRA_INFO << "Elapsed times of individual tasks -------------- \n";
    TIBRA_INFO << "Detection of trimmed elements: --- " << et_check_intersect / ((double) num_procs) << '\n';
    TIBRA_INFO << "Compute intersection: ------------ " << et_compute_intersection / ((double) num_procs) << "\n";
    TIBRA_INFO << "Moment fitting: ------------------ " << et_moment_fitting / ((double) num_procs) << "\n";
    TIBRA_INFO << "------------------------------------------------ \n";
  }

}

} // End namespace tibra
