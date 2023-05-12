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
    const IndexType number_elements_x = mParameters.NumberOfElements()[0];
    const IndexType number_elements_y = mParameters.NumberOfElements()[1];
    const IndexType number_elements_z = mParameters.NumberOfElements()[2];

    const IndexType global_number_of_elements = number_elements_x * number_elements_y * number_elements_z;
    mpElementContainer->reserve(global_number_of_elements);

    // Time Variables
    double et_check_intersect = 0.0;
    double et_compute_intersection = 0.0;
    double et_moment_fitting = 0.0;

    // Classify all elements.
    Unique<BRepOperatorBase::StatusVectorType> p_classifications = nullptr;
    if( mParameters.Get<bool>("embedding_flag") ){
        Timer timer_check_intersect{};
        p_classifications = mpBRepOperator->pGetElementClassifications();
        et_check_intersect += timer_check_intersect.Measure();
    }

    #pragma omp parallel for reduction(+ : et_compute_intersection) reduction(+ : et_check_intersect) reduction(+ : et_moment_fitting) schedule(dynamic)
    for( int index = 0; index < static_cast<int>(global_number_of_elements); ++index) {
        // Check classification status
        IntersectionStatus status{};
        if( mParameters.Get<bool>("embedding_flag") ){
            status = (*p_classifications)[index];
        }
        else { // If flag is false, consider all knotspans/ elements as inside
            status = IntersectionStatus::Inside;
        }

        if( status == IntersectionStatus::Inside || status == IntersectionStatus::Trimmed ) {
            // Get bounding box of element
            const auto bounding_box = mMapper.GetBoundingBoxFromIndex(index);

            // Map points into parametric/local space
            const auto bb_lower_bound_param = mMapper.PointFromGlobalToParam(bounding_box.first);
            const auto bb_upper_bound_param = mMapper.PointFromGlobalToParam(bounding_box.second);

            // Construct element and check status
            Shared<Element> new_element = MakeShared<Element>(index+1, bb_lower_bound_param, bb_upper_bound_param, mParameters);
            bool valid_element = false;

            // Distinguish between trimmed and non-trimmed elements.
            if( status == IntersectionStatus::Trimmed) {
                new_element->SetIsTrimmed(true);
                Timer timer_compute_intersection{};
                auto p_trimmed_domain = mpBRepOperator->pGetTrimmedDomain(bounding_box.first, bounding_box.second);
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
                    }
                }
            }
            else if( status == IntersectionStatus::Inside){
                // Get standard gauss legendre points
                if( !mParameters.GGQRuleIsUsed() ){
                    QuadratureSingleElement::AssembleIPs(*new_element, mParameters);
                }
                valid_element = true;
            }

            if( valid_element ){
                #pragma omp critical
                mpElementContainer->AddElement(new_element); // After this new_element is a null_ptr. Is std::moved to container.
            }
        }
    }

    if( mParameters.GGQRuleIsUsed() ){
        #pragma omp single
        QuadratureMultipleElements::AssembleIPs(*mpElementContainer, mParameters);
    }

    // Average time spent for each task
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
