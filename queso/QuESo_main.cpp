// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <fstream>
#include <thread>
#include <omp.h>

//// Project includes
#include "QuESo_main.h"
#include "io/io_utilities.h"
#include "utilities/mesh_utilities.h"
#include "embedding/brep_operator_factory.h"
#include "quadrature/single_element.h"
#include "quadrature/trimmed_element.h"
#include "quadrature/multiple_elements.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {

void QuESo::Run()
{
    Timer timer{};
    QuESo_INFO_IF(mParameters.EchoLevel() > 0) << "\nQuESo ------------------------------------------ START" << std::endl;

    // Allocate element/knotspans container
    mpElementContainer = MakeUnique<ElementContainer>(mParameters);

    // Loop over all conditions
    for( const auto& p_condition :  mParameters.GetConditions() ){
        const std::string& r_filename = p_condition->GetFilename();
        Unique<TriangleMesh> p_new_mesh = MakeUnique<TriangleMesh>();
        IO::ReadMeshFromSTL(*p_new_mesh, r_filename.c_str());
        // Condition owns triangle mesh.
        mConditions.push_back( ConditionFactory::New(p_condition, p_new_mesh) );
        mpBrepOperatorsBC.push_back( MakeUnique<BRepOperator>(mConditions.back()->GetTriangleMesh(), mParameters) );
    }

    // Read geometry
    double volume_brep = 0.0;
    if( mParameters.Get<bool>("embedding_flag") ) {
        // Read mesh
        const auto& r_filename = mParameters.Get<std::string>("input_filename");
        IO::ReadMeshFromSTL(mTriangleMesh, r_filename.c_str());
        // Write Surface Mesh to vtk file if eco_level > 0
        if( mParameters.EchoLevel() > 0){
            IO::WriteMeshToVTK(mTriangleMesh, "output/geometry.vtk", true);
        }
        // Construct BRepOperator
        mpBRepOperator = BRepOperatorFactory::New(mTriangleMesh, mParameters);

        // Compute volume
        volume_brep = MeshUtilities::VolumeOMP(mTriangleMesh);

        QuESo_INFO_IF(mParameters.EchoLevel() > 0) << "Read file: '" << r_filename << "'\n";
        QuESo_INFO_IF(mParameters.EchoLevel() > 0) << "Volume of B-Rep model: " << volume_brep << '\n';
    }

    // Start computation
    Compute();

    // Count number of trimmed elements
    SizeType number_of_trimmed_elements = 0;
    std::for_each(mpElementContainer->begin(), mpElementContainer->end(), [&number_of_trimmed_elements] (auto& el_it)
        { if( el_it->IsTrimmed() ) { number_of_trimmed_elements++; } });

    if( mParameters.EchoLevel() > 0) {
        // Write vtk files (binary = true)
        IO::WriteElementsToVTK(*mpElementContainer, "output/elements.vtk", true);
        IO::WritePointsToVTK(*mpElementContainer, "All", "output/integration_points.vtk", true);

        for( const auto& r_condition : mConditions ){
            std::string bc_filename = "output/BC_" + std::to_string(r_condition->GetId()) + ".stl";
            IO::WriteMeshToSTL(r_condition->GetConformingMesh(), bc_filename.c_str(), true);
        }

        QuESo_INFO << "Number of active elements: " << mpElementContainer->size() << std::endl;
        QuESo_INFO << "Number of trimmed elements: " << number_of_trimmed_elements << std::endl;

        if( mParameters.EchoLevel() > 1 ) {
            const double volume_ips = mpElementContainer->GetVolumeOfAllIPs();
            QuESo_INFO << "The computed quadrature represents " << volume_ips/volume_brep * 100
                << "% of the volume of the BRep model.\n";
        }

        QuESo_INFO << "Elapsed time: " << timer.Measure() << std::endl;
        QuESo_INFO << "QuESo ------------------------------------------- END\n" << std::endl;
    }
}

void QuESo::Compute(){
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
            const auto bounding_box_xyz = mMapper.GetBoundingBoxXYZFromIndex(index);
            const auto bounding_box_uvw = mMapper.GetBoundingBoxUVWFromIndex(index);

            // Construct element and check status
            Shared<Element> new_element = MakeShared<Element>(index+1, bounding_box_xyz, bounding_box_uvw, mParameters);
            bool valid_element = false;

            // Distinguish between trimmed and non-trimmed elements.
            if( status == IntersectionStatus::Trimmed) {
                new_element->SetIsTrimmed(true);
                Timer timer_compute_intersection{};
                auto p_trimmed_domain = mpBRepOperator->pGetTrimmedDomain(bounding_box_xyz.first, bounding_box_xyz.second);
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

    // Treat conditions
    #pragma omp parallel for
    for( int index = 0; index < static_cast<int>(global_number_of_elements); ++index ) {
        for( IndexType i = 0; i < mConditions.size(); ++i ){
            const auto bounding_box_xyz = mMapper.GetBoundingBoxXYZFromIndex(index);
            const auto p_new_mesh = mpBrepOperatorsBC[i]->pClipTriangleMeshUnique(bounding_box_xyz.first, bounding_box_xyz.second);
            if( p_new_mesh->NumOfTriangles() > 0 ){
                #pragma omp critical
                mConditions[i]->AddToConformingMesh(*p_new_mesh);
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
        QuESo_INFO << "Elapsed times of individual tasks -------------- \n";
        QuESo_INFO << "Detection of trimmed elements: --- " << et_check_intersect / ((double) num_procs) << '\n';
        QuESo_INFO << "Compute intersection: ------------ " << et_compute_intersection / ((double) num_procs) << "\n";
        QuESo_INFO << "Moment fitting: ------------------ " << et_moment_fitting / ((double) num_procs) << "\n";
        QuESo_INFO << "------------------------------------------------ \n";
    }

}

} // End namespace queso
