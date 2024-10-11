//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// STL includes
#include <fstream>
#include <thread>
#include <omp.h>

//// Project includes
#include "queso/QuESo_main.h"
#include "queso/io/io_utilities.h"
#include "queso/utilities/mesh_utilities.h"
#include "queso/embedding/brep_operator.h"
#include "queso/quadrature/single_element.hpp"
#include "queso/quadrature/trimmed_element.hpp"
#include "queso/quadrature/multiple_elements.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {

void QuESo::Run()
{
    Timer timer{};
    const auto& r_general_settings = mSettings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    QuESo_INFO_IF(echo_level > 0) << "QuESo: Run -------------------------------------- START" << std::endl;

    const auto& r_filename = r_general_settings.GetValue<std::string>(GeneralSettings::input_filename);
    QuESo_INFO_IF(echo_level > 0) << " o Read file: '" << r_filename << "'\n";

    double volume_brep = 0.0;
    // Compute volume
    volume_brep = MeshUtilities::VolumeOMP(*mpTriangleMesh);
    QuESo_INFO_IF(echo_level > 0) << " o Volume of STL model: " << volume_brep << '\n';
    // Construct BRepOperator
    mpBRepOperator = MakeUnique<BRepOperator>(*mpTriangleMesh);
    // Allocate background grid
    mpBackgroundGrid = MakeUnique<BackgroundGridType>(mSettings);

    // Start computation
    std::array<double, 5> elapsed_times = Compute();

    // Count number of trimmed elements
    SizeType number_of_trimmed_elements = 0;
    std::for_each(mpBackgroundGrid->ElementsBegin(), mpBackgroundGrid->ElementsEnd(), [&number_of_trimmed_elements] (auto& el_it)
        { if( el_it.IsTrimmed() ) { number_of_trimmed_elements++; } });

    if( echo_level > 0) {
        QuESo_INFO << " o Number of active elements: " << mpBackgroundGrid->NumberOfActiveElements() << std::endl;
        QuESo_INFO << " o Number of trimmed elements: " << number_of_trimmed_elements << std::endl;

        if( echo_level > 1 ) {
            const double volume_ips = mpBackgroundGrid->GetVolumeOfAllIPs();
            QuESo_INFO << " o The computed quadrature represents " << volume_ips/volume_brep * 100.0
                << "%\n   of the volume of the STL model.\n";

            // Average time spent for each task
            QuESo_INFO << " o Elapsed time (total): " << timer.Measure() << " sec\n";
            QuESo_INFO << " o Elapsed times of individual tasks:\n";
            QuESo_INFO << "   - Classification of elements:    " << elapsed_times[0] << " sec\n";
            QuESo_INFO << "   - Computation of intersections:  " << elapsed_times[1] << " sec\n";
            QuESo_INFO << "   - Solution of moment fitting eq: " << elapsed_times[2] << " sec\n";
            QuESo_INFO << "   - Construction of GGQ rules:     " << elapsed_times[3] << " sec\n";
            QuESo_INFO << "   - Processing of BC's             " << elapsed_times[4] << " sec\n";
        } else {
            QuESo_INFO << " o Elapsed time: " << timer.Measure() << " sec\n";
        }
        QuESo_INFO << "QuESo: Run ---------------------------------------- END\n" << std::endl;
    }

    if( r_general_settings.GetValue<bool>(GeneralSettings::write_output_to_file) ) {
        Timer timer_output{};
        const std::string output_directory_name = r_general_settings.GetValue<std::string>(GeneralSettings::output_directory_name);
        if( echo_level > 0) {
            QuESo_INFO << "QuESo: Write output to file --------------------- START\n";
            QuESo_INFO << " o Output directory: '" << output_directory_name << "'\n";
        }
        // Write vtk files (binary = true)
        IO::WriteMeshToVTK(*mpTriangleMesh, (output_directory_name + "/geometry.vtk" ).c_str(), true);
        IO::WriteElementsToVTK(*mpBackgroundGrid, (output_directory_name + "/elements.vtk").c_str(), true);
        IO::WritePointsToVTK(*mpBackgroundGrid, "All", (output_directory_name + "/integration_points.vtk").c_str(), true);
        IndexType cond_index = 0;
        for( const auto& p_condition : mpBackgroundGrid->GetConditions() ){
            const std::string bc_filename = output_directory_name + '/'
                + p_condition->GetSettings().GetValue<std::string>(ConditionSettings::condition_type)
                + '_' + std::to_string(++cond_index) + ".stl";
            IO::WriteConditionToSTL(*p_condition,  bc_filename.c_str(), true);
        }
        if( echo_level > 0) {
            QuESo_INFO << " o Elapsed time: " << timer_output.Measure() << " sec\n";
            QuESo_INFO << "QuESo: Write output to file ----------------------- End\n" << std::endl;
        }
    }

}

std::array<double,5> QuESo::Compute(){

    // Reserve element container
    const IndexType global_number_of_elements = mGridIndexer.NumberOfElements();
    mpBackgroundGrid->ReserveElements(global_number_of_elements);

    // Time Variables
    double et_check_intersect = 0.0;
    double et_compute_intersection = 0.0;
    double et_moment_fitting = 0.0;

    // Get neccessary settings
    const IntegrationMethod integration_method = mSettings[MainSettings::non_trimmed_quadrature_rule_settings]
        .GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method);
    const bool ggq_rule_ise_used =  static_cast<int>(integration_method) >= 3;
    const auto& r_trimmed_quad_rule_settings = mSettings[MainSettings::trimmed_quadrature_rule_settings];
    const double min_vol_element_ratio = std::max(r_trimmed_quad_rule_settings.GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-10);
    const IndexType num_boundary_triangles = r_trimmed_quad_rule_settings.GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles);
    const double moment_fitting_residual = r_trimmed_quad_rule_settings.GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
    const bool neglect_elements_if_stl_is_flawed = r_trimmed_quad_rule_settings.GetValue<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed);
    const auto& r_grid_settings = mSettings[MainSettings::background_grid_settings];
    const Vector3i polynomial_order = r_grid_settings.GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);
    const Vector3i number_of_elements = r_grid_settings.GetValue<Vector3i>(BackgroundGridSettings::number_of_elements);

    const IndexType echo_level = mSettings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level);

    // Classify all elements.
    Unique<BRepOperator::StatusVectorType> p_classifications = nullptr;

    Timer timer_check_intersect{};
    p_classifications = mpBRepOperator->pGetElementClassifications(mSettings);
    et_check_intersect = timer_check_intersect.Measure();

    #pragma omp parallel for reduction(+ : et_compute_intersection) reduction(+ : et_moment_fitting) schedule(dynamic)
    for( int index = 0; index < static_cast<int>(global_number_of_elements); ++index) {
        // Check classification status
        const IntersectionState status = (*p_classifications)[index];

        if( status == IntersectionState::inside || status == IntersectionState::trimmed ) {
            // Get bounding box of element
            const auto bounding_box_xyz = mGridIndexer.GetBoundingBoxXYZFromIndex(index);
            const auto bounding_box_uvw = mGridIndexer.GetBoundingBoxUVWFromIndex(index);

            // Construct element and check status:
            Unique<ElementType> new_element = MakeUnique<ElementType>(index+1, bounding_box_xyz, bounding_box_uvw);
            bool valid_element = false;

            // Distinguish between trimmed and non-trimmed elements.
            if( status == IntersectionState::trimmed) {
                new_element->SetIsTrimmed(true);
                Timer timer_compute_intersection{};
                auto p_trimmed_domain = mpBRepOperator->pGetTrimmedDomain(bounding_box_xyz.first, bounding_box_xyz.second,
                                                                          min_vol_element_ratio, num_boundary_triangles, neglect_elements_if_stl_is_flawed);
                if( p_trimmed_domain ){
                    new_element->pSetTrimmedDomain(p_trimmed_domain);
                    valid_element = true;
                }
                et_compute_intersection += timer_compute_intersection.Measure();

                // If valid solve moment fitting equation
                if( valid_element ){
                    Timer timer_moment_fitting{};
                    QuadratureTrimmedElement<ElementType>::AssembleIPs(*new_element, polynomial_order, moment_fitting_residual, echo_level);
                    et_moment_fitting += timer_moment_fitting.Measure();

                    if( new_element->GetIntegrationPoints().size() == 0 ){
                        valid_element = false;
                    }
                }
            }
            else if( status == IntersectionState::inside){
                // Get standard gauss legendre points
                if( !ggq_rule_ise_used ){
                    QuadratureSingleElement<ElementType>::AssembleIPs(*new_element, polynomial_order, integration_method);
                }
                valid_element = true;
            }

            if( valid_element ){
                #pragma omp critical // TODO: improve this.
                mpBackgroundGrid->AddElement(new_element); // After this new_element is a null_ptr. Is std::moved to container.
            }
        }

    }

    double et_ggq_rules = 0.0;
    if( ggq_rule_ise_used ){
        Timer timer_ggq_rules;
        QuadratureMultipleElements<ElementType>::AssembleIPs(*mpBackgroundGrid, number_of_elements, polynomial_order, integration_method);
        et_ggq_rules = timer_ggq_rules.Measure();
    }

    // Treat conditions
    Timer timer_conditions;

    const auto& r_conditions_settings_list =  mSettings[MainSettings::conditions_settings_list];
    const IndexType number_of_conditions = r_conditions_settings_list.NumberOfSubDictionaries();
    mpBackgroundGrid->ReserveConditions(number_of_conditions);
    for( IndexType i_c = 0; i_c < number_of_conditions; ++i_c ) {
        const SettingsBaseType& r_condition_settings = r_conditions_settings_list[i_c];

        const std::string& r_filename = r_condition_settings.GetValue<std::string>(ConditionSettings::input_filename);
        Unique<TriangleMeshInterface> p_new_mesh = MakeUnique<TriangleMesh>();
        IO::ReadMeshFromSTL(*p_new_mesh, r_filename.c_str());

        Unique<ConditionType> p_new_condition = MakeUnique<ConditionType>(p_new_mesh, r_condition_settings);
        BRepOperator brep_operator_c(p_new_condition->GetTriangleMesh());

        #pragma omp parallel for
        for( int index = 0; index < static_cast<int>(global_number_of_elements); ++index ) {

            const auto bounding_box_xyz = mGridIndexer.GetBoundingBoxXYZFromIndex(index);
            auto p_new_mesh = brep_operator_c.pClipTriangleMeshUnique(bounding_box_xyz.first, bounding_box_xyz.second);
            if( p_new_mesh->NumOfTriangles() > 0 ) {
                const auto p_el = mpBackgroundGrid->pGetElement(index+1);
                auto p_new_segment = MakeUnique<ConditionType::ConditionSegmentType>(p_el, p_new_mesh);
                #pragma omp critical
                p_new_condition->AddSegment(p_new_segment);
            }
        }
        mpBackgroundGrid->AddCondition(p_new_condition);
    }
    double et_conditions = (mpBackgroundGrid->NumberOfConditions() > 0) ? timer_conditions.Measure() : 0.0;

    IndexType num_threads = 1;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
    return {et_check_intersect, et_compute_intersection / ((double) num_threads),
        et_moment_fitting / ((double) num_threads), et_ggq_rules, et_conditions};

}

void QuESo::Check() const {
    // Check if bounding box fully contains the triangle mesh.
    if( mSettings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level) > 0 ){
        PointType lower_bound = mSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
        PointType upper_bound = mSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz);
        auto bb_mesh = MeshUtilities::BoundingBox(*mpTriangleMesh);
        if( lower_bound[0] > bb_mesh.first[0]  ||
            lower_bound[1] > bb_mesh.first[1]  ||
            lower_bound[2] > bb_mesh.first[2]  ||
            upper_bound[0] < bb_mesh.second[0] ||
            upper_bound[1] < bb_mesh.second[1] ||
            upper_bound[2] < bb_mesh.second[2] ) {
                QuESo_INFO << "Warning :: The given bounding box: 'lower_bound_xyz' : " << lower_bound
                    << ", 'upper_bound_xyz:' " << upper_bound << " does not fully contain the bounding box of "
                    << "the input STL: 'lower_bound_xyz' : " << bb_mesh.first << ", 'upper_bound_xyz:' "
                    << bb_mesh.second << '\n';
        }
    }

}

} // End namespace queso
