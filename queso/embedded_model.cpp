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
#include <omp.h>

//// Project includes
#include "queso/embedded_model.h"
#include "queso/io/io_utilities.h"
#include "queso/utilities/mesh_utilities.h"
#include "queso/embedding/brep_operator.h"
#include "queso/quadrature/single_element.hpp"
#include "queso/quadrature/trimmed_element.hpp"
#include "queso/quadrature/multiple_elements.hpp"

namespace queso {

void EmbeddedModel::ComputeVolume(const TriangleMeshInterface& rTriangleMesh){

    CheckIfMeshIsWithinBoundingBox(rTriangleMesh);

    // Set ModelInfo.
    auto& r_model_info = GetModelInfo();
    // EmbeddedGeometryInfo.
    const double volume = MeshUtilities::VolumeOMP(rTriangleMesh);
    r_model_info[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::volume, volume);
    const bool is_closed = MeshUtilities::EstimateQuality(rTriangleMesh) < 1e-10;
    r_model_info[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::is_closed, is_closed);

    /* Begin: Get neccessary settings */
    /// @todo Pass r_settings to all functions.

    const auto& r_settings = GetSettings();
    // MainSettings::general_settings.
    const IndexType echo_level = r_settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level);

    // MainSettings::non_trimmed_quadrature_rule_settings.
    const IntegrationMethod integration_method = r_settings[MainSettings::non_trimmed_quadrature_rule_settings].
        GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method);
    const bool ggq_rule_is_used =  static_cast<int>(integration_method) >= 3;

    // MainSettings::trimmed_quadrature_rule_settings.
    const auto& r_trimmed_quad_rule_settings = r_settings[MainSettings::trimmed_quadrature_rule_settings];
    const double min_vol_element_ratio = std::max<double>(r_trimmed_quad_rule_settings.
        GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-10);
    const IndexType num_boundary_triangles = r_trimmed_quad_rule_settings.
        GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles);
    const double moment_fitting_residual = r_trimmed_quad_rule_settings.
        GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
    const bool neglect_elements_if_stl_is_flawed = r_trimmed_quad_rule_settings.
        GetValue<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed);

    // MainSettings::background_grid_settings.
    const auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    const Vector3i polynomial_order = r_grid_settings.GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);

    /* End: Get neccessary settings */

    // Start timer.
    Timer timer_total{};

    // Reserve element container.
    const IndexType global_number_of_elements = mGridIndexer.NumberOfElements();
    mBackgroundGrid.ReserveElements(global_number_of_elements);

    // Construct BRepOperator.
    BRepOperator brep_operator(rTriangleMesh);

    // Classify all elements.
    Timer timer_check_intersect{};
    Unique<BRepOperator::StatusVectorType> p_classifications = brep_operator.pGetElementClassifications(r_settings);
    auto& r_volume_time_info = r_model_info[MainInfo::elapsed_time_info][ElapsedTimeInfo::volume_time_info];
    r_volume_time_info.SetValue(VolumeTimeInfo::classification_of_elements, timer_check_intersect.Measure());

    /* Begin: Initialize info variables */
    // Variables for performance measurements.
    double et_compute_intersection = 0.0;
    double et_moment_fitting = 0.0;

    // Grid information counters
    SizeType num_active_elements = 0;
    SizeType num_trimmed_elements = 0;
    /* End: Initialize info variables */

    // Loop over all elements.
    #pragma omp parallel for reduction(+ : et_compute_intersection, et_moment_fitting, num_active_elements, num_trimmed_elements) schedule(dynamic)
    for( int for_index = 0; for_index < static_cast<int>(global_number_of_elements); ++for_index) {
		const auto index = static_cast<IndexType>(for_index);
        // Check classification status.
        const IntersectionState status = (*p_classifications)[index];

        if( status == IntersectionState::inside || status == IntersectionState::trimmed ) {
            // Get bounding box of element.
            const auto bounding_box_xyz = mGridIndexer.GetBoundingBoxXYZFromIndex(index);
            const auto bounding_box_uvw = mGridIndexer.GetBoundingBoxUVWFromIndex(index);

            // Construct element.
            Unique<ElementType> p_new_element = MakeUnique<ElementType>(index+1, bounding_box_xyz, bounding_box_uvw);
            bool valid_element_flag = false;

            // Distinguish between trimmed and non-trimmed elements.
            if( status == IntersectionState::trimmed) {
                // Compute intersection between element and embedded geometry.
                Timer timer_compute_intersection{};
                auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(bounding_box_xyz.first, bounding_box_xyz.second,
                                                                        min_vol_element_ratio, num_boundary_triangles, neglect_elements_if_stl_is_flawed);
                if( p_trimmed_domain ){
                    p_new_element->pSetTrimmedDomain(p_trimmed_domain);
                    valid_element_flag = true;
                }
                et_compute_intersection += timer_compute_intersection.Measure();

                // If valid, solve moment fitting equation.
                if( valid_element_flag ){
                    Timer timer_moment_fitting{};
                    QuadratureTrimmedElement<ElementType>::AssembleIPs(*p_new_element, polynomial_order, moment_fitting_residual, echo_level);
                    et_moment_fitting += timer_moment_fitting.Measure();

                    if( p_new_element->GetIntegrationPoints().size() == 0 ){
                        valid_element_flag = false;
                    }
                }
            }
            else { // status == IntersectionState::inside.
                if( !ggq_rule_is_used ){
                    // Get standard gauss legendre points.
                    QuadratureSingleElement<ElementType>::AssembleIPs(*p_new_element, polynomial_order, integration_method);
                }
                valid_element_flag = true;
            }

            if( valid_element_flag ){
                // Update background grid counters.
                ++num_active_elements;
                if( p_new_element->IsTrimmed() ) { ++num_trimmed_elements; }

                // Add p_new_element to the background grid.
                #pragma omp critical
                mBackgroundGrid.AddElement(std::move(p_new_element));
            }
        }
    } // #pragma parallel omp for reduction.

    // Assmble Generalized Gaussian quadrature rules (if enabled).
    double et_ggq_rules = 0.0;
    if( ggq_rule_is_used ){
        // Construct Generalized Gaussian quadrature rules.
        Timer timer_ggq_rules{};
        QuadratureMultipleElements<ElementType>::AssembleIPs(mBackgroundGrid, polynomial_order, integration_method);
        et_ggq_rules = timer_ggq_rules.Measure();
    }

    // Stop timer (total time).
    const double elapsed_time_total = timer_total.Measure();

    /* Begin: Write model to mModeInfo */

    // ElpasedTimeInfo.
    auto& r_elapsed_time_info = r_model_info[MainInfo::elapsed_time_info];
    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::total, elapsed_time_total);
    const double total_time = (r_elapsed_time_info.IsSet(ElapsedTimeInfo::total)) ?
        r_elapsed_time_info.GetValue<double>(ElapsedTimeInfo::total) : 0.0;
    r_elapsed_time_info.SetValue(ElapsedTimeInfo::total, (total_time+elapsed_time_total) );

    // Get num of threads.
    const auto num_threads = static_cast<IndexType>(omp_get_max_threads());
    r_volume_time_info.SetValue(VolumeTimeInfo::computation_of_intersections, et_compute_intersection / static_cast<double>(num_threads) );
    r_volume_time_info.SetValue(VolumeTimeInfo::solution_of_moment_fitting_eqs, et_moment_fitting / static_cast<double>(num_threads) );
    r_volume_time_info.SetValue(VolumeTimeInfo::construction_of_ggq_rules, et_ggq_rules);

    // BackgroundGridInfo.
    r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_active_elements, num_active_elements );
    r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_trimmed_elements, num_trimmed_elements );
    const IndexType num_full_elements = (num_active_elements-num_trimmed_elements);
    r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_full_elements, num_full_elements );
    const IndexType num_inactive_elements =  ( mGridIndexer.NumberOfElements() - num_active_elements);
    r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_inactive_elements, num_inactive_elements );

    // QuadratureInfo.
    double represented_volume = 0.0;
    SizeType tot_num_points_full = 0;
    SizeType tot_num_points_trimmed = 0;
    const auto el_it_ptr_begin = mBackgroundGrid.ElementsBegin();
    #pragma omp parallel for reduction(+ : represented_volume, tot_num_points_full, tot_num_points_trimmed) schedule(guided)
    for( int i = 0; i < static_cast<int>(num_active_elements); ++i ){
        const auto el_ptr = (el_it_ptr_begin + i);
        const double det_j = el_ptr->DetJ();
        const auto& r_points = el_ptr->GetIntegrationPoints();
        represented_volume += std::accumulate(r_points.begin(), r_points.end(), 0.0,
            [det_j](double sum, const auto& point) {
                return sum + point.Weight() * det_j; }
        );
        auto& target_counter = el_ptr->IsTrimmed() ? tot_num_points_trimmed : tot_num_points_full;
        target_counter += r_points.size();
    }
    const SizeType tot_num_points = tot_num_points_trimmed + tot_num_points_full;
    r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::tot_num_points, tot_num_points);
    r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::represented_volume, represented_volume);
    r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::percentage_of_geometry_volume, represented_volume/volume*100.0);
    const double num_of_points_per_full_element = (num_full_elements > 0) ?
        static_cast<double>(tot_num_points_full)/static_cast<double>(num_full_elements) : 0.0;
    r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_full_element, num_of_points_per_full_element);
    const double num_of_points_per_trimmed_element = (num_trimmed_elements > 0) ?
        static_cast<double>(tot_num_points_trimmed)/static_cast<double>(num_trimmed_elements) : 0.0;
    r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_trimmed_element, num_of_points_per_trimmed_element);

    /* End: Write model to mModeInfo */

    // Print some info to console.
    PrintVolumeInfo();
}

void EmbeddedModel::ComputeCondition(const TriangleMeshInterface& rTriangleMesh, const MainDictionaryType& rConditionSettings) {

    CheckIfMeshIsWithinBoundingBox(rTriangleMesh);

    // Start timer.
    Timer timer_condition{};

    // Get model info.
    auto& r_model_info = GetModelInfo();

    /* Begin: Write to r_model_info */

    // Create condition info.
    auto p_new_cond_info = DictionaryFactory<key::MainValuesTypeTag>::Create("ConditionInfo");

    // ConditionInfo
    const IndexType condition_id = rConditionSettings.GetValue<IndexType>(ConditionSettings::condition_id);
    p_new_cond_info->SetValue(ConditionInfo::condition_id, condition_id);
    const double surface_area = MeshUtilities::AreaOMP(rTriangleMesh);
    p_new_cond_info->SetValue(ConditionInfo::surf_area, surface_area);

    // Add new p_new_cond_info to r_model_info.
    r_model_info.GetList(MainInfo::conditions_infos_list).push_back(std::move(p_new_cond_info));

    /* End: Write to r_model_info */

    // Get again the reference to the just created condition info.
    auto& r_new_condition_info = *(r_model_info.GetList(MainInfo::conditions_infos_list).back());

    // Create new condition and brep_operator.
    Unique<ConditionType> p_new_condition = MakeUnique<ConditionType>(rConditionSettings, r_new_condition_info);
    BRepOperator brep_operator(rTriangleMesh);

    /// Initialize info variables.
    double surf_area_in_active_domain = 0.0;

    // Loop over background grid.
    #pragma omp parallel for reduction(+ : surf_area_in_active_domain) schedule(dynamic)
    for( int for_index = 0; for_index < static_cast<int>(mGridIndexer.NumberOfElements()); ++for_index ) {
		const auto index = static_cast<IndexType>(for_index);
        // Clip embedded geoemetry with the current element (bounding box).
        const auto bounding_box_xyz = mGridIndexer.GetBoundingBoxXYZFromIndex(index);
        auto p_new_mesh = brep_operator.pClipTriangleMeshUnique(bounding_box_xyz.first, bounding_box_xyz.second);

        if( p_new_mesh->NumOfTriangles() > 0 ) {
            const auto p_el = mBackgroundGrid.pGetElement(index+1);
            const double surf_area_segment = MeshUtilities::Area(*p_new_mesh);

            // If p_el != nullptr, the current condition is within an active element.
            // Otherwise, the current condition is outside the active domain.
            auto p_new_segment = p_el ? MakeUnique<ConditionType::ConditionSegmentType>(index, p_el, p_new_mesh)
                : MakeUnique<ConditionType::ConditionSegmentType>(index, p_new_mesh);
            surf_area_in_active_domain += p_el ?  surf_area_segment : 0.0;

            // Add condition segment to condition.
            #pragma omp critical
            p_new_condition->AddSegment(p_new_segment);
        }
    }

    // Add condition to background grid.
    mBackgroundGrid.AddCondition(std::move(p_new_condition));

    /* Begin: Write to ModelInfo */

    // ConditionInfo::
    r_new_condition_info.SetValue(ConditionInfo::perc_surf_area_in_active_domain, surf_area_in_active_domain/surface_area*100.0);

    // ElapsedTimeInfo::
    const double measured_time = timer_condition.Measure();
    auto& r_time_info = r_model_info[MainInfo::elapsed_time_info];
    auto& r_condition_time_info = r_time_info[ElapsedTimeInfo::conditions_time_info];
    const double total_time = (r_time_info.IsSet(ElapsedTimeInfo::total)) ?
        r_time_info.GetValue<double>(ElapsedTimeInfo::total) : 0.0;
    r_time_info.SetValue(ElapsedTimeInfo::total, (total_time+measured_time) );

    // ConditionsTimeInfo::
    const double total_time_conditions = (r_condition_time_info.IsSet(ConditionsTimeInfo::total)) ?
        r_condition_time_info.GetValue<double>(ConditionsTimeInfo::total) : 0.0;
    r_condition_time_info.SetValue(ConditionsTimeInfo::total, (total_time_conditions+measured_time) );

    /* End: Write to ModelInfo */

    // Print condition indo to console.
    PrintConditionInfo(r_new_condition_info);
}

void EmbeddedModel::WriteModelToFile() const {
    const auto& r_settings = GetSettings();
    const auto& r_general_settings = r_settings[MainSettings::general_settings];
    if( r_general_settings.GetValue<bool>(GeneralSettings::write_output_to_file) ) {
        // Start timer.
        Timer timer_output{};

        // Get output_directory_name.
        const std::string output_directory_name = r_general_settings.GetValue<std::string>(GeneralSettings::output_directory_name);
        const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
        QuESo_INFO_IF(echo_level > 0) << ":: WriteFileInfo :: Output directory: '" << output_directory_name << "'\n";

        // Write vtk files (binary = true).
        IO::WriteElementsToVTK(mBackgroundGrid, (output_directory_name + "/elements.vtk"), IO::EncodingType::binary);
        IO::WritePointsToVTK(mBackgroundGrid, (output_directory_name + "/integration_points.vtk"), IO::EncodingType::binary);
        std::for_each(mBackgroundGrid.ConditionsBegin(), mBackgroundGrid.ConditionsEnd(),
            [&output_directory_name](const auto& r_condition){
                IndexType condition_id = r_condition.GetSettings().template GetValue<IndexType>(ConditionSettings::condition_id);
                const std::string bc_filename = output_directory_name + "/condition_id_" + std::to_string(condition_id) + ".stl";
                IO::WriteConditionToSTL(r_condition, bc_filename, IO::EncodingType::binary);
        });

        /* Begin: Write to r_model_info */
        auto& r_model_info = GetModelInfoMutable();

        // ElapsedTimeInfo::
        const double measured_time = timer_output.Measure();
        auto& r_time_info = r_model_info[MainInfo::elapsed_time_info];
        const double total_time = (r_time_info.IsSet(ElapsedTimeInfo::total)) ?
            r_time_info.GetValue<double>(ElapsedTimeInfo::total) : 0.0;
        r_time_info.SetValue(ElapsedTimeInfo::total, (total_time+measured_time) );

        // WriteFilesTimeInfo::
        auto& r_write_files_time_info = r_time_info[ElapsedTimeInfo::write_files_time_info];
        r_write_files_time_info.SetValue(WriteFilesTimeInfo::total, measured_time);

        /* End: Write to r_model_info */

        // Write r_model_info to JSON file.
        IO::WriteDictionaryToJSON(r_model_info, (output_directory_name + "/model_info.json"));

        // Final print to console.
        QuESo_INFO_IF(echo_level > 0) << ":: ElapsedTimeInfo :: Elapsed time: " << measured_time << " sec\n";
    }
}

void EmbeddedModel::CheckIfMeshIsWithinBoundingBox(const TriangleMeshInterface& rTriangleMesh) const {
    const auto& r_settings = GetSettings();
    // Check if bounding box fully contains the triangle mesh.
    if( r_settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level) > 0 ){
        // Get grid dimensions.
        PointType lower_bound = r_settings[MainSettings::background_grid_settings].
            GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
        PointType upper_bound = r_settings[MainSettings::background_grid_settings].
            GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz);

        // Get geometry dimensions.
        auto bb_mesh = MeshUtilities::BoundingBox(rTriangleMesh);

        // Check dimensions
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

void EmbeddedModel::PrintVolumeInfo() const {
    const auto& r_settings = GetSettings();
    const auto& r_model_info = GetModelInfo();
    const auto& r_general_settings = r_settings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    if( echo_level > 0) {
        const auto& r_grid_info = r_model_info[MainInfo::background_grid_info];
        const IndexType num_active_elements = r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_active_elements);
        const IndexType num_trimmed_elements = r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements);
        QuESo_INFO << ":: BackgroundGridInfo :: Number of active elements: " << num_active_elements << std::endl;
        QuESo_INFO << ":: BackgroundGridInfo :: Number of trimmed elements: " << num_trimmed_elements << std::endl;
        const auto& r_quad_info = r_model_info[MainInfo::quadrature_info];
        const IndexType num_quadrature_points = r_quad_info.GetValue<IndexType>(QuadratureInfo::tot_num_points);
        QuESo_INFO << ":: QuadratureRuleInfo :: Number of integration points: " << num_quadrature_points << std::endl;
        if( echo_level > 1 ) {
            const double percentage_of_geometry_volume = r_quad_info.GetValue<double>(QuadratureInfo::percentage_of_geometry_volume);
            QuESo_INFO << ":: QuadratureRuleInfo :: The computed quadrature represents " << percentage_of_geometry_volume
                << "% of the volume of the STL model. Note that this number can depend on the setting 'min_element_volume_ratio'"
                << " in 'QuESoSettings.json'.\n";
        }
    }
}

void EmbeddedModel::PrintVolumeElapsedTimeInfo() const {
    const auto& r_settings = GetSettings();
    const auto& r_model_info = GetModelInfo();
    const auto& r_general_settings = r_settings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    const auto& r_volume_time_info = r_model_info[MainInfo::elapsed_time_info][ElapsedTimeInfo::volume_time_info];
    if( echo_level > 1 ) {
        // Average time spent for each task
        QuESo_INFO << ":: ElpasedTimeInfo (Create Volume) :: Total time: " << r_volume_time_info.GetValue<double>(VolumeTimeInfo::total) << " sec\n";
        QuESo_INFO << ":: ElpasedTimeInfo (Create Volume) :: Individual tasks:\n";
        QuESo_INFO << "   -- Classification of elements:    " <<
            r_volume_time_info.GetValue<double>(VolumeTimeInfo::classification_of_elements) << " sec\n";
        QuESo_INFO << "   -- Computation of intersections:  " <<
            r_volume_time_info.GetValue<double>(VolumeTimeInfo::computation_of_intersections) << " sec\n";
        QuESo_INFO << "   -- Solution of moment fitting eq: " <<
            r_volume_time_info.GetValue<double>(VolumeTimeInfo::solution_of_moment_fitting_eqs) << " sec\n";
        QuESo_INFO << "   -- Construction of GGQ rules:     " <<
            r_volume_time_info.GetValue<double>(VolumeTimeInfo::construction_of_ggq_rules) << " sec\n";
    } else {
        QuESo_INFO_IF(echo_level > 0) << ":: ElpasedTimeInfo (Create Volume) :: Elapsed time: " << r_volume_time_info.GetValue<double>(VolumeTimeInfo::total) << " sec\n";
    }
}

void EmbeddedModel::PrintConditionInfo(const MainDictionaryType& rConditionInfo) const {
    const auto& r_settings = GetSettings();
    const auto& r_general_settings = r_settings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    if( echo_level > 0) {
        IndexType condition_id = rConditionInfo.GetValue<IndexType>(ConditionInfo::condition_id);
        QuESo_INFO << ":: ConditionInfo :: Creating condition with id: " << condition_id << std::endl;
        if( echo_level > 1 ) {
            const auto& surf_area = rConditionInfo.GetValue<double>(ConditionInfo::surf_area);
            QuESo_INFO << "   -- Surface area: " << surf_area << '\n';
            const auto& perc_area = rConditionInfo.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain);
            QuESo_INFO << "   -- " << perc_area
                << "% of the surface area lies within an active element. Note"
                << " that this number can depend on the setting 'min_element_volume_ratio'"
                << " in 'QuESoSettings.json'.\n";
        }
    }
}

void EmbeddedModel::PrintConditionsElapsedTimeInfo() const {
    const auto& r_settings = GetSettings();
    const auto& r_model_info = GetModelInfo();
    const auto& r_general_settings = r_settings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    if( echo_level > 0 ){
        const auto& r_conditions_time_info = r_model_info[MainInfo::elapsed_time_info][ElapsedTimeInfo::conditions_time_info];
        QuESo_INFO << ":: ElpasedTimeInfo (Create Conditions) :: Elapsed time: " << r_conditions_time_info.GetValue<double>(ConditionsTimeInfo::total) << " sec\n";
    }
}

} // End namespace queso
