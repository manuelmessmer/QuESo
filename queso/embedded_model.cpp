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
#include "queso/embedded_model.h"
#include "queso/io/io_utilities.h"
#include "queso/utilities/mesh_utilities.h"
#include "queso/embedding/brep_operator.h"
#include "queso/quadrature/single_element.hpp"
#include "queso/quadrature/trimmed_element.hpp"
#include "queso/quadrature/multiple_elements.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {


void EmbeddedModel::ComputeVolume(const TriangleMeshInterface& rTriangleMesh){

    CheckIfMeshIsWithinBoundingBox(rTriangleMesh);

    /// Set ModelInfo
    // EmbeddedGeometryInfo
    const double volume = MeshUtilities::VolumeOMP(rTriangleMesh);
    mModelInfo[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::volume, volume);
    const bool is_closed = MeshUtilities::EstimateQuality(rTriangleMesh) < 1e-10;
    mModelInfo[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::is_closed, is_closed);

    // Get neccessary settings
    const IntegrationMethod integration_method = mSettings[MainSettings::non_trimmed_quadrature_rule_settings]
        .GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method);
    const bool ggq_rule_ise_used =  static_cast<int>(integration_method) >= 3;
    const auto& r_trimmed_quad_rule_settings = mSettings[MainSettings::trimmed_quadrature_rule_settings];
    const double min_vol_element_ratio = std::max<double>(r_trimmed_quad_rule_settings.GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-10);
    const IndexType num_boundary_triangles = r_trimmed_quad_rule_settings.GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles);
    const double moment_fitting_residual = r_trimmed_quad_rule_settings.GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
    const bool neglect_elements_if_stl_is_flawed = r_trimmed_quad_rule_settings.GetValue<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed);
    const auto& r_grid_settings = mSettings[MainSettings::background_grid_settings];
    const Vector3i polynomial_order = r_grid_settings.GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);
    const Vector3i number_of_elements = r_grid_settings.GetValue<Vector3i>(BackgroundGridSettings::number_of_elements);
    const IndexType echo_level = mSettings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level);

    // Start timer
    Timer timer_total{};

    // Reserve element container
    const IndexType global_number_of_elements = mGridIndexer.NumberOfElements();
    mBackgroundGrid.ReserveElements(global_number_of_elements);

    // Construct BRepOperator
    BRepOperator brep_operator(rTriangleMesh);

    // Classify all elements.
    Timer timer_check_intersect{};
    Unique<BRepOperator::StatusVectorType> p_classifications = brep_operator.pGetElementClassifications(mSettings[MainSettings::background_grid_settings]);
    auto& r_volume_time_info = mModelInfo[MainInfo::elapsed_time_info][ElapsedTimeInfo::volume_time_info];
    r_volume_time_info.SetValue(VolumeTimeInfo::classification_of_elements, timer_check_intersect.Measure());

    //// Info variables
    // TimeInfo
    double et_compute_intersection = 0.0;
    double et_moment_fitting = 0.0;
    // GridInfo
    SizeType num_active_elements = 0;
    SizeType num_trimmed_elements = 0;
    // Num of threads
    IndexType num_threads = 1;

    // Loop over all elements
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();

        #pragma omp for reduction(+ : et_compute_intersection, et_moment_fitting, num_active_elements, num_trimmed_elements) schedule(dynamic)
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
                    auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(bounding_box_xyz.first, bounding_box_xyz.second,
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
                    ++num_active_elements;
                    if( new_element->IsTrimmed() ) { ++num_trimmed_elements; }
                    #pragma omp critical // TODO: improve this.
                    mBackgroundGrid.AddElement(new_element); // After this new_element is a null_ptr. Is std::moved to container.
                }
            }
        } /// #pragma omp for reduction
    } /// End #pragma omp parallel

    /// Assmble Generalized Gaussian quadrature rules (if enabled).
    double et_ggq_rules = 0.0;
    if( ggq_rule_ise_used ){
        Timer timer_ggq_rules{};
        QuadratureMultipleElements<ElementType>::AssembleIPs(mBackgroundGrid, number_of_elements, polynomial_order, integration_method);
        et_ggq_rules = timer_ggq_rules.Measure();
    }
    const double elapsed_time_total = timer_total.Measure();

    /// Set ModelInfo
    // ElpasedTimeInfo
    auto& r_elapsed_time_info = mModelInfo[MainInfo::elapsed_time_info];
    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::total, elapsed_time_total);
    const double total_time = (r_elapsed_time_info.IsSet(ElapsedTimeInfo::total)) ?
        r_elapsed_time_info.GetValue<double>(ElapsedTimeInfo::total) : 0.0;
    r_elapsed_time_info.SetValue(ElapsedTimeInfo::total, (total_time+elapsed_time_total) );
    r_volume_time_info.SetValue(VolumeTimeInfo::computation_of_intersections, et_compute_intersection / ((double) num_threads) );
    r_volume_time_info.SetValue(VolumeTimeInfo::solution_of_moment_fitting_eqs, et_moment_fitting / ((double) num_threads) );
    r_volume_time_info.SetValue(VolumeTimeInfo::construction_of_ggq_rules, et_ggq_rules);

    // BackgroundGridInfo
    mModelInfo[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_active_elements, num_active_elements );
    mModelInfo[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_trimmed_elements, num_trimmed_elements );
    const IndexType num_full_elements = (num_active_elements-num_trimmed_elements);
    mModelInfo[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_full_elements, num_full_elements );
    const IndexType num_inactive_elements =  ( mGridIndexer.NumberOfElements() - num_active_elements);
    mModelInfo[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_inactive_elements, num_inactive_elements );

    // QuadratureInfo
    double represented_volume = 0.0;
    SizeType tot_num_points_full = 0;
    SizeType tot_num_points_trimmed = 0;
    const auto el_it_ptr_begin = mBackgroundGrid.ElementsBegin();
    #pragma omp parallel for reduction(+ : represented_volume) reduction(+ : tot_num_points_full) reduction(+ : tot_num_points_trimmed)
    for( int i = 0; i < static_cast<int>(num_active_elements); ++i ){
        const auto& el_ptr = *(el_it_ptr_begin + i);
        const double det_j = el_ptr->DetJ();
        const auto& r_points = el_ptr->GetIntegrationPoints();
        for( const auto& r_point : r_points ){
            represented_volume += r_point.Weight()*det_j;
        }
        if( el_ptr->IsTrimmed() ){
            tot_num_points_trimmed += r_points.size();
        } else {
            tot_num_points_full += r_points.size();
        }
    }
    const SizeType tot_num_points = tot_num_points_trimmed + tot_num_points_full;
    mModelInfo[MainInfo::quadrature_info].SetValue(QuadratureInfo::tot_num_points, tot_num_points);
    mModelInfo[MainInfo::quadrature_info].SetValue(QuadratureInfo::represented_volume, represented_volume);
    mModelInfo[MainInfo::quadrature_info].SetValue(QuadratureInfo::percentage_of_geometry_volume, represented_volume/volume*100.0);
    mModelInfo[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_full_element, ((double)tot_num_points_full)/num_full_elements);
    mModelInfo[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_trimmed_element, ((double)tot_num_points_trimmed)/num_trimmed_elements);

    PrintVolumeInfo();
}

void EmbeddedModel::ComputeCondition(const TriangleMeshInterface& rTriangleMesh, const SettingsBaseType& rConditionSettings) {

    CheckIfMeshIsWithinBoundingBox(rTriangleMesh);

    Timer timer_condition{};
    auto& r_new_cond_info = mModelInfo.CreateNewConditionInfo();
    /// Set ModelInfo
    // ConditionInfo
    const IndexType condition_id = rConditionSettings.GetValue<IndexType>(ConditionSettings::condition_id);
    r_new_cond_info.SetValue(ConditionInfo::condition_id, condition_id);
    const double surface_area = MeshUtilities::AreaOMP(rTriangleMesh);
    r_new_cond_info.SetValue(ConditionInfo::surf_area, surface_area);

    // Create new condition and brep_operator
    Unique<ConditionType> p_new_condition = MakeUnique<ConditionType>(rConditionSettings, r_new_cond_info);
    BRepOperator brep_operator(rTriangleMesh);

    /// Info variables
    double surf_area_segments = 0.0;
    double surf_area_in_active_domain = 0.0;

    // Loop over background grid
    #pragma omp parallel for reduction(+ : surf_area_segments, surf_area_in_active_domain)
    for( int index = 0; index < static_cast<int>(mGridIndexer.NumberOfElements()); ++index ) {
        const auto bounding_box_xyz = mGridIndexer.GetBoundingBoxXYZFromIndex(index);
        auto p_new_mesh = brep_operator.pClipTriangleMeshUnique(bounding_box_xyz.first, bounding_box_xyz.second);
        if( p_new_mesh->NumOfTriangles() > 0 ) {
            const auto p_el = mBackgroundGrid.pGetElement(index+1);
            const double surf_area_segment = MeshUtilities::Area(*p_new_mesh);
            surf_area_segments += surf_area_segment;
            auto p_new_segment = p_el ? MakeUnique<ConditionType::ConditionSegmentType>(index, p_el, p_new_mesh)
                : MakeUnique<ConditionType::ConditionSegmentType>(index, p_new_mesh);
            surf_area_in_active_domain += p_el ?  surf_area_segment : 0.0;
            #pragma omp critical
            p_new_condition->AddSegment(p_new_segment);
        }
    }
    mBackgroundGrid.AddCondition(p_new_condition);

    /// Set ModelInfo
    // ConditionInfo
    r_new_cond_info.SetValue(ConditionInfo::perc_surf_area_in_active_domain, surf_area_in_active_domain/surface_area*100.0);
    auto& r_condition_time_info = mModelInfo[MainInfo::elapsed_time_info][ElapsedTimeInfo::conditions_time_info];
    // ElapsedTimeInfo
    const double total_time = (r_condition_time_info.IsSet(ConditionsTimeInfo::total)) ?
        r_condition_time_info.GetValue<double>(ConditionsTimeInfo::total) : 0.0;
    r_condition_time_info.SetValue(ConditionsTimeInfo::total, (total_time+timer_condition.Measure()) );

    PrintConditionInfo(r_new_cond_info);
}

void EmbeddedModel::CheckIfMeshIsWithinBoundingBox(const TriangleMeshInterface& rTriangleMesh) const {
    // Check if bounding box fully contains the triangle mesh.
    if( mSettings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level) > 0 ){
        PointType lower_bound = mSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
        PointType upper_bound = mSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz);
        auto bb_mesh = MeshUtilities::BoundingBox(rTriangleMesh);
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
    const auto& r_general_settings = mSettings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    if( echo_level > 0) {
        const auto& r_grid_info = mModelInfo[MainInfo::background_grid_info];
        const IndexType num_active_elements = r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_active_elements);
        const IndexType num_trimmed_elements = r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements);
        QuESo_INFO << ":: BackgroundGridInfo :: Number of active elements: " << num_active_elements << std::endl;
        QuESo_INFO << ":: BackgroundGridInfo :: Number of trimmed elements: " << num_trimmed_elements << std::endl;
        const auto& r_quad_info = mModelInfo[MainInfo::quadrature_info];
        const IndexType num_quadrature_points = r_quad_info.GetValue<IndexType>(QuadratureInfo::tot_num_points);
        QuESo_INFO << ":: QuadratureRuleInfo :: Number of integration points: " << num_quadrature_points << std::endl;
        if( echo_level > 1 ) {
            const auto& r_quad_info = mModelInfo[MainInfo::quadrature_info];
            const double percentage_of_geometry_volume = r_quad_info.GetValue<double>(QuadratureInfo::percentage_of_geometry_volume);
            QuESo_INFO << ":: QuadratureRuleInfo :: The computed quadrature represents " << percentage_of_geometry_volume
                << "% of the volume of the STL model. Note that this number can depend on the setting 'min_element_volume_ratio'"
                << " in 'QuESoSettings.json'.\n";
        }
    }
}

void EmbeddedModel::PrintVolumeElapsedTimeInfo() const {
    const auto& r_general_settings = mSettings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    const auto& r_volume_time_info = mModelInfo[MainInfo::elapsed_time_info][ElapsedTimeInfo::volume_time_info];
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

void EmbeddedModel::PrintConditionInfo(const ModelInfoBaseType& rConditionInfo) const {
    const auto& r_general_settings = mSettings[MainSettings::general_settings];
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
    const auto& r_general_settings = mSettings[MainSettings::general_settings];
    const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
    if( echo_level > 0 ){
        const auto& r_conditions_time_info = mModelInfo[MainInfo::elapsed_time_info][ElapsedTimeInfo::conditions_time_info];
        QuESo_INFO << ":: ElpasedTimeInfo (Create Conditions) :: Elapsed time: " << r_conditions_time_info.GetValue<double>(ConditionsTimeInfo::total) << " sec\n";
    }
}

} // End namespace queso
