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

#ifndef REGISTER_KEYS_HPP
#define REGISTER_KEYS_HPP

//// Project includes
#include "queso/includes/keys.hpp"


namespace queso {

// Create value type lists
QuESo_CREATE_VALUE_TYPE_LIST(MainValueTypeList,
    PointType, Vector3i, bool, double, IndexType, std::string, IntegrationMethodType, GridTypeType);

// Create type tags. These can be used with 'static_assert' to check if a given key is of correct type.
// We could create different ValuesTypeTag's here. However, for now, one "Main" seems sufficient.
QuESo_CREATE_KEY_SET_TO_VALUE_TYPE_TAG(MainValuesTypeTag, MainValueTypeList);
QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(ListTypeTag);
QuESo_CREATE_KEY_SET_TO_OBJECT_TYPE_TAG(SubDictTypeTag);

/* --- Keys for Settings --- */

// --- MainSettings ---
// SubDicts
QuESo_DEFINE_KEY_SET( MainSettings, SubDictTypeTag,
    QuESo_KEY_LIST(general_settings, background_grid_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings) );
QuESo_DEFINE_KEY_TO_OBJECT(MainSettings, general_settings, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(MainSettings, background_grid_settings, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(MainSettings, trimmed_quadrature_rule_settings, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(MainSettings, non_trimmed_quadrature_rule_settings, SubDictTypeTag);
QuESo_REGISTER_KEY_SET(MainSettings, SubDictTypeTag,
    QuESo_KEY(MainSettings::general_settings),
    QuESo_KEY(MainSettings::background_grid_settings),
    QuESo_KEY(MainSettings::trimmed_quadrature_rule_settings),
    QuESo_KEY(MainSettings::non_trimmed_quadrature_rule_settings)
);
// Lists
QuESo_DEFINE_KEY_SET( MainSettings, ListTypeTag,
    QuESo_KEY_LIST(conditions_settings_list) );
QuESo_DEFINE_KEY_TO_OBJECT(MainSettings, conditions_settings_list, ListTypeTag);
QuESo_REGISTER_KEY_SET(MainSettings, ListTypeTag,
    QuESo_KEY(MainSettings::conditions_settings_list),
);

// --- GeneralSettings ---
// Values
QuESo_DEFINE_KEY_SET( GeneralSettings, MainValuesTypeTag,
    QuESo_KEY_LIST(input_filename, output_directory_name, echo_level, write_output_to_file) );
QuESo_DEFINE_KEY_TO_VALUE(GeneralSettings, input_filename, MainValuesTypeTag, std::string);
QuESo_DEFINE_KEY_TO_VALUE(GeneralSettings, output_directory_name, MainValuesTypeTag, std::string);
QuESo_DEFINE_KEY_TO_VALUE(GeneralSettings, echo_level, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(GeneralSettings, write_output_to_file, MainValuesTypeTag, bool);
QuESo_REGISTER_KEY_SET(GeneralSettings, MainValuesTypeTag,
    QuESo_KEY(GeneralSettings::input_filename),
    QuESo_KEY(GeneralSettings::output_directory_name),
    QuESo_KEY(GeneralSettings::echo_level),
    QuESo_KEY(GeneralSettings::write_output_to_file)
);

// --- BackgroundGridSettings ---
// Values
QuESo_DEFINE_KEY_SET( BackgroundGridSettings, MainValuesTypeTag,
    QuESo_KEY_LIST(grid_type, lower_bound_xyz, upper_bound_xyz, lower_bound_uvw, upper_bound_uvw, polynomial_order, number_of_elements) );
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, grid_type, MainValuesTypeTag, GridTypeType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, lower_bound_xyz, MainValuesTypeTag, PointType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, upper_bound_xyz, MainValuesTypeTag, PointType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, lower_bound_uvw, MainValuesTypeTag, PointType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, upper_bound_uvw, MainValuesTypeTag, PointType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, polynomial_order, MainValuesTypeTag, Vector3i);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridSettings, number_of_elements, MainValuesTypeTag, Vector3i);
QuESo_REGISTER_KEY_SET(BackgroundGridSettings, MainValuesTypeTag,
    QuESo_KEY(BackgroundGridSettings::grid_type),
    QuESo_KEY(BackgroundGridSettings::lower_bound_xyz),
    QuESo_KEY(BackgroundGridSettings::upper_bound_xyz),
    QuESo_KEY(BackgroundGridSettings::lower_bound_uvw),
    QuESo_KEY(BackgroundGridSettings::upper_bound_uvw),
    QuESo_KEY(BackgroundGridSettings::polynomial_order),
    QuESo_KEY(BackgroundGridSettings::number_of_elements)
);

// --- TrimmedQuadratureRuleSettings ---
// Values
QuESo_DEFINE_KEY_SET( TrimmedQuadratureRuleSettings, MainValuesTypeTag,
    QuESo_KEY_LIST(moment_fitting_residual, min_element_volume_ratio, min_num_boundary_triangles, neglect_elements_if_stl_is_flawed) );
QuESo_DEFINE_KEY_TO_VALUE(TrimmedQuadratureRuleSettings, moment_fitting_residual, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(TrimmedQuadratureRuleSettings, min_element_volume_ratio, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(TrimmedQuadratureRuleSettings, min_num_boundary_triangles, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(TrimmedQuadratureRuleSettings, neglect_elements_if_stl_is_flawed, MainValuesTypeTag, bool);
QuESo_REGISTER_KEY_SET(TrimmedQuadratureRuleSettings, MainValuesTypeTag,
    QuESo_KEY(TrimmedQuadratureRuleSettings::moment_fitting_residual),
    QuESo_KEY(TrimmedQuadratureRuleSettings::min_element_volume_ratio),
    QuESo_KEY(TrimmedQuadratureRuleSettings::min_num_boundary_triangles),
    QuESo_KEY(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed)
);

// --- NonTrimmedQuadratureRuleSettings ---
// Values
QuESo_DEFINE_KEY_SET( NonTrimmedQuadratureRuleSettings, MainValuesTypeTag,
    QuESo_KEY_LIST(integration_method) );
QuESo_DEFINE_KEY_TO_VALUE(NonTrimmedQuadratureRuleSettings, integration_method, MainValuesTypeTag, IntegrationMethodType);
QuESo_REGISTER_KEY_SET(NonTrimmedQuadratureRuleSettings, MainValuesTypeTag,
    QuESo_KEY(NonTrimmedQuadratureRuleSettings::integration_method)
);

// --- ConditionSettings ---
// Values
QuESo_DEFINE_KEY_SET( ConditionSettings, MainValuesTypeTag,
    QuESo_KEY_LIST(condition_id, condition_type, input_filename, modulus, direction, value, penalty_factor) );
QuESo_DEFINE_KEY_TO_VALUE(ConditionSettings, condition_id, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(ConditionSettings, condition_type, MainValuesTypeTag, std::string);
QuESo_DEFINE_KEY_TO_VALUE(ConditionSettings, modulus, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(ConditionSettings, direction, MainValuesTypeTag, PointType);
QuESo_DEFINE_KEY_TO_VALUE(ConditionSettings, value, MainValuesTypeTag, PointType);
QuESo_DEFINE_KEY_TO_VALUE(ConditionSettings, penalty_factor, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(ConditionSettings, MainValuesTypeTag,
    QuESo_KEY(ConditionSettings::condition_id),
    QuESo_KEY(ConditionSettings::condition_type),
    QuESo_KEY(ConditionSettings::modulus),
    QuESo_KEY(ConditionSettings::direction),
    QuESo_KEY(ConditionSettings::value),
    QuESo_KEY(ConditionSettings::penalty_factor)
);

/* --- Keys for ModelInfo --- */

// --- MainInfo ---
// SubDicts
QuESo_DEFINE_KEY_SET( MainInfo, SubDictTypeTag,
    QuESo_KEY_LIST(embedded_geometry_info, quadrature_info, background_grid_info, elapsed_time_info) );
QuESo_DEFINE_KEY_TO_OBJECT(MainInfo, embedded_geometry_info, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(MainInfo, quadrature_info, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(MainInfo, background_grid_info, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(MainInfo, elapsed_time_info, SubDictTypeTag);
QuESo_REGISTER_KEY_SET(MainInfo, SubDictTypeTag,
    QuESo_KEY(MainInfo::embedded_geometry_info),
    QuESo_KEY(MainInfo::quadrature_info),
    QuESo_KEY(MainInfo::background_grid_info),
    QuESo_KEY(MainInfo::elapsed_time_info)
);
// Lists
QuESo_DEFINE_KEY_SET( MainInfo, ListTypeTag,
    QuESo_KEY_LIST(conditions_infos_list) );
QuESo_DEFINE_KEY_TO_OBJECT(MainInfo, conditions_infos_list, ListTypeTag);
QuESo_REGISTER_KEY_SET(MainInfo, ListTypeTag,
    QuESo_KEY(MainInfo::conditions_infos_list)
);

// --- EmbeddedGeometryInfo ---
// Values
QuESo_DEFINE_KEY_SET( EmbeddedGeometryInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(is_closed, volume) );
QuESo_DEFINE_KEY_TO_VALUE(EmbeddedGeometryInfo, is_closed, MainValuesTypeTag, bool);
QuESo_DEFINE_KEY_TO_VALUE(EmbeddedGeometryInfo, volume, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(EmbeddedGeometryInfo, MainValuesTypeTag,
    QuESo_KEY(EmbeddedGeometryInfo::is_closed),
    QuESo_KEY(EmbeddedGeometryInfo::volume)
);

// --- QuadratureInfo ---
// Values
QuESo_DEFINE_KEY_SET( QuadratureInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(represented_volume, percentage_of_geometry_volume, tot_num_points, num_of_points_per_full_element, num_of_points_per_trimmed_element) );
QuESo_DEFINE_KEY_TO_VALUE(QuadratureInfo, represented_volume, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(QuadratureInfo, percentage_of_geometry_volume, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(QuadratureInfo, tot_num_points, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(QuadratureInfo, num_of_points_per_full_element, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(QuadratureInfo, num_of_points_per_trimmed_element, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(QuadratureInfo, MainValuesTypeTag,
    QuESo_KEY(QuadratureInfo::represented_volume),
    QuESo_KEY(QuadratureInfo::percentage_of_geometry_volume),
    QuESo_KEY(QuadratureInfo::tot_num_points),
    QuESo_KEY(QuadratureInfo::num_of_points_per_full_element),
    QuESo_KEY(QuadratureInfo::num_of_points_per_trimmed_element)
);

// --- BackgroundGridInfo ---
// Values
QuESo_DEFINE_KEY_SET( BackgroundGridInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(num_active_elements, num_trimmed_elements, num_full_elements, num_inactive_elements) );
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridInfo, num_active_elements, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridInfo, num_trimmed_elements, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridInfo, num_full_elements, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(BackgroundGridInfo, num_inactive_elements, MainValuesTypeTag, IndexType);
QuESo_REGISTER_KEY_SET(BackgroundGridInfo, MainValuesTypeTag,
    QuESo_KEY(BackgroundGridInfo::num_active_elements),
    QuESo_KEY(BackgroundGridInfo::num_trimmed_elements),
    QuESo_KEY(BackgroundGridInfo::num_full_elements),
    QuESo_KEY(BackgroundGridInfo::num_inactive_elements)
);

// --- ConditionInfo ---
// Values
QuESo_DEFINE_KEY_SET( ConditionInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(condition_id, surf_area, perc_surf_area_in_active_domain) );
QuESo_DEFINE_KEY_TO_VALUE(ConditionInfo, condition_id, MainValuesTypeTag, IndexType);
QuESo_DEFINE_KEY_TO_VALUE(ConditionInfo, surf_area, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(ConditionInfo, perc_surf_area_in_active_domain, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(ConditionInfo, MainValuesTypeTag,
    QuESo_KEY(ConditionInfo::condition_id),
    QuESo_KEY(ConditionInfo::surf_area),
    QuESo_KEY(ConditionInfo::perc_surf_area_in_active_domain)
);

// --- ElapsedTimeInfo ---
// Values
QuESo_DEFINE_KEY_SET( ElapsedTimeInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(total) );
QuESo_DEFINE_KEY_TO_VALUE(ElapsedTimeInfo, total, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(ElapsedTimeInfo, MainValuesTypeTag,
    QuESo_KEY(ElapsedTimeInfo::total)
);
// SubDicts
QuESo_DEFINE_KEY_SET( ElapsedTimeInfo, SubDictTypeTag,
    QuESo_KEY_LIST(volume_time_info, conditions_time_info, write_files_time_info) );
QuESo_DEFINE_KEY_TO_OBJECT(ElapsedTimeInfo, volume_time_info, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(ElapsedTimeInfo, conditions_time_info, SubDictTypeTag);
QuESo_DEFINE_KEY_TO_OBJECT(ElapsedTimeInfo, write_files_time_info, SubDictTypeTag);
QuESo_REGISTER_KEY_SET(ElapsedTimeInfo, SubDictTypeTag,
    QuESo_KEY(ElapsedTimeInfo::volume_time_info),
    QuESo_KEY(ElapsedTimeInfo::conditions_time_info),
    QuESo_KEY(ElapsedTimeInfo::write_files_time_info)
);

// --- VolumeTimeInfo ---
// Values
QuESo_DEFINE_KEY_SET( VolumeTimeInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(total, classification_of_elements, computation_of_intersections, solution_of_moment_fitting_eqs, construction_of_ggq_rules) );
QuESo_DEFINE_KEY_TO_VALUE(VolumeTimeInfo, total, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(VolumeTimeInfo, classification_of_elements, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(VolumeTimeInfo, computation_of_intersections, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(VolumeTimeInfo, solution_of_moment_fitting_eqs, MainValuesTypeTag, double);
QuESo_DEFINE_KEY_TO_VALUE(VolumeTimeInfo, construction_of_ggq_rules, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(VolumeTimeInfo, MainValuesTypeTag,
    QuESo_KEY(VolumeTimeInfo::total),
    QuESo_KEY(VolumeTimeInfo::classification_of_elements),
    QuESo_KEY(VolumeTimeInfo::computation_of_intersections),
    QuESo_KEY(VolumeTimeInfo::solution_of_moment_fitting_eqs),
    QuESo_KEY(VolumeTimeInfo::construction_of_ggq_rules)
);

// --- ConditionsTimeInfo ---
// Values
QuESo_DEFINE_KEY_SET( ConditionsTimeInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(total) );
QuESo_DEFINE_KEY_TO_VALUE(ConditionsTimeInfo, total, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(ConditionsTimeInfo, MainValuesTypeTag,
    QuESo_KEY(ConditionsTimeInfo::total)
);

// --- WriteFilesTimeInfo ---
// Values
QuESo_DEFINE_KEY_SET( WriteFilesTimeInfo, MainValuesTypeTag,
    QuESo_KEY_LIST(total) );
QuESo_DEFINE_KEY_TO_VALUE(WriteFilesTimeInfo, total, MainValuesTypeTag, double);
QuESo_REGISTER_KEY_SET(WriteFilesTimeInfo, MainValuesTypeTag,
    QuESo_KEY(WriteFilesTimeInfo::total)
);

/* --- Keys for Testing --- */
#if QUESO_BUILD_TESTING
namespace Testing {

    // TestKeys1 :: SubDictTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys1, SubDictTypeTag, QuESo_KEY_LIST(zero, one, two, three, four) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, zero, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, one, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, two, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, three, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys1, four, SubDictTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys1, SubDictTypeTag,
        QuESo_KEY(TestKeys1::zero),
        QuESo_KEY(TestKeys1::one),
        QuESo_KEY(TestKeys1::two),
        QuESo_KEY(TestKeys1::three),
        QuESo_KEY(TestKeys1::four)
    );

    // TestKeys2 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys2, ListTypeTag, QuESo_KEY_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, zero, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, one, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, two, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys2, three, ListTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys2, ListTypeTag,
        QuESo_KEY(TestKeys2::zero),
        QuESo_KEY(TestKeys2::one),
        QuESo_KEY(TestKeys2::two),
        QuESo_KEY(TestKeys2::three)
    );

    // TestKeys3 :: MainValuesTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys3, MainValuesTypeTag, QuESo_KEY_LIST(zero, one, two, three, four, five, six, seven) );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, zero, MainValuesTypeTag, PointType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, one, MainValuesTypeTag, Vector3i);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, two, MainValuesTypeTag, bool);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, three, MainValuesTypeTag, double);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, four, MainValuesTypeTag, IndexType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, five, MainValuesTypeTag, std::string);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, six, MainValuesTypeTag, IntegrationMethodType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys3, seven, MainValuesTypeTag, GridTypeType);
    QuESo_REGISTER_KEY_SET( TestKeys3, MainValuesTypeTag,
        QuESo_KEY(TestKeys3::zero),
        QuESo_KEY(TestKeys3::one),
        QuESo_KEY(TestKeys3::two),
        QuESo_KEY(TestKeys3::three),
        QuESo_KEY(TestKeys3::four),
        QuESo_KEY(TestKeys3::five),
        QuESo_KEY(TestKeys3::six),
        QuESo_KEY(TestKeys3::seven)
    );

    // TestKeys4 :: SubDictTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys4, SubDictTypeTag, QuESo_KEY_LIST(zero, one, two, three) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, zero, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, one, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, two, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, three, SubDictTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys4, SubDictTypeTag,
        QuESo_KEY(TestKeys4::zero),
        QuESo_KEY(TestKeys4::one),
        QuESo_KEY(TestKeys4::two),
        QuESo_KEY(TestKeys4::three)
    );
    // TestKeys4 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys4, ListTypeTag, QuESo_KEY_LIST(five, six, seven) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, five, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, six, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys4, seven, ListTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys4, ListTypeTag,
        QuESo_KEY(TestKeys4::five),
        QuESo_KEY(TestKeys4::six),
        QuESo_KEY(TestKeys4::seven)
    );

    // TestKeys5 :: SubDictTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys5, SubDictTypeTag, QuESo_KEY_LIST(zero, one) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, zero, SubDictTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, one, SubDictTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys5, SubDictTypeTag,
        QuESo_KEY(TestKeys5::zero),
        QuESo_KEY(TestKeys5::one)
    );
    // TestKeys5 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys5, ListTypeTag, QuESo_KEY_LIST(five, six) );
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, five, ListTypeTag);
    QuESo_DEFINE_KEY_TO_OBJECT(TestKeys5, six, ListTypeTag);
    QuESo_REGISTER_KEY_SET(TestKeys5, ListTypeTag,
        QuESo_KEY(TestKeys5::five),
        QuESo_KEY(TestKeys5::six)
    );
    // TestKeys5 :: ListTypeTag
    QuESo_DEFINE_KEY_SET( TestKeys5, MainValuesTypeTag, QuESo_KEY_LIST(seven, eight, nine, ten) );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, seven, MainValuesTypeTag, double);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, eight, MainValuesTypeTag, IndexType);
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, nine, MainValuesTypeTag, PointType );
    QuESo_DEFINE_KEY_TO_VALUE(TestKeys5, ten, MainValuesTypeTag, IndexType );
    QuESo_REGISTER_KEY_SET(TestKeys5, MainValuesTypeTag,
        QuESo_KEY(TestKeys5::seven),
        QuESo_KEY(TestKeys5::eight),
        QuESo_KEY(TestKeys5::nine),
        QuESo_KEY(TestKeys5::ten)
    );

}
#endif // End QUESO_BUILD_TESTING

} // End namespace queso

#endif // REGISTER_KEYS_HPP