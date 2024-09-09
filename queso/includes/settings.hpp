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

#ifndef SETTINGS_INCLUDE_HPP
#define SETTINGS_INCLUDE_HPP


/// Project includes
#include "queso/containers/dictionary.hpp"


namespace queso {

/// Definition of dictionary keys
enum class Root {main};
enum class Main {general_settings, mesh_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings, conditions_settings};
enum class GeneralSettings {input_filename, output_directory_name, echo_level, embedding_flag};
enum class MeshSettings {lower_bound_xyz, upper_bound_xyz, lower_bound_uvw, upper_bound_uvw, polynomial_order, number_of_elements};
enum class TrimmedQuadratureRuleSettings {moment_fitting_residual, min_element_volume_ratio, min_num_boundary_triangles, use_customized_trimmed_points };
enum class NonTrimmedQuadratureRuleSettings {integration_method};
enum class ConditionSettings {condition_id, condition_type, input_filename, modulus, direction, value, penalty_factor};

typedef Dictionary<Root, Main, GeneralSettings, MeshSettings, TrimmedQuadratureRuleSettings, NonTrimmedQuadratureRuleSettings, ConditionSettings, IndexType> SettingsBaseType;

class Settings : public SettingsBaseType
{
public:
    typedef SettingsBaseType BaseType;
    typedef std::string Str;
    using BaseType::GetValue;
    using BaseType::SetValue;

    bool Set = true;
    bool DontSet = false;

    Settings() : BaseType(Root::main, Str("settings")) {
        /// GeneralSettings
        auto& r_general_settings = AddEmptySubDictionary(Main::general_settings, Str("general_settings"));
        r_general_settings.AddValues(std::make_tuple(
            std::make_tuple(GeneralSettings::input_filename, Str("input_filename"), Str("dummy"), DontSet ),
            std::make_tuple(GeneralSettings::output_directory_name, Str("output_directory_name"), Str("queso_output"), Set ),
            std::make_tuple(GeneralSettings::echo_level, Str("echo_level"), IndexType(1), Set),
            std::make_tuple(GeneralSettings::embedding_flag, Str("embedding_flag"), true, Set)
        ));

        /// MeshSettings
        auto& r_mesh_settings = AddEmptySubDictionary(Main::mesh_settings, Str("mesh_settings"));
        r_mesh_settings.AddValues(std::make_tuple(
            std::make_tuple(MeshSettings::lower_bound_xyz, Str("lower_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet ),
            std::make_tuple(MeshSettings::upper_bound_xyz, Str("upper_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(MeshSettings::lower_bound_uvw, Str("lower_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(MeshSettings::upper_bound_uvw, Str("upper_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(MeshSettings::polynomial_order, Str("polynomial_order"), Vector3i{0, 0, 0}, DontSet ),
            std::make_tuple(MeshSettings::number_of_elements, Str("number_of_elements"), Vector3i{0, 0, 0}, DontSet )
        ));

        /// TrimmedQuadratureRuleSettings
        auto& trimmed_quadrature_rule_settings = AddEmptySubDictionary(Main::trimmed_quadrature_rule_settings, Str("trimmed_quadrature_rule_settings"));
        trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(TrimmedQuadratureRuleSettings::moment_fitting_residual, Str("moment_fitting_residual"), 1.0e-10, Set ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_element_volume_ratio, Str("min_element_volume_ratio"), 1.0e-3, Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, Str("min_num_boundary_triangles"), IndexType(100), Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::use_customized_trimmed_points, Str("use_customized_trimmed_points"), false, Set  )
        ));

        /// NonTrimmedQuadratureRuleSettings
        auto& non_trimmed_quadrature_rule_settings = AddEmptySubDictionary(Main::non_trimmed_quadrature_rule_settings, Str("non_trimmed_quadrature_rule_settings"));
        non_trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(NonTrimmedQuadratureRuleSettings::integration_method, Str("integration_method"), IntegrationMethod::Gauss, Set )
        ));

        /// ConditionSettings
        AddEmptySubDictionary(Main::conditions_settings, Str("conditions_settings"));
    }

    SettingsBaseType& CreateNewConditionSettings() {
        auto& r_conditions_settings = (*this)[Main::conditions_settings];
        IndexType size = r_conditions_settings.NumberOfSubDictionaries();
        std::string condition_name = "condition_" + std::to_string(size);
        auto& r_new_condition_settings = r_conditions_settings.AddEmptySubDictionary(size, condition_name);
        r_new_condition_settings.AddValues(std::make_tuple(
            std::make_tuple(ConditionSettings::condition_id, Str("condition_id"), IndexType(0), DontSet ),
            std::make_tuple(ConditionSettings::condition_type, Str("condition_type"), std::string("dummy"), DontSet  ),
            std::make_tuple(ConditionSettings::input_filename, Str("input_filename"), std::string("dummy"), DontSet  ),
            std::make_tuple(ConditionSettings::modulus, Str("modulus"), 0.0, DontSet  ),
            std::make_tuple(ConditionSettings::direction, Str("direction"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(ConditionSettings::value, Str("value"), 0.0, DontSet  ),
            std::make_tuple(ConditionSettings::penalty_factor, Str("penalty_factor"), 0.0, DontSet  )
        ));

        return r_new_condition_settings;
    }

private:
    using BaseType::AddValues;
    using BaseType::AddEmptyDictionary;
};




///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP