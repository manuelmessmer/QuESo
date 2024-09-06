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

/// Definition of dictionart keys
enum class Root {main};
enum class Main {general_settings, mesh_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings};
enum class GeneralSettings {input_filename, output_directory_name, echo_level, embedding_flag};
enum class MeshSettings {lower_bound_xyz, upper_bound_xyz, lower_bound_uvw, upper_bound_uvw, polynomial_order, number_of_elements};
enum class TrimmedQuadratureRuleSettings {moment_fitting_residual, min_element_volume_ratio, min_num_boundary_triangles, use_customized_trimmed_points };
enum class NonTrimmedQuadratureRuleSettings {integration_method};

class Settings : public Dictionary<Root, Main, GeneralSettings, MeshSettings, TrimmedQuadratureRuleSettings, NonTrimmedQuadratureRuleSettings>
{
public:
    typedef Dictionary<Root, Main, GeneralSettings, MeshSettings, TrimmedQuadratureRuleSettings, NonTrimmedQuadratureRuleSettings> BaseType;
    typedef std::string Str;
    using BaseType::GetValue;
    using BaseType::SetValue;

    bool Set = true;
    bool DontSet = false;

    Settings() : BaseType(Root::main, Str("Settings")) {
        /// GeneralSettings
        auto& r_general_settings = AddEmptyDictionary(Main::general_settings, Str("GeneralSettings"));
        r_general_settings.AddValues(std::make_tuple(
            std::make_tuple(GeneralSettings::input_filename, Str("input_filename"), Str("dummy"), DontSet ),
            std::make_tuple(GeneralSettings::output_directory_name, Str("output_directory_name"), Str("queso_output"), Set ),
            std::make_tuple(GeneralSettings::echo_level, Str("echo_level"), IndexType(1), Set),
            std::make_tuple(GeneralSettings::embedding_flag, Str("embedding_flag"), true, Set)
        ));

        /// MeshSettings
        auto& r_mesh_settings = AddEmptyDictionary(Main::mesh_settings, Str("mesh_settings"));
        r_mesh_settings.AddValues(std::make_tuple(
            std::make_tuple(MeshSettings::lower_bound_xyz, Str("lower_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet ),
            std::make_tuple(MeshSettings::upper_bound_xyz, Str("upper_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(MeshSettings::lower_bound_uvw, Str("lower_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(MeshSettings::upper_bound_uvw, Str("upper_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(MeshSettings::polynomial_order, Str("polynomial_order"), Vector3i{0, 0, 0}, DontSet ),
            std::make_tuple(MeshSettings::number_of_elements, Str("number_of_elements"), Vector3i{0, 0, 0}, DontSet )
        ));

        /// TrimmedQuadratureRuleSettings
        auto& trimmed_quadrature_rule_settings = AddEmptyDictionary(Main::trimmed_quadrature_rule_settings, Str("trimmed_quadrature_rule_settings"));
        trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(TrimmedQuadratureRuleSettings::moment_fitting_residual, Str("moment_fitting_residual"), 1.0e-10, Set ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_element_volume_ratio, Str("min_element_volume_ratio"), 1.0e-3, Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, Str("min_num_boundary_triangles"), IndexType(100), Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::use_customized_trimmed_points, Str("use_customized_trimmed_points"), false, Set  )
        ));

        /// NonTrimmedQuadratureRuleSettings
        auto& non_trimmed_quadrature_rule_settings = AddEmptyDictionary(Main::non_trimmed_quadrature_rule_settings, Str("non_trimmed_quadrature_rule_settings"));
         non_trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(NonTrimmedQuadratureRuleSettings::integration_method, Str("integration_method"), IntegrationMethod::Gauss, Set )
        ));
    }

private:
    using BaseType::AddValues;
    using BaseType::AddEmptyDictionary;
};

///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP