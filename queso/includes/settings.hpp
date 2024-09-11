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
enum class Main {general_settings, background_grid_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings, conditions_settings_list, testing_settings};
enum class GeneralSettings {input_filename, output_directory_name, echo_level, write_output_to_file};
enum class BackgroundGridSettings {grid_type, lower_bound_xyz, upper_bound_xyz, lower_bound_uvw, upper_bound_uvw, polynomial_order, number_of_elements};
enum class TrimmedQuadratureRuleSettings {moment_fitting_residual, min_element_volume_ratio, min_num_boundary_triangles };
enum class NonTrimmedQuadratureRuleSettings {integration_method};
enum class TestingSettings {use_customized_trimmed_points, embedding_flag};
enum class ConditionSettings {condition_id, condition_type, input_filename, modulus, direction, value, penalty_factor};

typedef Dictionary<Root, Main, GeneralSettings, BackgroundGridSettings, TrimmedQuadratureRuleSettings, NonTrimmedQuadratureRuleSettings, ConditionSettings, IndexType, TestingSettings> SettingsBaseType;

class Settings : public SettingsBaseType
{
public:
    typedef SettingsBaseType BaseType;
    typedef std::string Str;

    bool Set = true;
    bool DontSet = false;

    Settings(bool TestingFlag = false) : BaseType(Root::main, Str("settings")) {

        /// GeneralSettings
        auto& r_general_settings = AddEmptySubDictionary(Main::general_settings, Str("general_settings"));
        r_general_settings.AddValues(std::make_tuple(
            std::make_tuple(GeneralSettings::input_filename, Str("input_filename"), Str("dummy"), DontSet ),
            std::make_tuple(GeneralSettings::output_directory_name, Str("output_directory_name"), Str("queso_output"), Set ),
            std::make_tuple(GeneralSettings::echo_level, Str("echo_level"), IndexType(1), Set),
            std::make_tuple(GeneralSettings::write_output_to_file, Str("write_output_to_file"), true, Set)

        ));

        /// BackgroundGridSettings
        auto& r_background_grid_settings = AddEmptySubDictionary(Main::background_grid_settings, Str("background_grid_settings"));
        r_background_grid_settings.AddValues(std::make_tuple(
            std::make_tuple(BackgroundGridSettings::grid_type, Str("grid_type"), BackgroundGridType::b_spline_grid, DontSet ),
            std::make_tuple(BackgroundGridSettings::lower_bound_xyz, Str("lower_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::upper_bound_xyz, Str("upper_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::lower_bound_uvw, Str("lower_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::upper_bound_uvw, Str("upper_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::polynomial_order, Str("polynomial_order"), Vector3i{0, 0, 0}, DontSet ),
            std::make_tuple(BackgroundGridSettings::number_of_elements, Str("number_of_elements"), Vector3i{0, 0, 0}, DontSet )
        ));

        /// TrimmedQuadratureRuleSettings
        auto& trimmed_quadrature_rule_settings = AddEmptySubDictionary(Main::trimmed_quadrature_rule_settings, Str("trimmed_quadrature_rule_settings"));
        trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(TrimmedQuadratureRuleSettings::moment_fitting_residual, Str("moment_fitting_residual"), 1.0e-10, Set ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_element_volume_ratio, Str("min_element_volume_ratio"), 1.0e-3, Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, Str("min_num_boundary_triangles"), IndexType(100), Set  )
        ));

        /// NonTrimmedQuadratureRuleSettings
        auto& non_trimmed_quadrature_rule_settings = AddEmptySubDictionary(Main::non_trimmed_quadrature_rule_settings, Str("non_trimmed_quadrature_rule_settings"));
        non_trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(NonTrimmedQuadratureRuleSettings::integration_method, Str("integration_method"), IntegrationMethod::Gauss, Set )
        ));

        if( TestingFlag ) {
            auto& r_testing_settings = AddEmptySubDictionary(Main::testing_settings, Str("testing_settings"));
            r_testing_settings.AddValues(std::make_tuple(
                std::make_tuple(TestingSettings::use_customized_trimmed_points, Str("use_customized_trimmed_points"), false, Set  ),
                std::make_tuple(TestingSettings::embedding_flag, Str("embedding_flag"), true, Set)
            ));
        }

        /// ConditionSettings
        AddEmptySubDictionary(Main::conditions_settings_list, Str("conditions_settings_list"));
    }

    SettingsBaseType& CreateNewConditionSettings() {
        auto& r_conditions_settings = (*this)[Main::conditions_settings_list];
        IndexType size = r_conditions_settings.NumberOfSubDictionaries();
        std::string condition_name = "condition_" + std::to_string(size);
        auto& r_new_condition_settings = r_conditions_settings.AddEmptySubDictionary(size, condition_name);
        r_new_condition_settings.AddValues(std::make_tuple(
            std::make_tuple(ConditionSettings::condition_id, Str("condition_id"), IndexType(0), DontSet ),
            std::make_tuple(ConditionSettings::condition_type, Str("condition_type"), std::string("dummy"), DontSet  ),
            std::make_tuple(ConditionSettings::input_filename, Str("input_filename"), std::string("dummy"), DontSet  ),
            std::make_tuple(ConditionSettings::modulus, Str("modulus"), 0.0, DontSet  ),
            std::make_tuple(ConditionSettings::direction, Str("direction"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(ConditionSettings::value, Str("value"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(ConditionSettings::penalty_factor, Str("penalty_factor"), 0.0, DontSet  )
        ));

        return r_new_condition_settings;
    }

    /// @brief Check values.
    void Check() {
        // Make sure this value is not numerically zero.
        const double value = (*this)[Main::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio);
        if( value < 0.9e-10 ){
            (*this)[Main::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 1e-10);
        }

        (*this)[Main::general_settings].CheckIfValuesAreSet();
        (*this)[Main::background_grid_settings].CheckIfValuesAreSet();
        (*this)[Main::trimmed_quadrature_rule_settings].CheckIfValuesAreSet();
        (*this)[Main::non_trimmed_quadrature_rule_settings].CheckIfValuesAreSet();
    }

private:
    // Hide the following functions
    using BaseType::AddValues;
    using BaseType::AddEmptySubDictionary;
};




///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP