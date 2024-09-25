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
// IMPORTANT: If key is added here, it must also be added to add_settings_to_python.cpp and to json_import.py.
enum class Root {main_settings};
enum class MainSettings {general_settings, background_grid_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings, conditions_settings_list, testing_settings};
enum class GeneralSettings {input_filename, output_directory_name, echo_level, write_output_to_file};
enum class BackgroundGridSettings {grid_type, lower_bound_xyz, upper_bound_xyz, lower_bound_uvw, upper_bound_uvw, polynomial_order, number_of_elements};
enum class TrimmedQuadratureRuleSettings {moment_fitting_residual, min_element_volume_ratio, min_num_boundary_triangles, neglect_elements_if_stl_is_flawed };
enum class NonTrimmedQuadratureRuleSettings {integration_method};
enum class TestingSettings {use_customized_trimmed_points, embedding_flag};
enum class ConditionSettings {condition_id, condition_type, input_filename, modulus, direction, value, penalty_factor};

typedef Dictionary<Root, MainSettings, GeneralSettings, BackgroundGridSettings, TrimmedQuadratureRuleSettings, NonTrimmedQuadratureRuleSettings, ConditionSettings, IndexType, TestingSettings> SettingsBaseType;

///@name QuESo Classes
///@{

/**
 * @class  Settings.
 * @author Manuel Messmer.
 * @brief  Derived version from Dictionary, which specifies all available Keys, ValueTypes and default values.
 * @see    Dictionary.
**/
class Settings : public SettingsBaseType
{
public:
    ///@}
    ///@name Type definitions
    ///@{

    typedef SettingsBaseType BaseType;
    typedef std::string Str;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor. Sets up the default dictionary.
    /// @param TestingFlag If true, TestingSettings dictionary is added. However, in normal usecase, TestingFlag = false.
    Settings(bool TestingFlag = false) : BaseType(Root::main_settings, Str("settings")) {

        /// Let's define the default dictionary...
        bool Set = true; // Given values are set as default values.
        bool DontSet = false; // Given values are only dummy values used to deduced the associated type.
                            // However, calling IsSet() will return false.

        /// GeneralSettings
        auto& r_general_settings = AddEmptySubDictionary(MainSettings::general_settings, Str("general_settings"));
        r_general_settings.AddValues(std::make_tuple(
            std::make_tuple(GeneralSettings::input_filename, Str("input_filename"), Str("dummy"), DontSet ),
            std::make_tuple(GeneralSettings::output_directory_name, Str("output_directory_name"), Str("queso_output"), Set ),
            std::make_tuple(GeneralSettings::echo_level, Str("echo_level"), IndexType(1), Set),
            std::make_tuple(GeneralSettings::write_output_to_file, Str("write_output_to_file"), true, Set)

        ));

        /// BackgroundGridSettings
        auto& r_background_grid_settings = AddEmptySubDictionary(MainSettings::background_grid_settings, Str("background_grid_settings"));
        r_background_grid_settings.AddValues(std::make_tuple(
            std::make_tuple(BackgroundGridSettings::grid_type, Str("grid_type"), GridType::b_spline_grid, DontSet ),
            std::make_tuple(BackgroundGridSettings::lower_bound_xyz, Str("lower_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::upper_bound_xyz, Str("upper_bound_xyz"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::lower_bound_uvw, Str("lower_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::upper_bound_uvw, Str("upper_bound_uvw"), PointType{0.0, 0.0, 0.0}, DontSet  ),
            std::make_tuple(BackgroundGridSettings::polynomial_order, Str("polynomial_order"), Vector3i{0, 0, 0}, DontSet ),
            std::make_tuple(BackgroundGridSettings::number_of_elements, Str("number_of_elements"), Vector3i{0, 0, 0}, DontSet )
        ));

        /// TrimmedQuadratureRuleSettings
        auto& trimmed_quadrature_rule_settings = AddEmptySubDictionary(MainSettings::trimmed_quadrature_rule_settings, Str("trimmed_quadrature_rule_settings"));
        trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(TrimmedQuadratureRuleSettings::moment_fitting_residual, Str("moment_fitting_residual"), 1.0e-10, Set ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_element_volume_ratio, Str("min_element_volume_ratio"), 1.0e-3, Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, Str("min_num_boundary_triangles"), IndexType(100), Set  ),
            std::make_tuple(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed, Str("neglect_elements_if_stl_is_flawed"), true, Set  )
        ));

        /// NonTrimmedQuadratureRuleSettings
        auto& non_trimmed_quadrature_rule_settings = AddEmptySubDictionary(MainSettings::non_trimmed_quadrature_rule_settings, Str("non_trimmed_quadrature_rule_settings"));
        non_trimmed_quadrature_rule_settings.AddValues(std::make_tuple(
            std::make_tuple(NonTrimmedQuadratureRuleSettings::integration_method, Str("integration_method"), IntegrationMethod::Gauss, Set )
        ));

        /// Optional TestingSettings
        if( TestingFlag ) {
            auto& r_testing_settings = AddEmptySubDictionary(MainSettings::testing_settings, Str("testing_settings"));
            r_testing_settings.AddValues(std::make_tuple(
                std::make_tuple(TestingSettings::use_customized_trimmed_points, Str("use_customized_trimmed_points"), false, Set  ),
                std::make_tuple(TestingSettings::embedding_flag, Str("embedding_flag"), true, Set)
            ));
        }

        /// ConditionSettings
        AddEmptySubDictionary(MainSettings::conditions_settings_list, Str("conditions_settings_list"));
    }

    ///@}
    ///@name Member operations
    ///@{

    /// @brief Creates new condition settings.
    /// @return SettingsBaseType& Reference to dictionary that contains condition settings.
    SettingsBaseType& CreateNewConditionSettings() {
        bool DontSet = false; // Given values are only dummy values used to deduced the associated type.

        auto& r_conditions_settings = (*this)[MainSettings::conditions_settings_list];
        const IndexType size = r_conditions_settings.NumberOfSubDictionaries();
        const std::string condition_name = "condition_" + std::to_string(size);
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
        const double value = (*this)[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio);
        if( value < 0.9e-10 ){
            (*this)[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 1e-10);
        }

        (*this)[MainSettings::general_settings].CheckIfValuesAreSet();
        (*this)[MainSettings::background_grid_settings].CheckIfValuesAreSet();
        (*this)[MainSettings::trimmed_quadrature_rule_settings].CheckIfValuesAreSet();
        (*this)[MainSettings::non_trimmed_quadrature_rule_settings].CheckIfValuesAreSet();
    }

private:

    // Hide the following functions
    using BaseType::AddValues;
    using BaseType::AddEmptySubDictionary;

    ///@}
};

///@}
} // End queso namespace.

#endif // End SETTINGS_INCLUDE_HPP