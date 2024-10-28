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
#include "queso/utilities/math_utilities.hpp"

namespace queso {

/// Definition of dictionary keys
// IMPORTANT: If key is added here, it must also be added to add_settings_to_python.cpp and to json_import.py.

enum class Root {main_settings=DictStarts::start_subdicts};
enum class MainSettings {
    general_settings=DictStarts::start_subdicts, background_grid_settings, trimmed_quadrature_rule_settings, non_trimmed_quadrature_rule_settings,
    conditions_settings_list=DictStarts::start_lists };
enum class GeneralSettings {
    input_filename=DictStarts::start_values, output_directory_name, echo_level, write_output_to_file};
enum class BackgroundGridSettings {
    grid_type=DictStarts::start_values, lower_bound_xyz, upper_bound_xyz, lower_bound_uvw, upper_bound_uvw, polynomial_order, number_of_elements};
enum class TrimmedQuadratureRuleSettings {
    moment_fitting_residual=DictStarts::start_values, min_element_volume_ratio, min_num_boundary_triangles, neglect_elements_if_stl_is_flawed };
enum class NonTrimmedQuadratureRuleSettings {
    integration_method=DictStarts::start_values};
enum class ConditionSettings {
    condition_id=DictStarts::start_values, condition_type, input_filename, modulus, direction, value, penalty_factor};

typedef Dictionary<Root, MainSettings, GeneralSettings, BackgroundGridSettings, TrimmedQuadratureRuleSettings, NonTrimmedQuadratureRuleSettings, ConditionSettings> SettingsBaseType;

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
    Settings() : BaseType(Root::main_settings, Str("settings")) {

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
            std::make_tuple(NonTrimmedQuadratureRuleSettings::integration_method, Str("integration_method"), IntegrationMethod::gauss, Set )
        ));

        /// ConditionSettingsList
        AddEmptyList(MainSettings::conditions_settings_list, Str("conditions_settings_list"));
    }

    ///@}
    ///@name Member operations
    ///@{

    /// @brief Creates new condition settings.
    /// @return SettingsBaseType& Reference to dictionary that contains condition settings.
    SettingsBaseType& CreateNewConditionSettings() {
        bool DontSet = false; // Given values are only dummy values used to deduce the associated type.

        auto& r_conditions_settings = GetListObject(MainSettings::conditions_settings_list);
        auto& r_new_condition_settings = r_conditions_settings.AddListItem(std::make_tuple(
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
    const Settings& Check() const {

        /// Check if all values are set.
        (*this)[MainSettings::general_settings].CheckIfValuesAreSet();
        (*this)[MainSettings::background_grid_settings].CheckIfValuesAreSet();
        (*this)[MainSettings::trimmed_quadrature_rule_settings].CheckIfValuesAreSet();
        (*this)[MainSettings::non_trimmed_quadrature_rule_settings].CheckIfValuesAreSet();

        // condition_id and input_filename must be set for all conditions.
        const auto& r_condition_settings_list = this->GetList(MainSettings::conditions_settings_list);
        IndexType i = 0;
        for( const auto& r_condition_settings : r_condition_settings_list ) {
            QuESo_ERROR_IF( !r_condition_settings.IsSet(ConditionSettings::condition_id) )
                << "'condition_id' of condition (" << i << ") is not set.\n";
            QuESo_ERROR_IF( !r_condition_settings.IsSet(ConditionSettings::input_filename) )
                << "'input_filename' of condition (" << i << ") is not set.\n";
            ++i;
        }

        const IndexType echo_level = (*this)[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level);
        // Orders
        const Vector3i order = (*this)[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);

        const IndexType min_order =  Math::Min( order );
        const IndexType max_order =  Math::Max( order );
        QuESo_ERROR_IF(min_order < 1) << "Invalid Input. The polynomial order must be p > 0. \n";
        QuESo_INFO_IF(max_order > 4) << "Warning :: QuESo is designed to construct efficient quadrature rules for 1 <= p <= 4. "
            << "For higher polynomial degrees, the process might become slow. It is recommended to use quadratic bases. Generally, they offer the best performance and accuracy.\n";

        QuESo_INFO_IF(min_order == 1 && echo_level > 0) << "Info :: When using LINEAR finite elements in combination with rather complex geometries, it can be beneficial to employ quadratic quadrature rules within cut elements. "
            << "Linear quadrature rules will converge to the correct solution when using fine discretizations. However, quadratic quadrature rules "
            << "are simply more suited to capture complex cut domains and can hence provide better results for coarse meshes. "
            << "Thus, if you are integrating LINEAR finite elements, consider using '\"polynomial_order\" : [2, 2, 2]' and '\"integration_method\" : \"Gauss_Reduced1\"'. This will generate quadratic quadrature rules in all cut elements "
            << "and linear Gauss rules for all full/interior elements.\n";

        // Number of elements
        const Vector3i num_elements = (*this)[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements);
        const IndexType tot_num_elements = num_elements[0]*num_elements[1]*num_elements[2];
        QuESo_INFO_IF( tot_num_elements < 2 && echo_level > 0 ) << "You are using only one single element.\n";

        const IntegrationMethod integration_method = (*this)[MainSettings::non_trimmed_quadrature_rule_settings]
            .GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method);

        QuESo_ERROR_IF( min_order < 2 && integration_method == IntegrationMethod::gauss_reduced_2)
            << "Gauss_Reduced2 is only applicable to background grids with at least p=2.\n";

        // Check if ggq_rules are feasible.
        const bool ggq_rule_ise_used =  static_cast<int>(integration_method) >= 3;
        QuESo_ERROR_IF(ggq_rule_ise_used && max_order > 4) << "Generalized Gauss Quadrature (GGQ) rules are only available for p <= 4.\n";

        const GridType grid_type = (*this)[MainSettings::background_grid_settings]
            .GetValue<GridType>(BackgroundGridSettings::grid_type);
        QuESo_ERROR_IF( ggq_rule_ise_used && (grid_type != GridType::b_spline_grid)) << "GGQ_Rules can only be used in combination with 'grid_type' : 'b_spline_grid'.\n";

        QuESo_ERROR_IF(ggq_rule_ise_used && min_order < 2) << "Generalized Gauss Quadrature (GGQ) rules are only applicable to B-Spline meshes with at least p=2.\n";

        return *this;
    }

private:

    /// Hide the following functions
    /// @note Note making these functions protected in the Base class does not help, sine they are called here on
    /// a reference to base class instance.
    using BaseType::AddValues;
    using BaseType::AddListItem;
    using BaseType::GetListObject;
    using BaseType::AddEmptyList;
    using BaseType::AddEmptySubDictionary;

    ///@}
};

///@}
} // End queso namespace.

#endif // End SETTINGS_INCLUDE_HPP