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

#ifndef CHECK_DICTIONARY_UTILITIES_INCLUDE_HPP
#define CHECK_DICTIONARY_UTILITIES_INCLUDE_HPP

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/dictionary.hpp"
#include "queso/utilities/math_utilities.hpp"

namespace queso {
namespace CheckDictionaryUtilities {
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    /// @brief Throws if any value in the data set on the current level
    ///        of the given dictionary is not set. Does not check sub dictionaries.
    /// @tparam TKeyInfoType
    /// @param rSettings
    template<typename TKeyInfoType>
    void CheckIfAllValuesInDataSetAreSet( const MainDictionaryType& rSettings ) {
        using StringAccess = DictionaryStringAccess<MainDictionaryType>;
        for( const auto& [name_view, _] : TKeyInfoType::msStringToKeyMap ) {
            std::string name{name_view};
            QuESo_ERROR_IF( !StringAccess::IsSet(rSettings, name) ) << "Setting '" << name << "' is not set.\n";
        }
    }

    /// @brief Throws if the input settings are not set properly.
    /// @param rSettings
    void CheckSettings( const MainDictionaryType& rSettings ) {
        const IndexType echo_level = rSettings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level);

        // Check if all values are set.
        CheckIfAllValuesInDataSetAreSet<key::detail::GeneralSettingsMainValuesTypeTagKeySetInfo>(
            rSettings[MainSettings::general_settings]);
        CheckIfAllValuesInDataSetAreSet<key::detail::BackgroundGridSettingsMainValuesTypeTagKeySetInfo>(
            rSettings[MainSettings::background_grid_settings]);
        CheckIfAllValuesInDataSetAreSet<key::detail::TrimmedQuadratureRuleSettingsMainValuesTypeTagKeySetInfo>(
            rSettings[MainSettings::trimmed_quadrature_rule_settings]);
        CheckIfAllValuesInDataSetAreSet<key::detail::NonTrimmedQuadratureRuleSettingsMainValuesTypeTagKeySetInfo>(
            rSettings[MainSettings::non_trimmed_quadrature_rule_settings]);

        // Check if 'condition_id' and 'input_filename' are set for all conditions.
        const auto& r_condition_settings_list = rSettings.GetList(MainSettings::conditions_settings_list);
        IndexType i = 0;
        for( const auto& p_condition_settings : r_condition_settings_list ) {
            QuESo_ERROR_IF( !p_condition_settings->IsSet(ConditionSettings::condition_id) )
                << "'condition_id' of condition (" << i << ") is not set.\n";
            QuESo_ERROR_IF( !p_condition_settings->IsSet(ConditionSettings::input_filename) )
                << "'input_filename' of condition (" << i << ") is not set.\n";
            ++i;
        }

        // Check 'polynomial_order' related settings.
        const Vector3i order = rSettings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);

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

        // Check 'number_of_elements' related settings.
        const Vector3i num_elements = rSettings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements);
        const IndexType tot_num_elements = num_elements[0]*num_elements[1]*num_elements[2];
        QuESo_INFO_IF( tot_num_elements < 2 && echo_level > 0 ) << "You are using only one single element.\n";

        const IntegrationMethod integration_method = rSettings[MainSettings::non_trimmed_quadrature_rule_settings]
            .GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method);

        QuESo_ERROR_IF( min_order < 2 && integration_method == IntegrationMethod::gauss_reduced_2)
            << "'Gauss_Reduced2' is only applicable to background grids with at least p=2.\n";

        // Check if ggq_rules are feasible.
        const bool ggq_rule_ise_used =  static_cast<int>(integration_method) >= 3;
        QuESo_ERROR_IF(ggq_rule_ise_used && max_order > 4) << "Generalized Gauss Quadrature (GGQ) rules are only available for p <= 4.\n";

        const GridType grid_type = rSettings[MainSettings::background_grid_settings]
            .GetValue<GridType>(BackgroundGridSettings::grid_type);
        QuESo_ERROR_IF( ggq_rule_ise_used && (grid_type != GridType::b_spline_grid)) << "GGQ_Rules can only be used in combination with 'grid_type' : 'b_spline_grid'.\n";

        QuESo_ERROR_IF(ggq_rule_ise_used && min_order < 2) << "Generalized Gauss Quadrature (GGQ) rules are only applicable to B-Spline meshes with at least p=2.\n";
    }

} // End namespace CheckSettingsUtilities
} // End namespace queso

#endif // CHECK_DICTIONARY_UTILITIES_INCLUDE_HPP