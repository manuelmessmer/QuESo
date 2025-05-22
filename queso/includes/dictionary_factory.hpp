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

#ifndef DEFINE_DICTIONARY_FACTORY_HPP
#define DEFINE_DICTIONARY_FACTORY_HPP

/// STL includes
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
/// Project includes
#include "queso/containers/dictionary.hpp"

namespace queso {

using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

namespace factories {

/// @brief Returns a Dictionary containing the default (main) settings.
/// @return Unique<MainDictionaryType>
inline Unique<MainDictionaryType> CreateSettings() {
    using DictionaryType = MainDictionaryType;

    /*--- Main ---*/
    auto p_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        DictionaryType::EmptyKeySetType,
        key::detail::MainSettingsSubDictTypeTagKeySetInfo,
        key::detail::MainSettingsListTypeTagKeySetInfo>{} );

    /*--- General settings ---*/
    auto p_general_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::GeneralSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Add defaults.
    p_general_settings->SetValue(GeneralSettings::output_directory_name, std::string("queso_output"));
    p_general_settings->SetValue(GeneralSettings::echo_level, 1u);
    p_general_settings->SetValue(GeneralSettings::write_output_to_file, true);

    // Add to Main.
    p_settings->SetSubDictionary(MainSettings::general_settings, std::move(p_general_settings));

    /*--- Background grid settings ---*/
    auto p_background_grid_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::BackgroundGridSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Add to Main.
    p_settings->SetSubDictionary(MainSettings::background_grid_settings, std::move(p_background_grid_settings));

    /*--- Background grid settings ---*/
    auto p_trimmed_quadrature_rule_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::TrimmedQuadratureRuleSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Add defaults.
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 1.0e-10);
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 1.0e-3);
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, 100u);
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed, true);

    // Add to Main.
    p_settings->SetSubDictionary(MainSettings::trimmed_quadrature_rule_settings, std::move(p_trimmed_quadrature_rule_settings));

    /*--- Non background grid settings ---*/
    auto p_non_trimmed_quadrature_rule_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::NonTrimmedQuadratureRuleSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Set Defaults
    p_non_trimmed_quadrature_rule_settings->SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::gauss);

    // Add to Main.
    p_settings->SetSubDictionary(MainSettings::non_trimmed_quadrature_rule_settings, std::move(p_non_trimmed_quadrature_rule_settings));

    return p_settings;
}

/// @brief Returns a Dictionary containing the default condition settings.
/// @return Unique<MainDictionaryType>
inline Unique<MainDictionaryType> CreateConditionSettings() {
    using DictionaryType = MainDictionaryType;

    /*--- Condition settings ---*/
    auto p_condition_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ConditionSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    return p_condition_settings;
}

/// @brief Returns a Dictionary containing the default model info.
/// @return Unique<MainDictionaryType>
inline Unique<MainDictionaryType> CreateModelInfo() {
    using DictionaryType = MainDictionaryType;

    /*--- Main info ---*/
    auto p_model_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        DictionaryType::EmptyKeySetType,
        key::detail::MainInfoSubDictTypeTagKeySetInfo,
        key::detail::MainInfoListTypeTagKeySetInfo>{} );

    /*--- Embedded geometry info ---*/
    auto p_embedded_geometry_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::EmbeddedGeometryInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Add to main.
    p_model_info->SetSubDictionary(MainInfo::embedded_geometry_info, std::move(p_embedded_geometry_info));

    /*--- Quadrature info ---*/
    auto p_quadrature_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::QuadratureInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Add to main.
    p_model_info->SetSubDictionary(MainInfo::quadrature_info, std::move(p_quadrature_info));

    /*--- Background grid info ---*/
    auto p_background_grid_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::BackgroundGridInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Add to main.
    p_model_info->SetSubDictionary(MainInfo::background_grid_info, std::move(p_background_grid_info));

    /*--- Elapsed time info ---*/
    auto p_elapsed_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ElapsedTimeInfoMainValuesTypeTagKeySetInfo,
        key::detail::ElapsedTimeInfoSubDictTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType>{} );

    // Set defaults.
    p_elapsed_time_info->SetValue(ElapsedTimeInfo::total, 0.0);

    /*--- Elapsed volume time info ---*/
    auto p_volume_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::VolumeTimeInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Set defaults.
    p_volume_time_info->SetValue(VolumeTimeInfo::total, 0.0);
    p_volume_time_info->SetValue(VolumeTimeInfo::classification_of_elements, 0.0);
    p_volume_time_info->SetValue(VolumeTimeInfo::computation_of_intersections, 0.0);
    p_volume_time_info->SetValue(VolumeTimeInfo::construction_of_ggq_rules, 0.0);
    p_volume_time_info->SetValue(VolumeTimeInfo::solution_of_moment_fitting_eqs, 0.0);

    // Add volume_time_info to elapsed_time_info.
    p_elapsed_time_info->SetSubDictionary(ElapsedTimeInfo::volume_time_info, std::move(p_volume_time_info));

    /*--- Conditions time info ---*/
    auto p_conditions_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ConditionsTimeInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Set defaults.
    p_conditions_time_info->SetValue(ConditionsTimeInfo::total, 0.0);

    // Add conditions_time_info to elapsed_time_info.
    p_elapsed_time_info->SetSubDictionary(ElapsedTimeInfo::conditions_time_info, std::move(p_conditions_time_info));

    /*--- Write files time info ---*/
    auto p_write_files_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::WriteFilesTimeInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    // Set defaults.
    p_write_files_time_info->SetValue(WriteFilesTimeInfo::total, 0.0);

    // Add write_files_time_info to elapsed_time_info.
    p_elapsed_time_info->SetSubDictionary(ElapsedTimeInfo::write_files_time_info, std::move(p_write_files_time_info));

    // Add elapsed_time_info to Main.
    p_model_info->SetSubDictionary(MainInfo::elapsed_time_info, std::move(p_elapsed_time_info));

    return p_model_info;
}

/// @brief Returns a Dictionary containing the default condition info.
/// @return Unique<MainDictionaryType>
inline Unique<MainDictionaryType> CreateConditionInfo() {
    using DictionaryType = MainDictionaryType;

    /*--- Conditions info ---*/
    auto p_condition_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ConditionInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    return p_condition_info;
}

namespace detail {

/// Typdef.
template<typename TKeySetValuesTypeTag>
using RegistryType = std::unordered_map<std::string, std::function<Unique<Dictionary<TKeySetValuesTypeTag>>()>>;

/// Template declaration.
template<typename TKeySetValuesTypeTag>
inline RegistryType<TKeySetValuesTypeTag> RegisterAll();



/// @brief Registers all possible dictionaries in DictionaryFactory<key::MainValuesTypeTag>.
///        Template specialization for key::MainValuesTypeTag. We may add other specialization,
///        if we decide add other ValuesTypeTag's.
/// @return RegistryType<key::MainValuesTypeTag>
template<>
inline RegistryType<key::MainValuesTypeTag> RegisterAll() {

    RegistryType<key::MainValuesTypeTag> creators;
    creators["Settings"] = &CreateSettings;
    creators["ConditionSettings"] = &CreateConditionSettings;
    creators["ModelInfo"] = &CreateModelInfo;
    creators["ConditionInfo"] = &CreateConditionInfo;

    return creators;
}

} // End namesapce detail.
} // End namespace factories.

///@name QuESo Classes
///@{

/// @brief Provides factory methods for all available dictionaries.
/// @tparam TKeySetValuesTypeTag
template<typename TKeySetValuesTypeTag>
class DictionaryFactory {
public:
    ///@name Type definitions
    ///@{

    using DictionaryType = Dictionary<TKeySetValuesTypeTag>;
    using CreatorFunc = std::function<Unique<DictionaryType>()>;

    ///@}
    ///@name Operations
    ///@{

    static Unique<DictionaryType> Create(const std::string& rName){
        auto it = mRegistry.find(rName);
        QuESo_ERROR_IF( it == mRegistry.end() ) << "Dictionary: '" << rName << "' is not registered.\n"
            << "Registered dictionaries: " << GetPossibleDictionaryNames() << "\n";
        return it->second();
    }

private:
    ///@}
    ///@name Private operations
    ///@{

    static std::string GetPossibleDictionaryNames() {
        std::ostringstream oss;
        oss << '[';
        bool first = true;
        for (const auto& [key, _] : mRegistry) {
            if (!first) oss << ", ";
            oss << '\'' << key << '\'';
            first = false;
        }
        oss << ']';
        return oss.str();
    }

    ///@}
    ///@name Private members
    ///@{

    inline static const std::unordered_map<std::string, CreatorFunc> mRegistry = factories::detail::RegisterAll<TKeySetValuesTypeTag>();
    ///@}
}; // End class DictionaryFactory.
///@}

} // End namespace queso.

#endif // DEFINE_DICTIONARY_FACTORY_HPP
