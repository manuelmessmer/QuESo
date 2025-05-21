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
#include "queso/containers/dictionary_tmp.hpp"

namespace queso {

using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

namespace factories {

Unique<MainDictionaryType> CreateSettings() {
    using DictionaryType = MainDictionaryType;

    auto p_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        DictionaryType::EmptyKeySetType,
        key::detail::MainSettingsSubDictTypeTagKeySetInfo,
        key::detail::MainSettingsListTypeTagKeySetInfo>{} );

    auto p_general_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::GeneralSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_general_settings->SetValue(GeneralSettings::output_directory_name, std::string("queso_output"));
    p_general_settings->SetValue(GeneralSettings::echo_level, 1u);
    p_general_settings->SetValue(GeneralSettings::write_output_to_file, true);

    p_settings->SetSubDictionary(MainSettings::general_settings, std::move(p_general_settings));

    auto p_background_grid_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::BackgroundGridInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_settings->SetSubDictionary(MainSettings::background_grid_settings, std::move(p_background_grid_settings));

    auto p_trimmed_quadrature_rule_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::TrimmedQuadratureRuleSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 1.0e-10);
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 1.0e-3);
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, 100u);
    p_trimmed_quadrature_rule_settings->SetValue(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed, true);

    p_settings->SetSubDictionary(MainSettings::trimmed_quadrature_rule_settings, std::move(p_trimmed_quadrature_rule_settings));

    auto p_non_trimmed_quadrature_rule_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::NonTrimmedQuadratureRuleSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_non_trimmed_quadrature_rule_settings->SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::gauss);

    p_settings->SetSubDictionary(MainSettings::non_trimmed_quadrature_rule_settings, std::move(p_non_trimmed_quadrature_rule_settings));

    return p_settings;
}

Unique<MainDictionaryType> CreateConditionSettings() {
    using DictionaryType = MainDictionaryType;

    auto p_condition_settings = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ConditionSettingsMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    return p_condition_settings;
}

Unique<MainDictionaryType> CreateModelInfo() {
    using DictionaryType = MainDictionaryType;

    auto p_model_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        DictionaryType::EmptyKeySetType,
        key::detail::MainInfoSubDictTypeTagKeySetInfo,
        key::detail::MainInfoListTypeTagKeySetInfo>{} );

    // EmbeddedGeometryInfo subdictionary
    auto p_embedded_geometry_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::EmbeddedGeometryInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_model_info->SetSubDictionary(MainInfo::embedded_geometry_info, std::move(p_embedded_geometry_info));

    // QuadratureInfo subdictionary
    auto p_quadrature_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::QuadratureInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_model_info->SetSubDictionary(MainInfo::quadrature_info, std::move(p_quadrature_info));

    // BackgroundGridInfo subdictionary
    auto p_background_grid_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::BackgroundGridInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_model_info->SetSubDictionary(MainInfo::background_grid_info, std::move(p_background_grid_info));

    // ElapsedTimeInfo subdictionary
    auto p_elapsed_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ElapsedTimeInfoMainValuesTypeTagKeySetInfo,
        key::detail::ElapsedTimeInfoSubDictTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType>{} );

    // VolumeTimeInfo subdictionary
    auto p_volume_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::VolumeTimeInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_elapsed_time_info->SetSubDictionary(ElapsedTimeInfo::volume_time_info, std::move(p_volume_time_info));

    // ConditionsTimeInfo subdictionary
    auto p_conditions_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ConditionsTimeInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_elapsed_time_info->SetSubDictionary(ElapsedTimeInfo::conditions_time_info, std::move(p_conditions_time_info));

    // WriteFilesTimeInfo subdictionary
    auto p_write_files_time_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::WriteFilesTimeInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    p_elapsed_time_info->SetSubDictionary(ElapsedTimeInfo::write_files_time_info, std::move(p_write_files_time_info));

    p_model_info->SetSubDictionary(MainInfo::elapsed_time_info, std::move(p_elapsed_time_info));

    return p_model_info;
}

Unique<MainDictionaryType> CreateConditionInfo() {
    using DictionaryType = MainDictionaryType;

    auto p_condition_info = MakeUnique<DictionaryType>( DictionaryType::KeySetInfosTypeTag<
        key::detail::ConditionInfoMainValuesTypeTagKeySetInfo,
        DictionaryType::EmptyKeySetType,
        DictionaryType::EmptyKeySetType>{} );

    return p_condition_info;
}
namespace detail {

template<typename TKeySetValuesTypeTag>
using RegistryType = std::unordered_map<std::string, std::function<Unique<Dictionary<TKeySetValuesTypeTag>>()>>;

template<typename TKeySetValuesTypeTag>
RegistryType<TKeySetValuesTypeTag> RegisterAll();

template<>
RegistryType<key::MainValuesTypeTag> RegisterAll() {
    using DictionaryType = MainDictionaryType;

    RegistryType<key::MainValuesTypeTag> creators;
    creators["Settings"] = &CreateSettings;
    creators["ConditionSettings"] = &CreateConditionSettings;
    creators["ModelInfo"] = &CreateModelInfo;
    creators["ConditionInfo"] = &CreateConditionInfo;

    return creators;
}

} // End namesapce detail
} // End namespace factories

template<typename TKeySetValuesTypeTag>
class DictionaryFactory {
public:
    using DictionaryType = Dictionary<TKeySetValuesTypeTag>;
    using CreatorFunc = std::function<Unique<Dictionary>()>;

    static Unique<DictionaryType> Create(const std::string& rName){
        auto it = mRegistry.find(rName);
        QuESo_ERROR_IF( it == mRegistry.end() ) << "Dictionary: '" << rName << "' is not registered.\n"
            << "Registered dictionaries: " << GetPossibleDictionaryNames() << "\n";
        return it->second();
    }

private:

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

    inline static const std::unordered_map<std::string, CreatorFunc> mRegistry = factories::detail::RegisterAll<TKeySetValuesTypeTag>();
};

}

#endif // DEFINE_DICTIONARY_FACTORY_HPP
