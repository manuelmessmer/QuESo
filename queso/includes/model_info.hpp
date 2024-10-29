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

#ifndef MODEL_INFO_INCLUDE_HPP
#define MODEL_INFO_INCLUDE_HPP

/// Project includes
#include "queso/containers/dictionary.hpp"

namespace queso {

/// Definition of ModelInfo keys
enum class RootInfo {main_info=DictStarts::start_subdicts};
enum class MainInfo {
    embedded_geometry_info=DictStarts::start_subdicts, quadrature_info, background_grid_info, elapsed_time_info,
    conditions_infos_list=DictStarts::start_lists};
enum class EmbeddedGeometryInfo {
    is_closed=DictStarts::start_values, volume};
enum class QuadratureInfo {
    represented_volume=DictStarts::start_values, percentage_of_geometry_volume, tot_num_points, num_of_points_per_full_element, num_of_points_per_trimmed_element};
enum class BackgroundGridInfo {
    num_active_elements=DictStarts::start_values, num_trimmed_elements, num_full_elements, num_inactive_elements};
enum class ConditionInfo {
    condition_id=DictStarts::start_values, surf_area, perc_surf_area_in_active_domain};

enum class ElapsedTimeInfo {
    total=DictStarts::start_values,
    volume_time_info=DictStarts::start_subdicts, conditions_time_info, write_files_time_info,
    };
enum class VolumeTimeInfo {
    total=DictStarts::start_values, classification_of_elements, computation_of_intersections, solution_of_moment_fitting_eqs, construction_of_ggq_rules};
enum class ConditionsTimeInfo {
    total=DictStarts::start_values};
enum class WriteFilesTimeInfo {
    total=DictStarts::start_values};

typedef Dictionary<RootInfo, MainInfo, EmbeddedGeometryInfo, QuadratureInfo, BackgroundGridInfo, ConditionInfo,
    ElapsedTimeInfo, VolumeTimeInfo, ConditionsTimeInfo, WriteFilesTimeInfo> ModelInfoBaseType;

///@name QuESo Classes
///@{

/**
 * @class  ModelInfo.
 * @author Manuel Messmer.
 * @brief  Derived version from Dictionary, which specifies all available Keys, ValueTypes and default values.
 * @see    Dictionary.
**/
class ModelInfo : public ModelInfoBaseType
{
public:
    ///@}
    ///@name Type definitions
    ///@{

    typedef ModelInfoBaseType BaseType;
    typedef std::string Str;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor. Sets up the default dictionary.
    ModelInfo() : BaseType(RootInfo::main_info, Str("model_info")) {

        /// Let's define the default model_info...
        bool DontSet = false; // Given values are only dummy values used to deduce the associated type.
                              // However, calling IsSet() will return false.

        /// EmbeddedGeometryInfo
        auto& r_embedded_geometry_info = AddEmptySubDictionary(MainInfo::embedded_geometry_info, Str("embedded_geometry_info"));
        r_embedded_geometry_info.AddValues(std::make_tuple(
            std::make_tuple(EmbeddedGeometryInfo::is_closed, Str("is_closed"), false, DontSet ),
            std::make_tuple(EmbeddedGeometryInfo::volume, Str("volume"), 0.0, DontSet )
        ));

        /// QuadratureRuleInfo
        auto& r_quadrature_info = AddEmptySubDictionary(MainInfo::quadrature_info, Str("quadrature_info"));
        r_quadrature_info.AddValues(std::make_tuple(
            std::make_tuple(QuadratureInfo::represented_volume, Str("represented_volume"), 0.0, DontSet ),
            std::make_tuple(QuadratureInfo::percentage_of_geometry_volume, Str("percentage_of_geometry_volume"), 0.0, DontSet),
            std::make_tuple(QuadratureInfo::tot_num_points, Str("tot_num_points"), IndexType(0), DontSet ),
            std::make_tuple(QuadratureInfo::num_of_points_per_full_element, Str("num_of_points_per_full_element"), 0.0, DontSet ),
            std::make_tuple(QuadratureInfo::num_of_points_per_trimmed_element, Str("num_of_points_per_trimmed_element"), 0.0, DontSet )
        ));

        /// BackgroundGridInfo
        auto& r_background_grid_info = AddEmptySubDictionary(MainInfo::background_grid_info, Str("background_grid_info"));
        r_background_grid_info.AddValues(std::make_tuple(
            std::make_tuple(BackgroundGridInfo::num_active_elements, Str("num_active_elements"), IndexType(0), DontSet ),
            std::make_tuple(BackgroundGridInfo::num_trimmed_elements, Str("num_trimmed_elements"), IndexType(0), DontSet ),
            std::make_tuple(BackgroundGridInfo::num_full_elements, Str("num_full_elements"), IndexType(0), DontSet ),
            std::make_tuple(BackgroundGridInfo::num_inactive_elements, Str("num_inactive_elements"), IndexType(0), DontSet )
        ));

        /// ElapsedTimeInfo
        auto& r_elapsed_time_info = AddEmptySubDictionary(MainInfo::elapsed_time_info, Str("elapsed_time_info"));
        r_elapsed_time_info.AddValues(std::make_tuple(
            std::make_tuple(ElapsedTimeInfo::total, Str("total"), 0.0, DontSet )
        ));

        auto& r_volume_time_info = r_elapsed_time_info.AddEmptySubDictionary(ElapsedTimeInfo::volume_time_info, Str("volume_time_info"));
        r_volume_time_info.AddValues(std::make_tuple(
            std::make_tuple(VolumeTimeInfo::total, Str("total"), 0.0, DontSet ),
            std::make_tuple(VolumeTimeInfo::classification_of_elements, Str("classification_of_elements"), 0.0, DontSet ),
            std::make_tuple(VolumeTimeInfo::computation_of_intersections, Str("computation_of_intersections"), 0.0, DontSet ),
            std::make_tuple(VolumeTimeInfo::solution_of_moment_fitting_eqs, Str("solution_of_moment_fitting_eqs"), 0.0, DontSet),
            std::make_tuple(VolumeTimeInfo::construction_of_ggq_rules, Str("construction_of_ggq_rules"), 0.0, DontSet)
        ));

        auto& r_conditions_time_info = r_elapsed_time_info.AddEmptySubDictionary(ElapsedTimeInfo::conditions_time_info, Str("conditions_time_info"));
        r_conditions_time_info.AddValues(std::make_tuple(
            std::make_tuple(ConditionsTimeInfo::total, Str("total"), 0.0, DontSet )
        ));

        auto& r_write_files_time_info = r_elapsed_time_info.AddEmptySubDictionary(ElapsedTimeInfo::write_files_time_info, Str("write_files_time_info"));
        r_write_files_time_info.AddValues(std::make_tuple(
            std::make_tuple(WriteFilesTimeInfo::total, Str("total"), 0.0, DontSet )
        ));

        /// ConditionInfos
        AddEmptyList(MainInfo::conditions_infos_list, Str("conditions_infos_list"));
    }

    ///@}
    ///@name Member operations
    ///@{

    /// @brief Creates new condition settings.
    /// @return SettingsBaseType& Reference to dictionary that contains condition settings.
    ModelInfoBaseType& CreateNewConditionInfo() {
        bool DontSet = false; // Given values are only dummy values used to deduced the associated type.

        auto& r_conditions_infos = GetListObject(MainInfo::conditions_infos_list);
        auto& r_new_condition_info = r_conditions_infos.AddListItem(std::make_tuple(
            std::make_tuple(ConditionInfo::condition_id, Str("condition_id"), IndexType(0), DontSet ),
            std::make_tuple(ConditionInfo::surf_area, Str("surf_area"), 0.0, DontSet ),
            std::make_tuple(ConditionInfo::perc_surf_area_in_active_domain, Str("perc_surf_area_in_active_domain"), 0.0, DontSet  )
        ));

        return r_new_condition_info;
    }

private:

    /// Hide the following functions
    /// @note Note making these functions protected in the Base class does not help, since they are called here on
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

#endif // End MODEL_INFO_INCLUDE_HPP