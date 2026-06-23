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

#pragma once

/// Project includes
#include "queso/containers/background_grid.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/utilities/check_dictionary_utilities.hpp"

namespace queso {

///@name QuESo Classes
///@{

////
/**
 * @class  EmbeddedComponent
 * @author Manuel Messmer
 * @brief  Main class of QuESo.
 *         Provides interface to create integration points for embedded volumes and embedded conditions (boundaries)
 *         based on the settings provided in mSettings. EmbeddedComponent aggregates a BackgroundGrid, which stores the
 *         created active elements and conditions. Each active element stores its integrations points. Each condition is
 *         split into ConditionSegment's (that conform to the boundaries of the elements in the background grid).
 *         The corresponding boundary integration points are stored on these ConditionSegments.
 *         EmbeddedComponent also stores some information regearding the created component in mComponentInfo.
 **/
class EmbeddedComponent
{
public:
    ///@name Type Definitions
    ///@{

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>;
    using ElementViewType = typename BackgroundGridType::ElementViewType;
    using ConditionType = Condition<ElementViewType>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Helper to create EmbeddedComponent.
    /// @param pSettings (EmbeddedComponent takes unique ownership).
    /// @return EmbeddedComponent.
    static EmbeddedComponent Create(Unique<MainDictionaryType>&& pSettings)
    {
        CheckDictionaryUtilities::CheckSettings(*pSettings);
        return EmbeddedComponent(std::move(pSettings));
    }

private:
    /// @brief Constructor
    /// @param pSettings (EmbeddedComponent takes unique ownership).
    EmbeddedComponent(Unique<MainDictionaryType>&& pSettings)
        : mpSettings(std::move(pSettings)), mGridIndexer(*mpSettings), mBackgroundGrid(*mpSettings),
          mpComponentInfo(DictionaryFactory<key::MainValuesTypeTag>::Create("ComponentInfo"))
    {}

public:
    /// Copy Constructor
    EmbeddedComponent(const EmbeddedComponent& rOther) = delete;
    /// Copy Assignement
    EmbeddedComponent& operator=(const EmbeddedComponent& rOther) = delete;
    /// Move constructor
    EmbeddedComponent(EmbeddedComponent&& rOther) noexcept = default;
    /// Move assignement operator
    EmbeddedComponent& operator=(EmbeddedComponent&& rOther) noexcept = default;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Main function to create embedded component.
    ///       Creates integration points for both the embedded volume and all embedded conditions.
    ///       The respective geometries (TriangleMeshes) are taken from input STL files specified in mSettings.
    ///@todo Add try{} catch{} plus error handler
    void CreateAllFromSettings()
    {
        // Create volume
        const IndexType echo_level = GetSettings().GetRequiredValue<IndexType>(MainSettings::echo_level);
        QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Volume -------------------------------------- START"
                                      << std::endl;

        const auto& r_filename = GetSettings().GetRequiredValue<std::string>(MainSettings::input_filename);
        TriangleMesh triangle_mesh{};
        IO::ReadMeshFromSTL(triangle_mesh, r_filename.c_str());

        ComputeVolume(triangle_mesh.View());
        PrintVolumeElapsedTimeInfo();

        QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Volume ---------------------------------------- End\n";

        // Create conditions
        const auto& r_conditions_settings_list = GetSettings().GetList(MainSettings::conditions_settings_list);
        if (r_conditions_settings_list.size() > 0) {

            QuESo_INFO_IF(echo_level > 0)
                << "QuESo: Create Conditions ---------------------------------- START" << std::endl;
            for (const auto& p_condition_settings : r_conditions_settings_list) {
                const auto& r_filename_cond =
                    p_condition_settings->GetRequiredValue<std::string>(ConditionSettings::input_filename);
                TriangleMesh triangle_meshs_cond{};
                IO::ReadMeshFromSTL(triangle_meshs_cond, r_filename_cond.c_str());
                ComputeCondition(triangle_meshs_cond.View(), *p_condition_settings);
            }
            PrintConditionsElapsedTimeInfo();
            QuESo_INFO_IF(echo_level > 0)
                << "QuESo: Create Conditions ------------------------------------ End" << std::endl;
        }

        QuESo_INFO_IF(echo_level > 0) << "QuESo: Write Model To File -------------------------------- START\n";
        WriteModelToFile();
        QuESo_INFO_IF(echo_level > 0) << "QuESo: Write Model To File ---------------------------------- End\n"
                                      << std::endl;
    }

    ///@brief Creates integration points for an embedded volume that is enclosed/defined by rTriangleMesh.
    ///       This interface enables to pass a TriangleMeshView and, hence, facilitates other applications to
    ///       use QuESo on C++ level, which do not want QuESo to read rTriangleMesh from an input file.
    ///@param rTriangleMesh
    ///@see CreateAllFromSettings <- Creates volume and condition directly from input files specified in mSettings.
    ///@todo Add try{} catch{} plus error handler
    void CreateVolume(const TriangleMeshView& rTriangleMesh)
    { ComputeVolume(rTriangleMesh); }

    ///@brief Creates integration points for an embedded condition defined by rTriangleMesh.
    ///       This interface enables to pass a TriangleMeshView and, hence, facilitates other applications to
    ///       use QuESo on C++ level, which do not want QuESo to read rTriangleMesh from an input file.
    ///@param rTriangleMesh
    ///@param rConditionSettings
    ///@see CreateAllFromSettings <- Creates volume and condition directly from input files specified in mSettings.
    ///@todo Add try{} catch{} plus error handler
    void CreateCondition(const TriangleMeshView& rTriangleMesh, const MainDictionaryType& rConditionSettings)
    { ComputeCondition(rTriangleMesh, rConditionSettings); }

    /// @brief Writes this component to file.
    ///        Elements and integrations points are written to VTK files.
    ///        Conditions are written to STL files.
    ///        mpComponentInfo is written to JSON file.
    void WriteModelToFile() const;

    /// @brief Returns a lazy range of ElementView over all matching elements.
    /// @details TFilter: all -> both, trimmed -> trimmed only, untrimmed -> untrimmed only.
    template<BackgroundGridType::ElementFilter TFilter = BackgroundGridType::ElementFilter::all>
    auto GetElementViews() const
    { return mBackgroundGrid.template GetElementViews<TFilter>(); }

    /// @brief Returns a span over the raw element container matching TFilter.
    /// @details ElementFilter::all is not valid - use GetElementViews() instead.
    template<BackgroundGridType::ElementFilter TFilter>
    [[nodiscard]] auto GetElements() const noexcept
    { return mBackgroundGrid.template GetElements<TFilter>(); }

    /// @brief Returns all conditions.
    /// @return const Reference to ConditionContainerType
    const BackgroundGridType::ConditionContainerType& GetConditions() const
    { return mBackgroundGrid.GetConditions(); }

    ///@brief Returns the Settings
    ///@return const MainDictionaryType&
    const MainDictionaryType& GetSettings() const
    { return *mpSettings; }

    ///@brief Returns the ComponentInfo (const version).
    ///@return const MainDictionaryType&
    const MainDictionaryType& GetComponentInfo() const
    { return *mpComponentInfo; }

    ///@brief Returns the ComponentInfo (non-const version).
    ///@return MainDictionaryType&
    MainDictionaryType& GetComponentInfo()
    { return *mpComponentInfo; }

    ///@}
private:
    ///@name Private Member Operations
    ///@{

    ///@brief Returns the ComponentInfo as non-const reference. May be called from const member funtions.
    ///@return MainDictionaryType&
    MainDictionaryType& GetComponentInfoMutable() const
    { return *mpComponentInfo; }

    ///@brief Main function to compute the integration points for a volume enclosed/defined by rTriangleMesh.
    ///@param rTriangleMesh
    void ComputeVolume(const TriangleMeshView& rTriangleMesh);

    ///@brief Main function to compute the integration points for a condition defined by rTriangleMesh.
    ///@param rTriangleMesh
    ///@param rConditionSettings
    void ComputeCondition(const TriangleMeshView& rTriangleMesh, const MainDictionaryType& rConditionSettings);

    ///@brief Prints a warning, if the rTriangleMesh is not fully contained within the bounding box defined
    ///       by 'lower_bound_xyz' and 'upper_bound_xyz' in mSettings.
    ///@param rTriangleMesh
    void CheckIfMeshIsWithinBoundingBox(const TriangleMeshView& rTriangleMesh) const;

    ///@brief Prints some info to the console regarding the computed volume.
    ///       Since only one volume per EmbeddedComponent can be created no arguments have to be passed.
    void PrintVolumeInfo() const;

    ///@brief Prints some info to the console regarding the computed condition.
    ///@param rConditionInfo Info to be printed.
    void PrintConditionInfo(const MainDictionaryType& rConditionInfo) const;

    ///@brief Prints some info to the console regarding the elapsed computing times required to create the volume.
    ///       Info to be printed is taken from mComponentInfo[MainInfo::elapsed_time_info].
    void PrintVolumeElapsedTimeInfo() const;

    ///@brief Prints some info to the console regarding the elapsed computing times required to create the conditions.
    ///       Info to be printed is taken from mComponentInfo[MainInfo::elapsed_time_info].
    void PrintConditionsElapsedTimeInfo() const;

    ///@}
    ///@name Private Members Variables
    ///@{
    Unique<const MainDictionaryType> mpSettings;
    GridIndexer mGridIndexer;
    BackgroundGridType mBackgroundGrid;
    Unique<MainDictionaryType> mpComponentInfo;
    ///@}
};
///@} End QuESo Classes
}// End namespace queso
