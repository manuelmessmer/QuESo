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

#ifndef EMBEDDED_MODEL_INCLUDE_H
#define EMBEDDED_MODEL_INCLUDE_H

/// STL includes

/// Project includes
#include "queso/containers/triangle_mesh_interface.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/io/io_utilities.h"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/utilities/check_dictionary_utilities.hpp"
#include "queso/containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

////
/**
 * @class  EmbeddedModel
 * @author Manuel Messmer
 * @brief  Main class of QuESo.
 *         Provides interface to create integration points for embedded volumes and embedded conditions (boundaries)
 *         based on the settings provided in mSettings. EmbeddedModel aggregates a BackgroundGrid, which stores the
 *         created active elements and conditions. Each element stores its integrations points. Each condition is
 *         split into ConditionSegment's (that conform to the boundaries of the elements in the background grid).
 *         The corresponding boundary integration points are stored on these ConditionSegments.
 *         EmbeddedModel also stores some information regearding the created model in mModelInfo.
**/
class EmbeddedModel
{
public:
    ///@name Type Definitions
    ///@{

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using ElementType = Element<IntegrationPointType, BoundaryIntegrationPointType>;
    using BackgroundGridType = BackgroundGrid<ElementType>;
    using ConditionType = Condition<ElementType>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Helper to create EmbeddedModel.
    /// @param pSettings (EmbeddedModel takes unique ownership).
    /// @return EmbeddedModel.
    static EmbeddedModel Create(Unique<MainDictionaryType> pSettings) {
        CheckDictionaryUtilities::CheckSettings(*pSettings);
        return EmbeddedModel(std::move(pSettings));
    }

private:
    /// @brief Constructor
    /// @param pSettings (EmbeddedModel takes unique ownership).
    EmbeddedModel(Unique<MainDictionaryType> pSettings) :
        mpSettings(std::move(pSettings)),
        mGridIndexer(*mpSettings),
        mBackgroundGrid(*mpSettings),
        mpModelInfo(DictionaryFactory<key::MainValuesTypeTag>::Create("ModelInfo"))
    {
    }

public:
    /// Copy Constructor
    EmbeddedModel(const EmbeddedModel &rOther) = delete;
    /// Copy Assignement
    EmbeddedModel& operator= (const EmbeddedModel &rOther) = delete;
    /// Move constructor
    EmbeddedModel(EmbeddedModel&& rOther) noexcept = default;
    /// Move assignement operator
    EmbeddedModel& operator=(EmbeddedModel&& rOther) noexcept = default;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Main function to create embedded model.
    ///       Creates integration points for both the embedded volume and all embedded conditions.
    ///       The respective geometries (TriangleMeshes) are taken from input STL files specified in mSettings.
    ///@todo Add try{} catch{} plus error handler
    void CreateAllFromSettings() {
        // Create volume
        const auto& r_general_settings = GetSettings()[MainSettings::general_settings];
        const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
        QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Volume -------------------------------------- START" << std::endl;

        const auto& r_filename = r_general_settings.GetValue<std::string>(GeneralSettings::input_filename);
        TriangleMesh triangle_mesh{};
        IO::ReadMeshFromSTL(triangle_mesh, r_filename.c_str());

        ComputeVolume(triangle_mesh);
        PrintVolumeElapsedTimeInfo();

        QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Volume ---------------------------------------- End\n";

        // Create conditions
        const auto& r_conditions_settings_list = GetSettings().GetList(MainSettings::conditions_settings_list);
        if( r_conditions_settings_list.size() > 0 ){

            QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Conditions ---------------------------------- START" << std::endl;
            for( const auto& p_condition_settings : r_conditions_settings_list ){
                const auto& r_filename = p_condition_settings->GetValue<std::string>(ConditionSettings::input_filename);
                TriangleMesh triangle_mesh{};
                IO::ReadMeshFromSTL(triangle_mesh, r_filename.c_str());
                ComputeCondition(triangle_mesh, *p_condition_settings);
            }
            PrintConditionsElapsedTimeInfo();
            QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Conditions ------------------------------------ End" << std::endl;
        }

        QuESo_INFO_IF(echo_level > 0) << "QuESo: Write Model To File -------------------------------- START\n";
        WriteModelToFile();
        QuESo_INFO_IF(echo_level > 0) << "QuESo: Write Model To File ---------------------------------- End\n" << std::endl;
    }

    ///@brief Creates integration points for an embedded volume that is enclosed/defined by rTriangleMesh.
    ///       This interface enables to pass a TriangleMeshInterface and, hence, facilitates other applications to
    ///       use QuESo on C++ level, which do not want QuESo to read rTriangleMesh from an input file.
    ///@param rTriangleMesh
    ///@see CreateAllFromSettings <- Creates volume and condition directly from input files specified in mSettings.
    ///@todo Add try{} catch{} plus error handler
    void CreateVolume(const TriangleMeshInterface& rTriangleMesh){
        ComputeVolume(rTriangleMesh);
    }

    ///@brief Creates integration points for an embedded condition defined by rTriangleMesh.
    ///       This interface enables to pass a TriangleMeshInterface and, hence, facilitates other applications to
    ///       use QuESo on C++ level, which do not want QuESo to read rTriangleMesh from an input file.
    ///@param rTriangleMesh
    ///@param rConditionSettings
    ///@see CreateAllFromSettings <- Creates volume and condition directly from input files specified in mSettings.
    ///@todo Add try{} catch{} plus error handler
    void CreateCondition(const TriangleMeshInterface& rTriangleMesh, const MainDictionaryType& rConditionSettings){
        ComputeCondition(rTriangleMesh, rConditionSettings);
    }

    /// @brief Writes this model to file.
    ///        Elements and integrations points are written to VTK files.
    ///        Conditions are written to STL files.
    ///        mpModelInfo is written to JSON file.
    void WriteModelToFile() const;

    /// @brief Returns all active elements.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType::ElementContainerType& GetElements() const {
        return mBackgroundGrid.GetElements();
    }

    /// @brief Returns all conditions.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType::ConditionContainerType& GetConditions() const {
        return mBackgroundGrid.GetConditions();
    }

    ///@brief Returns the Settings
    ///@return const MainDictionaryType&
    const MainDictionaryType& GetSettings() const {
        return *mpSettings;
    }

    ///@brief Returns the ModelInfo (const version).
    ///@return const MainDictionaryType&
    const MainDictionaryType& GetModelInfo() const {
        return *mpModelInfo;
    }

    ///@brief Returns the ModelInfo (const version).
    ///@return const MainDictionaryType&
    MainDictionaryType& GetModelInfo() {
        return *mpModelInfo;
    }

    ///@}
private:

    ///@name Private Member Operations
    ///@{


    ///@brief Returns the ModelInfo.
    ///@return const MainDictionaryType&
    MainDictionaryType& GetModelInfoMutable() const {
        return *mpModelInfo;
    }

    ///@brief Main function to compute the integration points for a volume enclosed/defined by rTriangleMesh.
    ///@param rTriangleMesh
    void ComputeVolume(const TriangleMeshInterface& rTriangleMesh);

    ///@brief Main function to compute the integration points for a condition defined by rTriangleMesh.
    ///@param rTriangleMesh
    ///@param rConditionSettings
    void ComputeCondition(const TriangleMeshInterface& rTriangleMesh, const MainDictionaryType& rConditionSettings);

    ///@brief Prints a warning, if the rTriangleMesh is not fully contained within the bounding box defined
    ///       by 'lower_bound_xyz' and 'upper_bound_xyz' in mSettings.
    ///@param rTriangleMesh
    void CheckIfMeshIsWithinBoundingBox(const TriangleMeshInterface& rTriangleMesh) const;

    ///@brief Prints some info to the console regarding the computed volume.
    ///       Since only one volume per EmbeddedModel can be created no arguments have to be passed.
    void PrintVolumeInfo() const;

    ///@brief Prints some info to the console regarding the computed condition.
    ///@param rConditionInfo Info to be printed.
    void PrintConditionInfo(const MainDictionaryType& rConditionInfo) const;

    ///@brief Prints some info to the console regarding the elapsed computing times required to create the volume.
    ///       Info to be printed is taken from mModelInfo[MainInfo::elapsed_time_info].
    void PrintVolumeElapsedTimeInfo() const;

    ///@brief Prints some info to the console regarding the elapsed computing times required to create the conditions.
    ///       Info to be printed is taken from mModelInfo[MainInfo::elapsed_time_info].
    void PrintConditionsElapsedTimeInfo() const;

    ///@}
    ///@name Private Members Variables
    ///@{
    Unique<const MainDictionaryType> mpSettings;
    const GridIndexer mGridIndexer;
    BackgroundGridType mBackgroundGrid;
    Unique<MainDictionaryType> mpModelInfo;
    ///@}
};
///@} End QuESo Classes
} // End namespace queso

#endif // EMBEDDED_MODEL_INCLUDE_H
