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
#include "queso/containers/background_grid.hpp"
#include "queso/io/io_utilities.h"
#include "queso/includes/model_info.hpp"

namespace queso {

///@name QuESo Classes
///@{

////
/**
 * @class  EmbeddedModel
 * @author Manuel Messmer
 * @brief  Main class of QuESo.
*/
class EmbeddedModel
{
public:
    ///@name Type Definitions
    ///@{

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
    typedef BackgroundGrid<ElementType> BackgroundGridType;
    typedef Condition<ElementType> ConditionType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rParameters
    EmbeddedModel(const Settings& rSettings) :
        mSettings(rSettings.Check()),
        mGridIndexer(mSettings[MainSettings::background_grid_settings]),
        mBackgroundGrid(mSettings[MainSettings::background_grid_settings]),
        mModelInfo{}
    {
    }

    /// Copy Constructor
    EmbeddedModel(const EmbeddedModel &m) = delete;
    /// Copy Assignement
    EmbeddedModel & operator= (const EmbeddedModel &) = delete;
    /// Move constructor
    EmbeddedModel(EmbeddedModel&& rOther) = delete;
    /// Move assignement operator
    EmbeddedModel& operator=(EmbeddedModel&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void CreateFromSettings() {
        try {
            // Create volume
            const auto& r_general_settings = mSettings[MainSettings::general_settings];
            const IndexType echo_level = r_general_settings.GetValue<IndexType>(GeneralSettings::echo_level);
            QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Volume -------------------------------------- START" << std::endl;

            const auto& r_filename = r_general_settings.GetValue<std::string>(GeneralSettings::input_filename);
            TriangleMesh triangle_mesh{};
            IO::ReadMeshFromSTL(triangle_mesh, r_filename.c_str());

            ComputeVolume(triangle_mesh);
            PrintVolumeElapsedTimeInfo();

            QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Volume ---------------------------------------- End\n";

            // Create conditions
            const auto& r_conditions_settings_list = mSettings.GetList(MainSettings::conditions_settings_list);
            if( r_conditions_settings_list.size() > 0 ){

                QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Conditions ---------------------------------- START" << std::endl;
                for( const auto& r_condition_settings : r_conditions_settings_list ){
                    const auto& r_filename = r_condition_settings.GetValue<std::string>(ConditionSettings::input_filename);
                    TriangleMesh triangle_mesh{};
                    IO::ReadMeshFromSTL(triangle_mesh, r_filename.c_str());
                    ComputeCondition(triangle_mesh, r_condition_settings);
                }
                PrintConditionsElapsedTimeInfo();
                QuESo_INFO_IF(echo_level > 0) << "QuESo: Create Conditions ------------------------------------ End" << std::endl;
            }
        }
        catch (const Exception& rException ){
            std::cerr << rException.what();
        }
    }

    void CreateVolume(const TriangleMeshInterface& rTriangleMesh){
        try {
            ComputeVolume(rTriangleMesh);
        }
        catch (const Exception& rException ){
            std::cerr << rException.what();
        }
    }

    void CreateCondition(const TriangleMeshInterface& rTriangleMesh, const SettingsBaseType& rConditionSettings){
        try {
            ComputeCondition(rTriangleMesh, rConditionSettings);
        }
        catch (const Exception& rException ){
            std::cerr << rException.what();
        }
    }

    /// @brief Get all active elements.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType::ElementContainerType& GetElements() const {
        return mBackgroundGrid.GetElements();
    }

    /// @brief Get all conditions.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType::ConditionContainerType& GetConditions() const {
        return mBackgroundGrid.GetConditions();
    }

    const ModelInfo& GetModelInfo() const {
        return mModelInfo;
    }

    ///@}

private:

    ///@name Private Member Operations
    ///@{

    void ComputeVolume(const TriangleMeshInterface& rTriangleMesh);

    void ComputeCondition(const TriangleMeshInterface& rTriangleMesh, const SettingsBaseType& rConditionSettings);

    void CheckIfMeshIsWithinBoundingBox(const TriangleMeshInterface& rTriangleMesh) const;

    void PrintVolumeInfo() const;

    void PrintConditionInfo(const ModelInfoBaseType& rConditionInfo) const;

    void PrintVolumeElapsedTimeInfo() const;

    void PrintConditionsElapsedTimeInfo() const;

    ///@}
    ///@name Private Members Variables
    ///@{
    const Settings mSettings;
    const GridIndexer mGridIndexer;
    BackgroundGridType mBackgroundGrid;
    ModelInfo mModelInfo;
    ///@}

};

///@}

} // End namespace queso

#endif // QuESo_INCLUDE_H