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

#ifndef QuESo_INCLUDE_H
#define QuESo_INCLUDE_H

/// STL includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>

/// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/condition.hpp"
#include "queso/utilities/mapping_utilities.h"
#include "queso/includes/settings.hpp"
#include "queso/embedding/brep_operator.h"

namespace queso {

///@name QuESo Classes
///@{

////
/**
 * @class  QuESo
 * @author Manuel Messmer
 * @brief  Main class of QuESo.
*/
class QuESo
{
public:
    ///@name Type Definitions
    ///@{

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
    typedef BackgroundGrid<ElementType> BackgroundGridType;

    typedef std::vector<ElementType> ElementVectorType;
    typedef std::vector<Condition> ConditionVectorType;
    typedef std::vector<Unique<BRepOperator>> BRepOperatorPtrVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rParameters
    QuESo(const Settings& rSettings ) : mSettings(rSettings), mMapper(mSettings) {

        const_cast<Settings&>(mSettings).Check();
        mpTriangleMesh = MakeUnique<TriangleMesh>();
        const auto& r_general_settings = mSettings[MainSettings::general_settings];

        // Read mesh
        const auto& r_filename = r_general_settings.GetValue<std::string>(GeneralSettings::input_filename);
        IO::ReadMeshFromSTL(*mpTriangleMesh, r_filename.c_str());

        const auto& r_conditions_settings_list =  mSettings[MainSettings::conditions_settings_list];
        const IndexType number_of_conditions = r_conditions_settings_list.NumberOfSubDictionaries();
        for( IndexType i = 0; i < number_of_conditions; ++i){
            CreateNewCondition(r_conditions_settings_list[i]);
        }

        Check();
    }

    /// Copy Constructor
    QuESo(const QuESo &m) = delete;

    /// Copy Assignement
    QuESo & operator= (const QuESo &) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Run QuESo.
    void Run();

    /// @brief Get all active elements.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType& GetElements() const {
        return *mpBackgroundGrid;
    }

    /// @brief Get all conditions.
    /// @return const Reference to ElementVectorPtrType
    const ConditionVectorType& GetConditions() const {
        return mConditions;
    }

    /// @brief Creates a new condition and adds to mConditions-Vector.
    /// @param rConditionSettings
    /// @return Condition&
    Condition& CreateNewCondition(const SettingsBaseType& rConditionSettings);

    /// @brief Performs necessary checks.
    void Check() const;

    /// @brief Clear all containers.
    void Clear() {
        mpTriangleMesh->Clear();
        mpBRepOperator = nullptr;
        mpBrepOperatorsBC.clear();
        mpBackgroundGrid = nullptr;
        mConditions.clear();
    }

    /// @brief Returns triangle mesh.
    /// @return const TriangleMesh&
    const TriangleMeshInterface& GetTriangleMesh() {
        return *mpTriangleMesh;
    }

    ///@}

private:

    ///@name Private Members Variables
    ///@{
    Unique<TriangleMeshInterface> mpTriangleMesh;
    Unique<BRepOperator> mpBRepOperator;
    BRepOperatorPtrVectorType mpBrepOperatorsBC;
    Unique<BackgroundGridType> mpBackgroundGrid;
    ConditionVectorType mConditions;
    const Settings mSettings;
    Mapper mMapper;
    ///@}

    ///@name Private Member Operations
    ///@{

    /// @brief Compute
    std::array<double,5> Compute();
    ///@}
};

///@}

} // End namespace queso

#endif // QuESo_INCLUDE_H