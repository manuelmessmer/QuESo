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
#include "queso/containers/grid_indexer.hpp"
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
    typedef Condition<ElementType> ConditionType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rParameters
    QuESo(const Settings& rSettings ) : mSettings(rSettings), mGridIndexer(mSettings) {

        const_cast<Settings&>(mSettings).Check();
        mpTriangleMesh = MakeUnique<TriangleMesh>();
        const auto& r_general_settings = mSettings[MainSettings::general_settings];

        // Read mesh
        const auto& r_filename = r_general_settings.GetValue<std::string>(GeneralSettings::input_filename);
        IO::ReadMeshFromSTL(*mpTriangleMesh, r_filename.c_str());

        Check();
    }

    /// Copy Constructor
    QuESo(const QuESo &m) = delete;

    /// Copy Assignement
    QuESo & operator= (const QuESo &) = delete;

    /// Move constructor
    QuESo(QuESo&& rOther) = delete;
    /// Move assignement operator
    QuESo& operator=(QuESo&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Run QuESo.
    void Run();

    /// @brief Get all active elements.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType::ElementContainerType& GetElements() const {
        return mpBackgroundGrid->GetElements();
    }

    /// @brief Get all conditions.
    /// @return const Reference to ElementVectorPtrType
    const BackgroundGridType::ConditionContainerType& GetConditions() const {
        return mpBackgroundGrid->GetConditions();
    }


    /// @brief Performs necessary checks.
    void Check() const;

    /// @brief Clear all containers.
    void Clear() {
        mpTriangleMesh->Clear();
        mpBRepOperator = nullptr;
        mpBackgroundGrid = nullptr;
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
    Unique<BackgroundGridType> mpBackgroundGrid;
    const Settings mSettings;
    GridIndexer mGridIndexer;
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