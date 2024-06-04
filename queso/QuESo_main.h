//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//			    See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
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
#include "queso/containers/element_container.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/condition.hpp"
#include "queso/utilities/mapping_utilities.h"
#include "queso/includes/parameters.h"
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
    typedef ElementContainer<ElementType> ElementContainerType;

    typedef std::vector<ElementType> ElementVectorType;
    typedef std::vector<Condition> ConditionVectorType;
    typedef std::vector<Unique<BRepOperator>> BRepOperatorPtrVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rParameters
    QuESo(const Parameters& rParameters ) : mParameters(rParameters), mMapper(mParameters) {

        mParameters.Check();
        mpTriangleMesh = MakeUnique<TriangleMesh>();
        if( mParameters.Get<bool>("embedding_flag") ) {
            // Read main mesh
            if( mParameters.Get<std::string>("input_type") == "stl_file" ){
                // Read mesh
                const auto& r_filename = mParameters.Get<std::string>("input_filename");
                IO::ReadMeshFromSTL(*mpTriangleMesh, r_filename.c_str());
                QuESo_INFO_IF(mParameters.EchoLevel() > 0) << "Read file: '" << r_filename << "'\n";
            }
            // Read conditions
            for( const auto& r_condition_settings :  mParameters.GetConditionsSettingsVector() ){
                CreateNewCondition(r_condition_settings);
            }
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
    const ElementContainerType& GetElements() const {
        return *mpElementContainer;
    }

    /// @brief Get all conditions.
    /// @return const Reference to ElementVectorPtrType
    const ConditionVectorType& GetConditions() const {
        return mConditions;
    }

    /// @brief Creates a new condition and adds to mConditions-Vector.
    /// @param rConditionParameters
    /// @return Condition&
    Condition& CreateNewCondition(const ConditionParameters& rConditionParameters);

    /// @brief Performs necessary checks.
    void Check() const;

    /// @brief Clear all containers.
    void Clear() {
        mpTriangleMesh->Clear();
        mpBRepOperator = nullptr;
        mpBrepOperatorsBC.clear();
        mpElementContainer = nullptr;
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
    Unique<ElementContainerType> mpElementContainer;
    ConditionVectorType mConditions;
    const Parameters mParameters;
    Mapper mMapper;
    ///@}

    ///@name Private Member Operations
    ///@{

    /// @brief Compute
    void Compute();
    ///@}
};

///@}

} // End namespace queso

#endif // QuESo_INCLUDE_H