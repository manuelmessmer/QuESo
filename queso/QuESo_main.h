// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef QuESo_INCLUDE_H
#define QuESo_INCLUDE_H

/// STL includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>

/// Project includes
#include "includes/define.hpp"
#include "containers/element_container.hpp"
#include "containers/element.hpp"
#include "containers/condition.hpp"
#include "utilities/mapping_utilities.h"
#include "includes/parameters.h"
#include "embedding/brep_operator.h"

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

    typedef std::vector<Element> ElementVectorType;
    typedef std::vector<Condition> ConditionVectorType;
    typedef std::vector<Unique<BRepOperator>> BRepOperatorPtrVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rParameters
    QuESo(const Parameters& rParameters ) : mParameters(rParameters), mMapper(mParameters) {
        if( mParameters.Get<bool>("embedding_flag") ) {
            // Read main mesh
            if( mParameters.Get<std::string>("input_type") == "stl_file" ){
                // Read mesh
                const auto& r_filename = mParameters.Get<std::string>("input_filename");
                IO::ReadMeshFromSTL(mTriangleMesh, r_filename.c_str());
                QuESo_INFO_IF(mParameters.EchoLevel() > 0) << "Read file: '" << r_filename << "'\n";
            }
            // Read conditions
            for( const auto& r_condition_settings :  mParameters.GetConditionsSettingsVector() ){
                CreateNewCondition(r_condition_settings);
            }
        }
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
    const ElementContainer::ElementVectorPtrType& GetElements() const {
        return mpElementContainer->GetElements();
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

    /// @brief Clear all containers.
    void Clear() {
        mTriangleMesh.Clear();
        mpBRepOperator = nullptr;
        mpBrepOperatorsBC.clear();
        mpElementContainer = nullptr;
        mConditions.clear();
    }

    /// @brief Returns triangle mesh.
    /// @return const TriangleMesh&
    const TriangleMesh& GetTriangleMesh() {
        return mTriangleMesh;
    }

    ///@}

private:

    ///@name Private Members Variables
    ///@{
    TriangleMesh mTriangleMesh;
    Unique<BRepOperator> mpBRepOperator;
    BRepOperatorPtrVectorType mpBrepOperatorsBC;
    Unique<ElementContainer> mpElementContainer;
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