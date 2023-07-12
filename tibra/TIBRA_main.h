// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TIBRA_INCLUDE_H
#define TIBRA_INCLUDE_H

/// STL includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>

/// Project includes
#include "containers/element_container.hpp"
#include "containers/element.hpp"
#include "containers/condition.hpp"
#include "utilities/mapping_utilities.h"
#include "utilities/parameters.h"
#include "embedding/brep_operator.h"

namespace tibra {

///@name TIBRA Classes
///@{

////
/**
 * @class  TIBRA
 * @author Manuel Messmer
 * @brief  Main class of TIBRA.
*/
class TIBRA
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<Element> ElementVectorType;
    typedef std::vector<Shared<Condition>> ConditionPtrVectorType;
    typedef std::vector<Unique<BRepOperator>> BRepOperatorPtrVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor. Initializes Parameters and Mapper.
    TIBRA(const Parameters& rParameters ) : mParameters(rParameters), mMapper(mParameters)
    {
    }

    /// Copy Constructor
    TIBRA(const TIBRA &m) = delete;

    /// Copy Assignement
    TIBRA & operator= (const TIBRA &) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Run TIBRA.
    void Run();

    /// @brief Get all active elements.
    /// @return const Reference to ElementVectorPtrType
    const ElementContainer::ElementVectorPtrType& GetElements() const {
        return mpElementContainer->GetElements();
    }

    /// @brief Get all conditions.
    /// @return const Reference to ElementVectorPtrType
    const ConditionPtrVectorType& GetConditions() const {
        return mConditions;
    }

    /// @brief Clear all containers.
    void Clear() {
        mTriangleMesh.Clear();
        mTriangleMeshPost.Clear();
        mpBRepOperator = nullptr;
        mpBrepOperatorsBC.clear();
        mpElementContainer = nullptr;
        mConditions.clear();
    }

    ///@}
    ///@name Temporary operations to perform PosProcessing after Kratos Analysis
    ///      This will be moved to TriangleMesh.
    ///@{

    /// @brief Reads Filename and writes mesh to output/results.vtk
    /// @param Filename
    void ReadWritePostMesh() {
        const auto& r_filename = mParameters.Get<std::string>("postprocess_filename");
        IO::ReadMeshFromSTL(mTriangleMeshPost, r_filename.c_str());
        IO::WriteMeshToVTK(mTriangleMeshPost, "output/results.vtk", true);
    }

    /// @brief  Get mesh for prosptrocessing
    /// @return const Reference to TriangleMesh
    const TriangleMesh& GetPostMesh() const {
        return mTriangleMeshPost;
    }
    ///@}

private:

    ///@name Private Members Variables
    ///@{
    TriangleMesh mTriangleMesh;
    TriangleMesh mTriangleMeshPost;
    Unique<BRepOperatorBase> mpBRepOperator;
    BRepOperatorPtrVectorType mpBrepOperatorsBC;
    Unique<ElementContainer> mpElementContainer;
    ConditionPtrVectorType mConditions;
    const Parameters mParameters;

    Mapper mMapper;
    ///@}

    ///@name Private Member Operations
    ///@{

    /// @brief Run TIBRA
    void Compute();
    ///@}
};

///@}

} // End namespace tibra

#endif // TIBRA_INCLUDE_H