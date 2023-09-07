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
#include "containers/element_container.hpp"
#include "containers/element.hpp"
#include "containers/condition.hpp"
#include "utilities/mapping_utilities.h"
#include "utilities/parameters.h"
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
    typedef std::vector<Shared<Condition>> ConditionPtrVectorType;
    typedef std::vector<Unique<BRepOperator>> BRepOperatorPtrVectorType;

    ///@}
    ///@name  Life Cycle
    ///@{

    /// @brief Constructor. Initializes Parameters and Mapper.
    QuESo(const Parameters& rParameters ) : mParameters(rParameters), mMapper(mParameters)
    {
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

    /// @brief Returns triangle mesh.
    /// @return const TriangleMesh&
    const TriangleMesh& GetTriangleMesh() {
        return mTriangleMesh;
    }


    Unique<std::vector<double>> ClosestDistances(std::vector<PointType>& rPoints,
                                                 std::vector<PointType>& rDirections) const {

        return mpBRepOperator->ClosestDistances(rPoints, rDirections);
    }

    Unique<std::vector<bool>> IsInside(std::vector<PointType>& rPoints) const {
        return mpBRepOperator->IsInside(rPoints);
    }

    ///@}
    ///@name Temporary operations to perform PosProcessing after Kratos Analysis
    ///      This will be moved to TriangleMesh.
    ///@{

    /// @brief Reads Filename and writes mesh to output/results.vtk
    /// @param Filename
    void ReadWritePostMesh() {
        QuESo_INFO << "Warning :: Postprocessing mesh is deprecated. Use VtkEmbeddedGeometryOutputProcess from Kratos. \n";
        const auto& r_filename = mParameters.Get<std::string>("postprocess_filename");
        IO::ReadMeshFromSTL(mTriangleMeshPost, r_filename.c_str());
        IO::WriteMeshToVTK(mTriangleMeshPost, "output/results.vtk", true);
    }

    /// @brief  Get mesh for prosptrocessing
    /// @return const Reference to TriangleMesh
    const TriangleMesh& GetPostMesh() const {
        QuESo_INFO << "Warning :: Postprocessing mesh is deprecated. Use VtkEmbeddedGeometryOutputProcess from Kratos. \n";
        return mTriangleMeshPost;
    }
    ///@}

private:

    ///@name Private Members Variables
    ///@{
    TriangleMesh mTriangleMesh;
    TriangleMesh mTriangleMeshPost; // Deprecated. Will be removed soon.
    Unique<BRepOperatorBase> mpBRepOperator;
    BRepOperatorPtrVectorType mpBrepOperatorsBC;
    Unique<ElementContainer> mpElementContainer;
    ConditionPtrVectorType mConditions;
    const Parameters mParameters;

    Mapper mMapper;
    ///@}

    ///@name Private Member Operations
    ///@{

    /// @brief Run QuESo
    void Compute();
    ///@}
};

///@}

} // End namespace queso

#endif // QuESo_INCLUDE_H