// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CONDITION_INCLUDE_HPP
#define CONDITION_INCLUDE_HPP

//// STL includes
#include "includes/parameters.h"
#include "containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Condition
 * @author Manuel Messmer
 * @brief  Interface for conditions. Stores condition settings and triangle mesh.
**/
class Condition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Condition(Unique<TriangleMesh>& pTriangleMesh, const ConditionParameters& rConditionParameters )
        : mpInitialTriangleMesh(std::move(pTriangleMesh)), mConditionParameters(rConditionParameters)
    {
        mConformingTriangleMesh.Reserve(10*mpInitialTriangleMesh->NumOfTriangles());
        mConformingTriangleMesh.ReserveEdgesOnPlane(10);
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns initial / non-conforming triangle mesh.
    /// @return const TriangleMesh&
    const TriangleMesh& GetTriangleMesh() const {
        return *mpInitialTriangleMesh;
    }

    /// @brief Adds mesh section to the conforming triangle mesh.
    /// @param rNewMesh
    void AddToConformingMesh(const TriangleMesh& rNewMesh){
        MeshUtilities::Append(mConformingTriangleMesh, rNewMesh);
    }

    /// @brief Returns the conforming triangle mesh.
    /// @return
    const TriangleMesh& GetConformingMesh() const {
        return mConformingTriangleMesh;
    }

    /// @brief Returns setting of condtiion
    /// @return const ConditionParameters&
    const ConditionParameters& GetSettings() const {
        return mConditionParameters;
    }

private:
    ///@}
    ///@name Private Members
    ///@{

    Unique<TriangleMesh> mpInitialTriangleMesh;
    TriangleMesh mConformingTriangleMesh;
    const ConditionParameters& mConditionParameters;

    ///@}
}; // End class Condition

} // End queso namespace.


#endif // End CONDITION_INCLUDE_HPP