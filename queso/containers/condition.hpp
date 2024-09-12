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

#ifndef CONDITION_INCLUDE_HPP
#define CONDITION_INCLUDE_HPP

//// STL includes
#include "queso/includes/settings.hpp"
#include "queso/containers/triangle_mesh.hpp"

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

    /// @brief Constructor
    /// @param pTriangleMesh Ptr is moved. Ownership is passed to Condition.
    /// @param rConditionSettings
    Condition(Unique<TriangleMeshInterface>& pTriangleMesh, const SettingsBaseType& rConditionSettings )
        : mpInitialTriangleMesh(std::move(pTriangleMesh)), mConditionSettings(rConditionSettings)
    {
        mConformingTriangleMesh.Reserve(10*mpInitialTriangleMesh->NumOfTriangles());
        mConformingTriangleMesh.ReserveEdgesOnPlane(10);
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns initial / non-conforming triangle mesh.
    /// @return const TriangleMesh&
    const TriangleMeshInterface& GetTriangleMesh() const {
        return *mpInitialTriangleMesh;
    }

    /// @brief Adds mesh section to the conforming triangle mesh.
    /// @param rNewMesh
    void AddToConformingMesh(const TriangleMeshInterface& rNewMesh){
        MeshUtilities::Append(mConformingTriangleMesh, rNewMesh);
    }

    /// @brief Returns the conforming triangle mesh.
    /// @return const TriangleMesh&
    const TriangleMeshInterface& GetConformingMesh() const {
        return mConformingTriangleMesh;
    }

    /// @brief Returns setting of condtiion
    /// @return const SettingsBaseType&
    const SettingsBaseType& GetSettings() const {
        return mConditionSettings;
    }

private:
    ///@}
    ///@name Private Members
    ///@{

    Unique<TriangleMeshInterface> mpInitialTriangleMesh;
    TriangleMesh mConformingTriangleMesh;
    const SettingsBaseType& mConditionSettings;

    ///@}
}; // End class Condition
///@}
} // End queso namespace.

#endif // End CONDITION_INCLUDE_HPP