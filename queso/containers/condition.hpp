// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CONDITION_INCLUDE_HPP
#define CONDITION_INCLUDE_HPP

//// STL includes
#include "utilities/parameters.h"
#include "containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Condition (Base class)
 * @author Manuel Messmer
 * @brief  Interface for NeumannCondition and DirichletCondition.
**/
class Condition : public ParamCondition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    Condition(IndexType Id, const std::string& rFilename, const Vector3d& rPrescribed, Unique<TriangleMesh>& pTriangleMesh )
        : ParamCondition(Id, rFilename, rPrescribed ), mpInitialTriangleMesh(std::move(pTriangleMesh))
    {
        mConformingTriangleMesh.Reserve(10*mpInitialTriangleMesh->NumOfTriangles());
        mConformingTriangleMesh.ReserveEdgesOnPlane(10);
    }

    /// Destructor
    virtual ~Condition() = default;

    virtual const TriangleMesh& GetTriangleMesh() const {
        return *mpInitialTriangleMesh;
    }

    virtual void AddToConformingMesh(const TriangleMesh& rNewMesh){
        MeshUtilities::Append(mConformingTriangleMesh, rNewMesh);
    }

    virtual const TriangleMesh& GetConformingMesh() const {
        return mConformingTriangleMesh;
    }

private:
    ///@}
    ///@name Private Members
    ///@{

    Unique<TriangleMesh> mpInitialTriangleMesh;
    TriangleMesh mConformingTriangleMesh;

    ///@}
}; // End class ParamCondition

/**
 * @class  ConditionNeumann
 * @author Manuel Messmer
 * @brief  Container for neumann condition.
**/
class ConditionNeumann : public Condition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ConditionNeumann(IndexType Id, const std::string& rFilename, const Vector3d& rPrescribed, Unique<TriangleMesh>& pTriangleMesh )
        : Condition(Id, rFilename, rPrescribed, pTriangleMesh )
    {
    }

    ConditionNeumann(const Shared<ParamCondition>& pCondition, Unique<TriangleMesh>& pTriangleMesh )
        : Condition(pCondition->GetId(), pCondition->GetFilename(), pCondition->GetPrescribed(), pTriangleMesh )
    {
    }
    ///@}
    ///@name Operations
    ///@{

    /// Returns type of condition.
    ConditionTypeType Type() const override {
        return ConditionType::Neumann;
    }
    ///@}

}; /// End of class ConditionNeumann

/**
 * @class  ConditionDirichlet
 * @author Manuel Messmer
 * @brief  Container for dirichlet condition related parameters.
**/
class ConditionDirichlet : public Condition {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ConditionDirichlet(IndexType Id, const std::string& rFilename, const Vector3d& rPrescribed, double PenaltyFactor, Unique<TriangleMesh>& pTriangleMesh )
        : Condition(Id, rFilename, rPrescribed, pTriangleMesh ), mPenaltyFactor(PenaltyFactor)
    {
    }

    ConditionDirichlet(const Shared<ParamCondition>& pCondition, Unique<TriangleMesh>& pTriangleMesh )
        : Condition(pCondition->GetId(), pCondition->GetFilename(), pCondition->GetPrescribed(), pTriangleMesh ), mPenaltyFactor(pCondition->GetPenaltyFactor())
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// Returns type of condition.
    ConditionTypeType Type() const override {
        return ConditionType::Dirichlet;
    }

    /// Returns penalty factor.
    double GetPenaltyFactor() const override {
        return mPenaltyFactor;
    }

private:
    ///@}
    ///@name Private members
    ///@{
    double mPenaltyFactor;
    ///@}
}; // End of ConditionDirichlet class.

class ConditionFactory {
public:
    static Shared<Condition> New( const Shared<ParamCondition>& pCondition, Unique<TriangleMesh>& pTriangleMesh ){
        if( pCondition->Type() == ConditionType::Neumann ){
            return MakeShared<ConditionNeumann>(pCondition, pTriangleMesh);
        } else if ( pCondition->Type() == ConditionType::Dirichlet ) {
            return MakeShared<ConditionDirichlet>(pCondition, pTriangleMesh);
        } else {
            QuESo_ERROR << "Conditition type does not exist.\n";
        }
        return nullptr;

    }
};

} // End queso namespace.


#endif // End CONDITION_INCLUDE_HPP