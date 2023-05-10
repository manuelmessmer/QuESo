// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CONDITION_INCLUDE_H
#define CONDITION_INCLUDE_H

//// STL includes

//// Project includes
#include "containers/triangle_mesh.hpp"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  Element
 * @author Manuel Messmer
 * @brief  Condition container for information required for boundary conditions.
*/
class Condition
{

public:


    typedef Unique<TriangleMesh> TriangleMeshPtrType;

    Condition(TriangleMeshPtrType pTriangleMesh, const Vector3d& rPrescribed)
        : mpTriangleMesh(std::move(pTriangleMesh)), mPrescribed(rPrescribed), mPenalty(0.0), mType(Neumann) {}

    Condition(TriangleMeshPtrType pTriangleMesh, const Vector3d& rPrescribed, double Penalty)
        : mpTriangleMesh(std::move(pTriangleMesh)), mPrescribed(rPrescribed), mPenalty(Penalty), mType(Dirichlet) {}

    ConditionValueType Type() const {
        return mType;
    }

    TriangleMesh& GetTriangleMesh() {
        return *mpTriangleMesh.get();
    }

    Vector3d GetPrescribed() const {
        return mPrescribed;
    }

    double GetPenalty() const {
        return mPenalty;
    }

private:
    Unique<TriangleMesh> mpTriangleMesh;
    Vector3d mPrescribed;
    double mPenalty;
    ConditionValueType mType;
};


} // End namespace tibra

#endif // ELEMENT_INCLUDE_H