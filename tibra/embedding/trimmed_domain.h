// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_DOMAIN_INCLUDE_H
#define TRIMMED_DOMAIN_INCLUDE_H

/// External includes
#include <memory>

/// Project includes
#include "containers/triangle_mesh.h"
#include "embedding/aabb_tree.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  TrimmedDomain
 * @author Manuel Messmer
 * @brief  Provides geometrical operations for clipped Brep models (clipped triangle meshes).
 * @details Uses AABB Tree for fast search.
*/
class TrimmedDomain {

public:
    ///@name Type Definitions
    ///@{

    enum IntersectionStatus {Inside, Outside, Trimmed};
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param pClippedTriangleMesh
    ///@note mpClippedTriangleMesh must be passed to mTree() and not pClippedTriangleMesh(), since better is moved before!
    TrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound, std::unique_ptr<TriangleMesh> pClippedTriangleMesh)
        : mpClippedTriangleMesh(std::move(pClippedTriangleMesh)), mTree(*mpClippedTriangleMesh)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    bool IsOnBoundedSided(const PointType& rPoint) const;

    ///@}

private:
    ///@name Private Operations
    ///@{

    ///@brief Returns true if point is inside AABB.
    ///@param rPoint Query point.
    ///@param rLowerBound of AABB.
    ///@param rUpperBound of AABB.
    ///@return bool
    inline bool IsContained(const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound) const{
        if(    rPoint[0] < rLowerBound[0]
            || rPoint[0] > rUpperBound[0]
            || rPoint[1] < rLowerBound[1]
            || rPoint[1] > rUpperBound[1]
            || rPoint[2] < rLowerBound[2]
            || rPoint[2] > rUpperBound[2] )
        {
            return false;
        }

        return true;
    }

    ///@}
    ///@name Private Members
    ///@{

    AABB_tree mTree;
    std::unique_ptr<TriangleMesh> mpClippedTriangleMesh;
    ///@}
};

///@}

} // End namespace tibra

#endif // TRIMMED_DOMAIN_INCLUDE_H