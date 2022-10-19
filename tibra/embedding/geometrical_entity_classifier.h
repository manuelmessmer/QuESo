// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef GEOMETRICAL_ENTITY_CLASSIFICATION_INCLUDE_H
#define GEOMETRICAL_ENTITY_CLASSIFICATION_INCLUDE_H

/// External includes
#include <memory>

/// Project includes
#include "geometries/triangle_mesh.h"
#include "geometries/element.h"
#include "embedding/aabb_tree.h"

///@name TIBRA Classes
///@{

/**
 * @class  GeometricalEntityClassifier
 * @author Manuel Messmer
 * @brief  Provides functions to classify geometrical entities such as element boxes and points as intersected/inside/outside.
 * @details Uses AABB Tree for fast entity classification.
*/
class GeometricalEntityClassifier {

public:
    ///@name Type Definitions
    ///@{
    typedef TriangleMesh::Vector3d PointType;

    enum IntersectionStatus {Inside, Outside, Trimmed};
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    ///@param rTriangleMesh
    GeometricalEntityClassifier(const TriangleMesh& rTriangleMesh) : mTriangleMesh(rTriangleMesh), mTree(rTriangleMesh)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TriangleMesh.
    ///@param rPoint
    ///@return bool
    bool IsInside(const PointType& rPoint) const;

    ///@brief Returns intersections state of element.
    ///@param rElement
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const Element& rElement) const;

    ///@brief Returns intersections state of element.
    ///@param rLowerBound
    ///@param rUpperBound
    ///@param Tolerance Tolerance reduces element slightly. If Tolerance=0 touch is detected as intersection.
    ///                 If Tolerance>0, touch is not detected as intersection.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatus GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance) const;
    ///@}

private:
    ///@name Private Operations
    ///@{

    AABB_tree mTree;
    const TriangleMesh& mTriangleMesh;
    ///@}
};

///@}

#endif // GEOMETRICAL_ENTITY_CLASSIFICATION_INCLUDE_H