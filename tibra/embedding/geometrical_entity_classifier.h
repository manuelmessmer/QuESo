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

///
/**
 * @class  GeometricalEntityClassifier
 * @author Manuel Messmer
 * @brief  Provides functions to classify geometrical entities such as element boxes and points as intersected/inside/outside.
 * @details Uses AABB Tree for fast entity classification.
 * @todo Add const.
*/
class GeometricalEntityClassifier {

public:
    ///@name Type Definitions
    ///@{
    typedef TriangleMesh::Vector3d PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ///@brief Builds AABB tree for given mesh.
    GeometricalEntityClassifier(const TriangleMesh& rTriangleMesh) : mTriangleMesh(rTriangleMesh), mTree(rTriangleMesh)
    {
    }

    ///@}
    ///@name Operations
    ///@{
    bool IsInside(const PointType& rPoints);

    bool GetIntersectionState(const Element& rElement);

    bool GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound);
    ///@}

private:
    ///@name Private Operations
    ///@{

    AABB_tree mTree;
    const TriangleMesh& mTriangleMesh;
    ///@}
};

///@}

#endif //GEOMETRICAL_ENTITY_CLASSIFICATION_INCLUDE_H