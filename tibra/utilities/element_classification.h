// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_CLASSIFICATION_INCLUDE_H
#define ELEMENT_CLASSIFICATION_INCLUDE_H

/// External includes

/// Project includes
#include "geometries/triangle_mesh.h"

///@name TIBRA Classes
///@{

///
/**
 * @class  ElementClassification
 * @author Manuel Messmer
 * @brief  Provides functions to classify elements as inside, outisde or trimmed.
*/
class ElementClassification {

public:
    ///@name Type Definitions
    ///@{
    typedef TriangleMesh::Vector3d PointType;

    ///@}
    ///@name Operations
    ///@{
    static bool PointIsInside(const TriangleMesh& rTriangleMesh, const PointType& rPoint);
    ///@}
};

///@}

#endif //ELEMENT_CLASSIFICATION_INCLUDE_H