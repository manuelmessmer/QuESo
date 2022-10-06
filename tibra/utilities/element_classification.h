// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_CLASSIFICATION_INCLUDE_H
#define ELEMENT_CLASSIFICATION_INCLUDE_H

/// External includes
#include <memory>

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
    static std::unique_ptr<std::vector<bool>> PointsAreInside(const TriangleMesh& rTriangleMesh, const std::vector<PointType>& rPoints);
    ///@}

private:
    ///@name Private Operations
    ///@{
    static bool Anorm2(std::vector<double>& rResult, const std::vector<PointType>& rX);
    ///@}
};

///@}

#endif //ELEMENT_CLASSIFICATION_INCLUDE_H