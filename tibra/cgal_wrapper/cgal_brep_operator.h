// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef EMBEDDING_UTILITIES_INCLUDE_H
#define EMBEDDING_UTILITIES_INCLUDE_H

/// CGAL includes
// Domain
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
// Project includes
#include "geometries/element.h"
#include "geometries/integration_point.h"
#include "utilities/parameters.h"

namespace cgal {

///@name TIBRA Classes
///@{

/**
 * @class  CGAL brep operator
 * @author Manuel Messmer
 * @brief Provides geometrical operations for Brep models.
*/
class BRepOperator{

public:
    ///@name Type Definitions
    ///@{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
    typedef Element::IntegrationPointVectorType IntegrationPointVectorType;

    ///@}
    ///@name Operatios
    ///@{

    ///@brief Computes intersection mesh between rGeometry and rCube.
    ///@param rGeometry
    ///@param rCube
    ///@param rElement
    ///@param rParam
    ///@return bool Success?
    static bool ComputeIntersectionMesh(const SurfaceMeshType& rGeometry, SurfaceMeshType& rCube,
                                        Element& rElement, const Parameters& rParam);

}; // End Class BRepOperator

} // End namespace cgal

#endif // EMBEDDING_UTILITIES_INCLUDE_H