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


namespace EmbeddingUtilities {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
typedef Element::IntegrationPointVectorType IntegrationPointVectorType;

bool ComputeIntersectionMesh(const SurfaceMeshType& rGeometry, SurfaceMeshType& rCube,
                                    Element& rElement, const Parameters& rParam);

} // End Namespace EmbeddingUtilities

#endif // EMBEDDING_UTILITIES_INCLUDE_H