// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CUBE_MODELER_INCLUDE_H
#define CUBE_MODELER_INCLUDE_H

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace cgal {

///@name TIBRA Classes
///@{

/**
 * @class  CuboidModeler
 * @author Manuel Messmer
*/
class CuboidModeler {

public:
    ///@name Type Definitions
    ///@{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernalType;
    typedef CGALKernalType::Point_3 CGALPointType;
    typedef CGAL::Surface_mesh<CGALPointType> CGALMeshType;
    typedef CGAL::Mesh_polyhedron_3<CGALKernalType>::type CGALPolyhedronMeshType;
    typedef std::unique_ptr<CGALMeshType> CGALMeshPtrType;
    typedef std::array<double,3> PointType;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Return cuboid CGAL mesh.
    ///@param rLowerPoint
    ///@param rUpperPoint
    ///@return CGALMeshPtrType
    static CGALMeshPtrType MakeCuboid( const PointType& rLowerPoint, const PointType& rUpperPoint);
    ///@}

}; // End class CuboidModeler

///@}

}// End cgal

#endif // CUBE_MODELER_INCLUDE_H