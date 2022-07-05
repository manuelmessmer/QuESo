
#ifndef MODELER_INCLUDE_H
#define MODELER_INCLUDE_H

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// VTK incudes
#include <vtkCubeSource.h>
#include <vtkHexahedron.h>
#include <vtkSmartPointer.h>

namespace Modeler {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type                 Mesh;
typedef CGAL::Surface_mesh<K::Point_3>        SurfaceMeshType;

typedef std::array<double,3> PointType;


std::unique_ptr<SurfaceMeshType> make_cube_3( std::array<double,3> lower_point, std::array<double,3> upper_point);

void make_cube_3( SurfaceMeshType& rSurfaceMesh, std::array<double,3> lower_point, std::array<double,3> upper_point);

vtkSmartPointer<vtkHexahedron> GetVTKHexahedron( std::array<double,3> lower_point, std::array<double,3> upper_point);

}// End Namespace TestHelper

#endif // CUBE_MODELER_INCLUDE_H