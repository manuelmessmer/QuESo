// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INSIDE_TEST_H
#define INSIDE_TEST_H

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

// Projec includes
#include "geometries/element.h"

// External includes
#include <iostream>
#include <array>
#include <string>
#include <map>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
typedef CGAL::Side_of_triangle_mesh<SurfaceMeshType, K> InsideTestType;
typedef std::array<double,3> PointType;

class InsideTest {

public:
    enum IntersectionStatus {Inside, Outside, Trimmed};

    InsideTest(SurfaceMeshType& Polyhedron, PointType LowerPoint, PointType UpperPoint) :
        mLowerPoint(LowerPoint), mUpperPoint(UpperPoint)
    {
        // Initialize member pointers
        mpInsideTest = std::make_unique<InsideTestType>(Polyhedron);
    }

    bool is_inside_local_coordinates(PointType& point);

    bool is_inside(PointType& point);

    bool is_inside(Point_3& point);

    IntersectionStatus check_intersection_status_via_element_vertices( Element& element);

private:
    std::unique_ptr<InsideTestType> mpInsideTest;
    PointType mLowerPoint;
    PointType mUpperPoint;

};

#endif // INSIDE_TEST_H