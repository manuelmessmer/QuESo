// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INSIDE_TEST_H
#define INSIDE_TEST_H

//VTK includes
#include <vtkSelectEnclosedPoints.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

// Projec includes
#include "geometries/element.h"

// External includes
#include <iostream>
#include <array>
#include <string>
#include <map>

class Element;

class InsideTest {

public:
    typedef std::array<double,3> PointType;

    enum IntersectionStatus {Inside, Outside, Trimmed};

    InsideTest(vtkSmartPointer<vtkPolyData>& Polyhedron, PointType LowerPoint, PointType UpperPoint) :
        mLowerPoint(LowerPoint), mUpperPoint(UpperPoint)
    {
        // Initialize member pointers
        mpInsideTest = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
        mpInsideTest->SetTolerance(1e-6);
        mpInsideTest->Initialize(Polyhedron);
    }

    ~InsideTest(){
        mpInsideTest->Complete();
    }

    bool is_inside_local_coordinates(PointType& point);

    bool is_inside(PointType& point);

    IntersectionStatus check_intersection_status_via_element_vertices( Element& element);

private:
    vtkSmartPointer<vtkSelectEnclosedPoints> mpInsideTest;
    PointType mLowerPoint;
    PointType mUpperPoint;

};

#endif // INSIDE_TEST_H