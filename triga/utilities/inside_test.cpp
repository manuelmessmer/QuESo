// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// Projecet includes
#include "utilities/inside_test.h"
#include "utilities/mapping_utilities.h"

bool InsideTest::is_inside_local_coordinates(PointType& point){
    PointType tmp_point = MappingUtilities::FromLocalToGlobalSpace(point, mLowerPoint, mUpperPoint);
    Point_3 point_3(tmp_point[0], tmp_point[1], tmp_point[2]);

    return InsideTest::is_inside(point_3);
}

bool InsideTest::is_inside(PointType& point){
    Point_3 point_3(point[0], point[1], point[2]);

    return InsideTest::is_inside(point_3);
}

bool InsideTest::is_inside(Point_3& point){
    // TODO: Add more test points, not only vertices of cube.
    CGAL::Bounded_side res = (*mpInsideTest)(point);
    if (res == CGAL::ON_BOUNDED_SIDE) { return true; }
    if (res == CGAL::ON_BOUNDARY) { return false; }

    return false;
}

InsideTest::IntersectionStatus InsideTest::check_intersection_status_via_element_vertices( Element& rElement ) {

    PointType lower_point = rElement.GetLocalLowerPoint();
    PointType upper_point = rElement.GetLocalUpperPoint();

    PointType point_1 = {upper_point[0], lower_point[1], lower_point[2]};
    PointType point_2 = {lower_point[0], lower_point[1], upper_point[2]};
    PointType point_3 = {lower_point[0], lower_point[1], lower_point[2]};
    PointType point_4 = {lower_point[0], upper_point[1], lower_point[2]};
    PointType point_5 = {upper_point[0], lower_point[1], upper_point[2]};
    PointType point_6 = {lower_point[0], upper_point[1], upper_point[2]};
    PointType point_7 = {upper_point[0], upper_point[1], lower_point[2]};
    PointType point_8 = {upper_point[0], upper_point[1], upper_point[2]};

    std::array<PointType,8> points = {point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8};

    int nb_inside = 0;
    for( auto point : points){
        if( is_inside_local_coordinates(point) ) { ++nb_inside; }
    }
    InsideTest::IntersectionStatus status;
    if(nb_inside == 0)
        status = Outside;
    else if(nb_inside == 8)
        status = Inside;
    else
        status = Trimmed;

    return status;
}
