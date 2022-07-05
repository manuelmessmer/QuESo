// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "BaseClassModule"

#include <boost/test/unit_test.hpp>

#include "utilities/mapping_utilities.h"
#include "io/io_utilities.h"
#include "TrIGA_main.hpp"

typedef std::array<double,3> PointType;

namespace Testing{

BOOST_AUTO_TEST_CASE(SurfaceMeshPoints) {
    std::cout << "Testing :: Test Embedding Operations :: Surface Mesh Points" << std::endl;

    std::array<double, 3> point_A = {0.0, 0.0, 0.0};
    std::array<double, 3> point_B = {2.0, 2.0, 1.0};
    std::array<int, 3> number_of_elements = {1, 1, 1};
    std::array<int, 3> order = {2, 2, 2};

    int point_distribution_factor = 3;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 5000;
    double moment_fitting_residual = 1e-8;
    std::string integration_method = "Gauss";
    int echo_level = 0;

    std::string filename = "examples/data/cylinder.stl";

    TrIGA triga(filename, point_A, point_B, number_of_elements, order,
                         initial_triangle_edge_length, minimum_number_of_triangles,
                         moment_fitting_residual, point_distribution_factor, integration_method, echo_level);

    auto elements = triga.GetElements();

    BOOST_CHECK_EQUAL(elements->size(), 1);

    const auto& points_reduced = (*elements->begin())->GetIntegrationPointsTrimmed();
    BOOST_CHECK_LT(points_reduced.size(), 28);

    auto& triangles = (*elements->begin())->GetTriangles();
    BOOST_CHECK_GT(triangles.size(), minimum_number_of_triangles);

    int large_radius_count = 0;
    double area = 0.0;
    auto triangles_it_begin = triangles.begin();
    for( std::size_t i = 0; i < triangles.size(); ++i){
        auto triangle_it = triangles_it_begin + i;

        PointType coordinates = triangle_it->Center();

        BOOST_CHECK_GT(coordinates[2], -1e-6);
        BOOST_CHECK_LT(coordinates[2], 1.0+1e-6);

        BOOST_CHECK_GT(coordinates[0], -1e-6);
        BOOST_CHECK_GT(coordinates[1], -1e-6);

        if( coordinates[0] > 1e-6 && coordinates[1] > 1e-6
            && coordinates[2] > 1e-6 && coordinates[2] < 1-1e-6){ // check in z richtung fehlt noch
            double radius = std::sqrt( std::pow(coordinates[0],2) + std::pow(coordinates[1], 2) );
            BOOST_CHECK_GT(radius, 0.998);
        }
        area += triangle_it->Area();
    }
    BOOST_CHECK_LT(area, 5.141592654);
    BOOST_CHECK_GT(area, 5.135);
}

} // End namespace Testing