// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utilities/mapping_utilities.h"
#include "io/io_utilities.h"
#include "TIBRA_main.hpp"

typedef std::array<double,3> PointType;

namespace Testing{

BOOST_AUTO_TEST_SUITE( EmbeddingOperationsTestSuite )

BOOST_AUTO_TEST_CASE(Intersection) {
    std::cout << "Testing :: Test Embedding Operations :: Intersected Knot Span" << std::endl;

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

    std::string filename = "tibra/tests/cpp_tests/data/cylinder.stl";

    TIBRA tibra(filename, point_A, point_B, number_of_elements, order,
                         initial_triangle_edge_length, minimum_number_of_triangles,
                         moment_fitting_residual, point_distribution_factor, integration_method, echo_level);

    auto elements = tibra.GetElementContainer();

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

void TestElephantLarge( std::string IntegrationMethod, int p, int NumPointsInside, double Tolerance){

    std::array<double, 3> point_A = {-0.37, -0.55, -0.31};
    std::array<double, 3> point_B = {0.37, 0.55, 0.31};
    std::array<int, 3> number_of_elements = {14, 22, 12};
    std::array<int, 3> order = {p, p, p};

    int point_distribution_factor = 2;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 1000;
    double moment_fitting_residual = 1e-10;
    std::string integration_method = IntegrationMethod;
    int echo_level = 0;

    std::string filename = "tibra/tests/cpp_tests/data/elephant.stl";

    TIBRA tibra(filename, point_A, point_B, number_of_elements, order,
                         initial_triangle_edge_length, minimum_number_of_triangles,
                         moment_fitting_residual, point_distribution_factor, integration_method, echo_level);

    auto elements = tibra.GetElementContainer();

    // Compute total weight
    double weigth_trimmed = 0.0;
    double weigth_inside = 0.0;
    const auto el_it_begin = elements->begin();
    int num_elements_inside = 0;
    int num_elements_trimmed = 0;
    int num_points_inside = 0;
    for( IndexType i = 0; i < elements->size(); ++i){
        auto el_it = *(el_it_begin+i);
        if( el_it->IsTrimmed() ){
            const auto& points_trimmed = el_it->GetIntegrationPointsTrimmed();
            BOOST_CHECK_GT(points_trimmed.size(), 0);
            BOOST_CHECK_LT(points_trimmed.size(), (p+1)*(p+1)*(p+1)+1);
            for( auto point : points_trimmed ){
                weigth_trimmed += point.GetWeight();
            }
            num_elements_trimmed++;
        } else {
            const auto& points_inside = el_it->GetIntegrationPointsInside();
            //BOOST_CHECK_GT(points_inside.size(), 0);
            for( auto point : points_inside ){
                weigth_inside += point.GetWeight();
                num_points_inside++;
            }
            num_elements_inside++;
        }
    }
    BOOST_REQUIRE_EQUAL(num_elements_inside, 108);
    BOOST_CHECK_EQUAL(num_elements_trimmed, 604);
    BOOST_CHECK_EQUAL(num_points_inside, NumPointsInside);

    // Check volume inside
    const double jacobian = 0.74*1.1*0.62;
    const double ref_volume_inside = jacobian/(14*22*12)*num_elements_inside;
    const double volume_inside = weigth_inside*jacobian;
    const double rel_error_volume_inside = std::abs(volume_inside-ref_volume_inside)/ref_volume_inside;
    BOOST_CHECK_LT(rel_error_volume_inside, 1e-13);

    // Check volume inside + trimmed
    const double volume_tot = (weigth_trimmed+weigth_inside)*jacobian;
    const double ref_volume_tot = 0.0462012;
    const double rel_error_tot = std::abs(volume_tot-ref_volume_tot)/ref_volume_tot;
    BOOST_CHECK_LT(rel_error_tot, Tolerance);
}

void TestElephantSmall( std::string IntegrationMethod, int p, int NumPointsInside, double Tolerance){

    std::array<double, 3> point_A = {-0.37, -0.55, -0.31};
    std::array<double, 3> point_B = {0.37, 0.55, 0.31};
    std::array<int, 3> number_of_elements = {7, 11, 6};
    std::array<int, 3> order = {p, p, p};

    int point_distribution_factor = 2;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 1000;
    double moment_fitting_residual = 1e-10;
    std::string integration_method = IntegrationMethod;
    int echo_level = 0;

    std::string filename = "tibra/tests/cpp_tests/data/elephant.stl";

    TIBRA tibra(filename, point_A, point_B, number_of_elements, order,
                         initial_triangle_edge_length, minimum_number_of_triangles,
                         moment_fitting_residual, point_distribution_factor, integration_method, echo_level);

    auto elements = tibra.GetElementContainer();

    // Compute total weight
    double weigth_trimmed = 0.0;
    double weigth_inside = 0.0;
    const auto el_it_begin = elements->begin();
    int num_elements_inside = 0;
    int num_elements_trimmed = 0;
    int num_points_inside = 0;
    for( IndexType i = 0; i < elements->size(); ++i){
        auto el_it = *(el_it_begin+i);
        if( el_it->IsTrimmed() ){
            const auto& points_trimmed = el_it->GetIntegrationPointsTrimmed();
            BOOST_CHECK_GT(points_trimmed.size(), 0);
            BOOST_CHECK_LT(points_trimmed.size(), (p+1)*(p+1)*(p+1)+1);
            for( auto point : points_trimmed ){
                weigth_trimmed += point.GetWeight();
            }
            num_elements_trimmed++;
        } else {
            const auto& points_inside = el_it->GetIntegrationPointsInside();
            //BOOST_CHECK_GT(points_inside.size(), 0);
            for( auto point : points_inside ){
                weigth_inside += point.GetWeight();
                num_points_inside++;
            }
            num_elements_inside++;
        }
    }
    BOOST_CHECK_EQUAL(num_elements_inside, 5);
    BOOST_CHECK_EQUAL(num_elements_trimmed, 143);
    BOOST_CHECK_EQUAL(num_points_inside, NumPointsInside);

    // Check volume inside
    const double jacobian = 0.74*1.1*0.62;
    const double ref_volume_inside = jacobian/(7*11*6)*num_elements_inside;
    const double volume_inside = weigth_inside*jacobian;
    const double rel_error_volume_inside = std::abs(volume_inside-ref_volume_inside)/ref_volume_inside;
    BOOST_CHECK_LT(rel_error_volume_inside, 1e-13);

    // Check volume inside + trimmed
    const double volume_tot = (weigth_trimmed+weigth_inside)*jacobian;
    const double ref_volume_tot = 0.0462012;
    const double rel_error_tot = std::abs(volume_tot-ref_volume_tot)/ref_volume_tot;
    BOOST_CHECK_LT(rel_error_tot, Tolerance);

}

// p=2
BOOST_AUTO_TEST_CASE(VolumeElephant1) {
    std::cout << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=2)" << std::endl;
    TestElephantLarge("Gauss", 2, 2916, 0.0002);
}

BOOST_AUTO_TEST_CASE(VolumeElephant2) {
    std::cout << "Testing :: Test Embedding Operations :: Volume Elephant Optimal (p=2)" << std::endl;
    TestElephantLarge("ReducedExact", 2, 1786, 0.0002);
}

BOOST_AUTO_TEST_CASE(VolumeElephant3) {
    std::cout << "Testing :: Test Embedding Operations :: Volume Elephant ReducedOrder1 (p=2)" << std::endl;
    TestElephantLarge("ReducedOrder1", 2, 673, 0.0002);
}

BOOST_AUTO_TEST_CASE(VolumeElephant4) {
    std::cout << "Testing :: Test Embedding Operations :: Volume Elephant ReducedOrder2 (p=2)" << std::endl;
    TestElephantLarge("ReducedOrder2", 2, 406, 0.0002);
}

// p=3
BOOST_AUTO_TEST_CASE(VolumeElephant5) {
    std::cout << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=3)" << std::endl;
    TestElephantSmall("Gauss", 3, 320, 0.0002);
}
// p=4
BOOST_AUTO_TEST_CASE(VolumeElephant6) {
    std::cout << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=4)" << std::endl;
    TestElephantSmall("Gauss", 4, 625, 0.0002);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing