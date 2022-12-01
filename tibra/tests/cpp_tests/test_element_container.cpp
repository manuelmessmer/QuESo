// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

// Project includes
#include "utilities/parameters.h"
#include "containers/element.h"
#include "containers/element_container.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( ElementContainerTestSuite )

std::unique_ptr<ElementContainer> CreateTestElementContainer(Vector3i rNumberOfElemnts){
    PointType point_A = {0.0, 0.0, 0.0};
    PointType point_B = {1.0, 1.0, 1.0};

    Vector3i order = {2, 2, 2};
    int point_distribution_factor = 3;
    double initial_triangle_edge_length = 1;
    int minimum_number_of_triangles = 1;
    double moment_fitting_residual = 1e-8;
    std::string integration_method = "Gauss";
    int echo_level = 0;

    Parameters param(point_A, point_B, rNumberOfElemnts, order, initial_triangle_edge_length,
        minimum_number_of_triangles, moment_fitting_residual, point_distribution_factor, integration_method, echo_level);

    ElementContainer container(param);

    std::size_t number_elements = rNumberOfElemnts[0]*rNumberOfElemnts[1]*rNumberOfElemnts[2];
    for( int i = 1; i <= number_elements; ++i){
        PointType tmp_point_A = {0.0, 0.0, 0.0};
        PointType tmp_point_B = {0.1, 0.1, 0.1};
        std::shared_ptr<Element> tmp_element = std::make_shared<Element>(i, tmp_point_A, tmp_point_B, param);
        if( i != 10)
            container.AddElement(tmp_element);
    }

    return std::make_unique<ElementContainer>(std::move(container));
}

BOOST_AUTO_TEST_CASE(TestElementContainerX) {
    std::cout << "Testing :: Test Element Container :: Element Container Walking along X Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto container = CreateTestElementContainer(number_of_elements);

    bool local_end;
    bool found;
    std::size_t next_id;
    std::size_t current_id = 1;
    BOOST_CHECK_EQUAL(container->size(), 23);
    std::size_t active_element_counter = 1;
    for( int i = 1; i < container->size() + 1; ++i){
        auto neighbour = container->pGetNextElementInX(current_id, next_id, found, local_end);
        if( found ){
            std::size_t reverse_id;
            bool dummy_found;
            bool dummy_local_end;
            auto test_reverse_neighbour = container->pGetPreviousElementInX(next_id, reverse_id, dummy_found, dummy_local_end);
            BOOST_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        BOOST_CHECK_EQUAL(next_id, i+1);
        if( next_id == 10 ){
            BOOST_CHECK_EQUAL(local_end, true);
            BOOST_CHECK_EQUAL(found, false);
        }
        else {
            BOOST_CHECK_EQUAL(neighbour->GetId(), i+1);
            BOOST_CHECK_EQUAL(found, true);

            if( neighbour->GetId() % 3 == 0 ){
                BOOST_CHECK_EQUAL(local_end, true);
            }
            else {
                BOOST_CHECK_EQUAL(local_end, false);
            }
        }
        current_id = next_id;
    }
    BOOST_CHECK_EQUAL(active_element_counter, 23);
} // End Testcase

bool contains(std::vector<int>& v, int test_value){
    if(std::find(v.begin(), v.end(), test_value) != v.end()) {
        return true;
    }
    return false;
}

BOOST_AUTO_TEST_CASE(TestElementContainerY) {
    std::cout << "Testing :: Test Element Container :: Element Container Walking along Y Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto container = CreateTestElementContainer(number_of_elements);

    bool local_end;
    bool found;
    std::size_t next_id;
    std::size_t current_id = 1;
    std::vector<int> test_next_ids {1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12, 13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24};
    std::vector<int> test_local_ends {10, 11, 12, 22, 23, 24};

    BOOST_CHECK_EQUAL(container->size(), 23);
    std::size_t active_element_counter = 1;
    for( int i = 1; i < container->size() + 1; ++i){
        auto neighbour = container->pGetNextElementInY(current_id, next_id, found, local_end);
        if( found ){
            std::size_t reverse_id;
            bool dummy_found;
            bool dummy_local_end;
            auto test_reverse_neighbour = container->pGetPreviousElementInY(next_id, reverse_id, dummy_found, dummy_local_end);
            BOOST_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        BOOST_CHECK_EQUAL(next_id, test_next_ids[i]);

        if( next_id == 10 ){
            BOOST_CHECK_EQUAL(local_end, true);
            BOOST_CHECK_EQUAL(found, false);
        }
        else {
            BOOST_CHECK_EQUAL(found, true);
            BOOST_CHECK_EQUAL(neighbour->GetId(), test_next_ids[i]);

            if( contains(test_local_ends, neighbour->GetId()) ){
                BOOST_CHECK_EQUAL(local_end, true);
            }
            else {
                BOOST_CHECK_EQUAL(local_end, false);
            }
        }
        current_id = next_id;
    }
    BOOST_CHECK_EQUAL(active_element_counter, 23);
} // End Testcase

BOOST_AUTO_TEST_CASE(TestElementContainerZ) {
    std::cout << "Testing :: Test Element Container :: Element Container Walking along Z Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto container = CreateTestElementContainer(number_of_elements);

    bool local_end;
    bool found;
    std::size_t next_id;
    std::size_t current_id = 1;
    std::vector<int> test_next_ids {1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24};
    std::vector<int> test_local_ends {13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

    BOOST_CHECK_EQUAL(container->size(), 23);
    std::size_t active_element_counter = 1;
    for( int i = 1; i < container->size() + 1; ++i){
        auto neighbour = container->pGetNextElementInZ(current_id, next_id, found, local_end);
        if( found ){
            std::size_t reverse_id;
            bool dummy_found;
            bool dummy_local_end;
            auto test_reverse_neighbour = container->pGetPreviousElementInZ(next_id, reverse_id, dummy_found, dummy_local_end);
            BOOST_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        BOOST_CHECK_EQUAL(next_id, test_next_ids[i]);

        if( next_id == 10 ){
            BOOST_CHECK_EQUAL(local_end, true);
            BOOST_CHECK_EQUAL(found, false);
        }
        else {
            BOOST_CHECK_EQUAL(found, true);
            BOOST_CHECK_EQUAL(neighbour->GetId(), test_next_ids[i]);

            if( contains(test_local_ends, neighbour->GetId()) ){
                BOOST_CHECK_EQUAL(local_end, true);
            }
            else {
                BOOST_CHECK_EQUAL(local_end, false);
            }
        }
        current_id = next_id;
    }
    BOOST_CHECK_EQUAL(active_element_counter, 23);
} // End Testcase

BOOST_AUTO_TEST_SUITE_END()

} // Namespace Testing