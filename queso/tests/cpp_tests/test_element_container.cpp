//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/parameters.h"
#include "queso/containers/element.hpp"
#include "queso/containers/element_container.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ElementContainerTestSuite )

typedef IntegrationPoint IntegrationPointType;
typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
typedef ElementContainer<ElementType> ElementContainerType;

Unique<ElementContainerType> CreateTestElementContainer(Vector3i rNumberOfElemnts){

    Parameters param( {Component("number_of_elements", rNumberOfElemnts) } );
    Unique<ElementContainerType> container = MakeUnique<ElementContainerType>(param);

    IndexType number_elements = rNumberOfElemnts[0]*rNumberOfElemnts[1]*rNumberOfElemnts[2];
    for( IndexType i = 1; i <= number_elements; ++i){
        PointType tmp_point_A = {0.0, 0.0, 0.0};
        PointType tmp_point_B = {0.1, 0.1, 0.1};
        Unique<ElementType> tmp_element = MakeUnique<ElementType>(i, MakeBox(tmp_point_A, tmp_point_B),
                                                             MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}) );
        if( i != 10)
            container->AddElement(tmp_element);
    }

    return container;
}

BOOST_AUTO_TEST_CASE(TestElementContainerX) {
    QuESo_INFO << "Testing :: Test Element Container :: Element Container Walking along X Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto container = CreateTestElementContainer(number_of_elements);

    bool local_end;
    bool found;
    IndexType next_id;
    IndexType current_id = 1;
    QuESo_CHECK_EQUAL(container->size(), 23);
    IndexType active_element_counter = 1;
    for( IndexType i = 1; i < container->size() + 1; ++i){
        auto neighbour = container->pGetNextElementInX(current_id, next_id, found, local_end);
        if( found ){
            IndexType reverse_id;
            bool dummy_found;
            bool dummy_local_end;
            container->pGetPreviousElementInX(next_id, reverse_id, dummy_found, dummy_local_end);
            QuESo_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(next_id, i+1);
        if( next_id == 10 ){
            QuESo_CHECK(local_end);
            QuESo_CHECK_IS_FALSE(found);
        }
        else {
            QuESo_CHECK_EQUAL(neighbour->GetId(), i+1);
            QuESo_CHECK(found);

            if( neighbour->GetId() % 3 == 0 ){
                QuESo_CHECK(local_end);
            }
            else {
                QuESo_CHECK_IS_FALSE(local_end);
            }
        }
        current_id = next_id;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23);
} // End Testcase

bool contains(std::vector<int>& v, int test_value){
    if(std::find(v.begin(), v.end(), test_value) != v.end()) {
        return true;
    }
    return false;
}

BOOST_AUTO_TEST_CASE(TestElementContainerY) {
    QuESo_INFO << "Testing :: Test Element Container :: Element Container Walking along Y Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto container = CreateTestElementContainer(number_of_elements);

    bool local_end;
    bool found;
    IndexType next_id;
    IndexType current_id = 1;
    std::vector<int> test_next_ids {1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12, 13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24};
    std::vector<int> test_local_ends {10, 11, 12, 22, 23, 24};

    QuESo_CHECK_EQUAL(container->size(), 23);
    IndexType active_element_counter = 1;
    for( IndexType i = 1; i < container->size() + 1; ++i){
        auto neighbour = container->pGetNextElementInY(current_id, next_id, found, local_end);
        if( found ){
            IndexType reverse_id;
            bool dummy_found;
            bool dummy_local_end;
            container->pGetPreviousElementInY(next_id, reverse_id, dummy_found, dummy_local_end);
            QuESo_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_id), test_next_ids[i]);

        if( next_id == 10 ){
            QuESo_CHECK(local_end);
            QuESo_CHECK_IS_FALSE(found);
        }
        else {
            QuESo_CHECK_EQUAL(found, true);
            QuESo_CHECK_EQUAL(static_cast<int>(neighbour->GetId()), test_next_ids[i]);

            if( contains(test_local_ends, neighbour->GetId()) ){
                QuESo_CHECK(local_end);
            }
            else {
                QuESo_CHECK_IS_FALSE(local_end);
            }
        }
        current_id = next_id;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23);
} // End Testcase

BOOST_AUTO_TEST_CASE(TestElementContainerZ) {
    QuESo_INFO << "Testing :: Test Element Container :: Element Container Walking along Z Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto container = CreateTestElementContainer(number_of_elements);

    bool local_end;
    bool found;
    IndexType next_id;
    IndexType current_id = 1;
    std::vector<int> test_next_ids {1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24};
    std::vector<int> test_local_ends {13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

    QuESo_CHECK_EQUAL(container->size(), 23);
    IndexType active_element_counter = 1;
    for( IndexType i = 1; i < container->size() + 1; ++i){
        auto neighbour = container->pGetNextElementInZ(current_id, next_id, found, local_end);
        if( found ){
            IndexType reverse_id;
            bool dummy_found;
            bool dummy_local_end;
            container->pGetPreviousElementInZ(next_id, reverse_id, dummy_found, dummy_local_end);
            QuESo_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_id), test_next_ids[i]);

        if( next_id == 10 ){
            QuESo_CHECK(local_end);
            QuESo_CHECK_IS_FALSE(found);
        }
        else {
            QuESo_CHECK(found);
            QuESo_CHECK_EQUAL(static_cast<int>(neighbour->GetId()), test_next_ids[i]);

            if( contains(test_local_ends, neighbour->GetId()) ){
                QuESo_CHECK(local_end);
            }
            else {
                QuESo_CHECK_IS_FALSE(local_end);
            }
        }
        current_id = next_id;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23);
} // End Testcase

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso