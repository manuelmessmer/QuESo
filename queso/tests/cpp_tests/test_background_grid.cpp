//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// External includes
#include <boost/test/unit_test.hpp>
// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/background_grid.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( BackgroundGridTestSuite )

typedef IntegrationPoint IntegrationPointType;
typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
typedef BackgroundGrid<ElementType> BackgroundGridType;

Unique<BackgroundGridType> CreateTestBackgroundGrid(Vector3i rNumberOfElemnts){

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_background_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_background_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElemnts);
    r_background_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-24, -43, 5});
    r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{85, 46, 115});
    r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1-0, 1.0});
    r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});

    Unique<BackgroundGridType> p_grid = MakeUnique<BackgroundGridType>(r_settings);

    IndexType number_elements = rNumberOfElemnts[0]*rNumberOfElemnts[1]*rNumberOfElemnts[2];
    for( IndexType i = 1; i <= number_elements; ++i){
        PointType tmp_point_A = {0.0, 0.0, 0.0};
        PointType tmp_point_B = {0.1, 0.1, 0.1};
        Unique<ElementType> tmp_element = MakeUnique<ElementType>(i, MakeBox(tmp_point_A, tmp_point_B),
                                                             MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}) );
        if( i != 17)
            p_grid->AddElement(std::move(tmp_element));
    }

    return p_grid;
}

BOOST_AUTO_TEST_CASE(TestBackgroundGridX) {
    QuESo_INFO << "Testing :: Test Background Grid :: Background Grid Walking along X Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto p_grid = CreateTestBackgroundGrid(number_of_elements);

    bool local_end;
    bool found;
    IndexType next_id;
    IndexType current_id = 1;
    QuESo_CHECK_EQUAL(p_grid->NumberOfActiveElements(), 23);
    IndexType active_element_counter = 1;
    for( IndexType i = 1; i < p_grid->NumberOfActiveElements() + 1; ++i){
        auto neighbour = p_grid->pGetNextElementInX(current_id, next_id, local_end);
        found = false;
        if( neighbour ){
            found = true;
            IndexType reverse_id;
            bool local_end_reversed;
            p_grid->pGetPreviousElementInX(next_id, reverse_id, local_end_reversed);
            if( local_end ) {
                QuESo_CHECK(local_end_reversed);
            }
            QuESo_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(next_id, i+1);
        if( next_id == 17 ){
            QuESo_CHECK(local_end);
            QuESo_CHECK_IS_FALSE(found);
        }
        else {
            QuESo_CHECK_EQUAL(neighbour->GetId(), i+1);
            QuESo_CHECK(found);

            if( current_id % 3 == 0 ){
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

BOOST_AUTO_TEST_CASE(TestBackgroundGridY) {
    QuESo_INFO << "Testing :: Test Background Grid :: Background Grid Walking along Y Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto p_grid = CreateTestBackgroundGrid(number_of_elements);

    bool local_end;
    bool found;
    IndexType next_id;
    IndexType current_id = 1;
    std::vector<int> test_next_ids {1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12, 13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24};
    std::vector<int> test_local_ends {10, 11, 12, 22, 23, 24};

    QuESo_CHECK_EQUAL(p_grid->NumberOfActiveElements(), 23);
    IndexType active_element_counter = 1;
    for( IndexType i = 1; i < p_grid->NumberOfActiveElements() + 1; ++i){
        found = false;
        auto neighbour = p_grid->pGetNextElementInY(current_id, next_id, local_end);
        if( neighbour ){
            found = true;
            IndexType reverse_id;
            bool local_end_reversed;
            p_grid->pGetPreviousElementInY(next_id, reverse_id, local_end_reversed);
            if( local_end ) {
                QuESo_CHECK(local_end_reversed);
            }
            QuESo_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_id), test_next_ids[i]);

        if( next_id == 17 ){
            QuESo_CHECK(local_end);
            QuESo_CHECK_IS_FALSE(found);
        }
        else {
            QuESo_CHECK_EQUAL(found, true);
            QuESo_CHECK_EQUAL(static_cast<int>(neighbour->GetId()), test_next_ids[i]);

            if( contains(test_local_ends, current_id) ){
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

BOOST_AUTO_TEST_CASE(TestBackgroundGridZ) {
    QuESo_INFO << "Testing :: Test Background Grid :: Background Grid Walking along Z Direction" << std::endl;

    Vector3i number_of_elements = {3, 4, 2};
    auto p_grid = CreateTestBackgroundGrid(number_of_elements);

    bool local_end;
    bool found;
    IndexType next_id;
    IndexType current_id = 1;
    std::vector<int> test_next_ids {1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24};
    std::vector<int> test_local_ends {13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

    QuESo_CHECK_EQUAL(p_grid->NumberOfActiveElements(), 23);
    IndexType active_element_counter = 1;
    for( IndexType i = 1; i < p_grid->NumberOfActiveElements() + 1; ++i){
        found = false;
        auto neighbour = p_grid->pGetNextElementInZ(current_id, next_id, local_end);
        if( neighbour ){
            found = true;
            IndexType reverse_id;
            bool local_end_reversed;

            p_grid->pGetPreviousElementInZ(next_id, reverse_id, local_end_reversed);
            if( local_end ) {
                QuESo_CHECK(local_end_reversed);
            }
            QuESo_CHECK_EQUAL(current_id, reverse_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_id), test_next_ids[i]);

        if( next_id == 17 ){
            QuESo_CHECK(local_end);
            QuESo_CHECK_IS_FALSE(found);
        }
        else {
            QuESo_CHECK(found);
            QuESo_CHECK_EQUAL(static_cast<int>(neighbour->GetId()), test_next_ids[i]);

            if( contains(test_local_ends, current_id) ){
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
