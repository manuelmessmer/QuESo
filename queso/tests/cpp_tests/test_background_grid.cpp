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

/// External includes
#include <boost/test/unit_test.hpp>

/// Project includes
#include "queso/containers/background_grid.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/element.hpp"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso::Testing {

BOOST_AUTO_TEST_SUITE(BackgroundGridTestSuite)

using IntegrationPointType = IntegrationPoint;
using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
using ElementType = Element<IntegrationPointType, BoundaryIntegrationPointType>;
using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>;

BackgroundGridType CreateTestBackgroundGrid(const Vector3i& rNumberOfElements)
{

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_background_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_background_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
    r_background_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -24, -43, 5 });
    r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 85, 46, 115 });
    r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1.0, 1.0 });
    r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
    r_background_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
    r_background_grid_settings.CheckRequired();

    BackgroundGridType grid(r_settings);

    IndexType number_elements = rNumberOfElements[0] * rNumberOfElements[1] * rNumberOfElements[2];
    for (IndexType elements_id = 1; elements_id <= number_elements; ++elements_id) {
        ElementType new_element(
            elements_id, MakeBox({ 0.0, 0.0, 0.0 }, { 0.1, 0.1, 0.1 }), MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }));
        if (elements_id != 17) grid.AddElement(std::move(new_element));
    }
    grid.LockElements();

    return grid;
}

BOOST_AUTO_TEST_CASE(TestBackgroundGridX)
{
    QuESo_INFO << "Testing :: Test Background Grid :: Background Grid Walking along X Direction" << std::endl;

    Vector3i number_of_elements = { 3, 4, 2 };
    auto grid = CreateTestBackgroundGrid(number_of_elements);

    IndexType current_id = 1;
    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 23);
    IndexType active_element_counter = 1;
    for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
        const auto next_result = grid.GetNextElement(current_id, BackgroundGridType::Direction::x_forward);
        const auto* p_neighbour = next_result.pElement;
        bool found = false;
        if (p_neighbour) {
            found = true;
            const auto reverse_result =
                grid.GetNextElement(next_result.nextId, BackgroundGridType::Direction::x_backward);
            if (next_result.isEnd) { QuESo_CHECK(reverse_result.isEnd); }
            QuESo_CHECK_EQUAL(current_id, reverse_result.nextId);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(next_result.nextId, i + 1);
        if (next_result.nextId == 17) {
            QuESo_CHECK(next_result.isEnd);
            QuESo_CHECK_IS_FALSE(found);
        } else {
            QuESo_CHECK_EQUAL(p_neighbour->GetId(), i + 1);
            QuESo_CHECK(found);

            if (current_id % 3 == 0) {
                QuESo_CHECK(next_result.isEnd);
            } else {
                QuESo_CHECK_IS_FALSE(next_result.isEnd);
            }
        }
        current_id = next_result.nextId;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23);
}// End Testcase

bool contains(const std::vector<int>& v, int test_value)
{ return std::find(v.begin(), v.end(), test_value) != v.end(); }

BOOST_AUTO_TEST_CASE(TestBackgroundGridY)
{
    QuESo_INFO << "Testing :: Test Background Grid :: Background Grid Walking along Y Direction" << std::endl;

    Vector3i number_of_elements = { 3, 4, 2 };
    auto grid = CreateTestBackgroundGrid(number_of_elements);

    IndexType current_id = 1;
    const std::vector<int> test_next_ids{
        1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12, 13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24
    };
    const std::vector<int> test_local_ends{ 10, 11, 12, 22, 23, 24 };

    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 23);
    IndexType active_element_counter = 1;
    for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
        const auto next_result = grid.GetNextElement(current_id, BackgroundGridType::Direction::y_forward);
        const auto* p_neighbour = next_result.pElement;
        bool found = false;
        if (p_neighbour) {
            found = true;
            const auto reverse_result =
                grid.GetNextElement(next_result.nextId, BackgroundGridType::Direction::y_backward);
            if (next_result.isEnd) { QuESo_CHECK(reverse_result.isEnd); }
            QuESo_CHECK_EQUAL(current_id, reverse_result.nextId);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_result.nextId), test_next_ids[i]);

        if (next_result.nextId == 17) {
            QuESo_CHECK(next_result.isEnd);
            QuESo_CHECK_IS_FALSE(found);
        } else {
            QuESo_CHECK(found);
            QuESo_CHECK_EQUAL(static_cast<int>(p_neighbour->GetId()), test_next_ids[i]);

            if (contains(test_local_ends, current_id)) {
                QuESo_CHECK(next_result.isEnd);
            } else {
                QuESo_CHECK_IS_FALSE(next_result.isEnd);
            }
        }
        current_id = next_result.nextId;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23);
}// End Testcase

BOOST_AUTO_TEST_CASE(TestBackgroundGridZ)
{
    QuESo_INFO << "Testing :: Test Background Grid :: Background Grid Walking along Z Direction" << std::endl;

    Vector3i number_of_elements = { 3, 4, 2 };
    auto grid = CreateTestBackgroundGrid(number_of_elements);

    IndexType current_id = 1;
    const std::vector<int> test_next_ids{
        1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24
    };
    const std::vector<int> test_local_ends{ 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };

    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 23);
    IndexType active_element_counter = 1;
    for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
        const auto next_result = grid.GetNextElement(current_id, BackgroundGridType::Direction::z_forward);
        const auto* p_neighbour = next_result.pElement;
        bool found = false;
        if (p_neighbour) {
            found = true;
            const auto reverse_result =
                grid.GetNextElement(next_result.nextId, BackgroundGridType::Direction::z_backward);
            if (next_result.isEnd) { QuESo_CHECK(reverse_result.isEnd); }
            QuESo_CHECK_EQUAL(current_id, reverse_result.nextId);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_result.nextId), test_next_ids[i]);

        if (next_result.nextId == 17) {
            QuESo_CHECK(next_result.isEnd);
            QuESo_CHECK_IS_FALSE(found);
        } else {
            QuESo_CHECK(found);
            QuESo_CHECK_EQUAL(static_cast<int>(p_neighbour->GetId()), test_next_ids[i]);

            if (contains(test_local_ends, current_id)) {
                QuESo_CHECK(next_result.isEnd);
            } else {
                QuESo_CHECK_IS_FALSE(next_result.isEnd);
            }
        }
        current_id = next_result.nextId;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23);
}// End Testcase

BOOST_AUTO_TEST_SUITE_END()

}// namespace queso::Testing
