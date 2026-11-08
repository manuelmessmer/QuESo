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
#include <doctest/doctest.h>

//// STL includes
#include <algorithm>
#include <cmath>
#include <span>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/embedding/cell_surface_section.h"
#include "queso/embedding/trimmed_domain.h"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::testing {

TEST_SUITE("BackgroundGridTestSuite")
{

    namespace {

        using IntegrationPointType = IntegrationPoint;
        using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
        using BackgroundGridType = BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>;
        using UntrimmedElementType = typename BackgroundGridType::UntrimmedElementType;
        using TrimmedElementType = typename BackgroundGridType::TrimmedElementType;
        using ElementViewType = typename BackgroundGridType::ElementViewType;
        using MainDictionaryType = typename BackgroundGridType::MainDictionaryType;
        using ElementFilter = typename BackgroundGridType::ElementFilter;
        using Direction = typename BackgroundGridType::Direction;

        [[nodiscard]] bool Contains(const std::vector<int>& rValues, int Value)
        { return std::find(rValues.begin(), rValues.end(), Value) != rValues.end(); }

        [[nodiscard]] MainDictionaryType MakeSettings(const Vector3i& rNumberOfElements)
        {
            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& r_settings = *p_settings;

            auto& r_background_grid_settings = r_settings[MainSettings::background_grid_settings];
            r_background_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
            r_background_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
            r_background_grid_settings.SetValue(
                BackgroundGridSettings::lower_bound_xyz, PointType{ -24.0, -43.0, 5.0 }
            );
            r_background_grid_settings.SetValue(
                BackgroundGridSettings::upper_bound_xyz, PointType{ 85.0, 46.0, 115.0 }
            );
            r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1.0, 1.0 });
            r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 3.0 });
            r_background_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
            r_background_grid_settings.CheckRequired();

            return std::move(r_settings);
        }

        struct NoOpBuilder
        {
            static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;
            [[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType CellIndex, const ElementBounds& rBounds)
            { return UntrimmedElementType(CellIndex + 1, rBounds); }
        };

        struct WrongIdBuilder
        {
            static constexpr ElementFilter Builds = ElementFilter::untrimmed;

            [[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType, const ElementBounds& rBounds)
            { return UntrimmedElementType(999, rBounds); }
        };

        struct CapturingBuilder
        {
            static constexpr ElementFilter Builds = ElementFilter::untrimmed;
            std::optional<ElementBounds> bounds;

            [[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType CellIndex, const ElementBounds& rBounds)
            {
                bounds = rBounds;
                return UntrimmedElementType(CellIndex + 1, rBounds);
            }
        };

        [[nodiscard]] BackgroundGridType CreateUntrimmedTraversalGrid(const Vector3i& rNumberOfElements)
        {
            auto settings = MakeSettings(rNumberOfElements);
            auto& r_settings = settings;

            BackgroundGridType grid(r_settings);
            const IndexType number_of_elements = rNumberOfElements[0] * rNumberOfElements[1] * rNumberOfElements[2];


            NoOpBuilder builder{};
            for (IndexType cell_index = 0; cell_index < number_of_elements; ++cell_index) {
                if (cell_index + 1 == 17) continue;
                grid.MakeElement(builder, cell_index);
            }

            grid.LockElements();
            return grid;
        }

        struct UntrimmedBuilder
        {
            static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;
            [[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType CellIndex, const ElementBounds& rBounds)
            { return UntrimmedElementType(CellIndex + 1, rBounds); }
        };

        [[nodiscard]] TrimmedDomain MakeClosedDomain(const BoundingBoxType& rBounds)
        {
            const PointType delta = rBounds.upper - rBounds.lower;
            const double length_scale = std::max({ delta[0], delta[1], delta[2] });
            const double coordinate_scale = std::max(
                { 1.0,
                  std::abs(rBounds.lower[0]),
                  std::abs(rBounds.lower[1]),
                  std::abs(rBounds.lower[2]),
                  std::abs(rBounds.upper[0]),
                  std::abs(rBounds.upper[1]),
                  std::abs(rBounds.upper[2]) }
            );
            const GeometryTolerance tolerance =
                GeometryTolerance::FromScale({ .length_scale = length_scale, .coordinate_scale = coordinate_scale });
            embedding::CellFaceContours contours;
            embedding::CellSurfaceSection section(
                MeshUtilities::MakeMeshBox(rBounds.lower, rBounds.upper), std::move(contours), rBounds, tolerance
            );
            return TrimmedDomain(std::move(section), [](PointView) { return false; }, 0);
        }

        struct TrimmedBuilder
        {
            static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::trimmed;
            [[nodiscard]] std::optional<TrimmedElementType> Build(IndexType CellIndex, const ElementBounds& rBounds)
            { return TrimmedElementType(CellIndex + 1, rBounds, MakeClosedDomain(rBounds.global)); }
        };

        [[nodiscard]] BackgroundGridType CreateMixedGrid(bool Lock = true)
        {
            constexpr Vector3i number_of_elements{ 2, 2, 2 };
            auto settings = MakeSettings(number_of_elements);
            auto& r_settings = settings;

            BackgroundGridType grid(r_settings);


            const std::vector<IndexType> untrimmed_ids{ 1, 3, 6 };
            const std::vector<IndexType> trimmed_ids{ 2, 5, 8 };

            UntrimmedBuilder untrimmed_builder{};
            TrimmedBuilder trimmed_builder{};

            for (IndexType element_id = 1; element_id <= 8; ++element_id) {
                if (element_id == 4 || element_id == 7) continue;
                if (std::ranges::find(untrimmed_ids, element_id) != untrimmed_ids.end()) {
                    grid.MakeElement(untrimmed_builder, element_id - 1);
                } else {
                    grid.MakeElement(trimmed_builder, element_id - 1);
                }
            }

            if (Lock) { grid.LockElements(); }
            return grid;
        }

        template<typename TRange>
        [[nodiscard]] std::vector<IndexType> CollectIds(TRange&& rRange)
        {
            std::vector<IndexType> ids{};
            for (const auto& rEntry : rRange) { ids.push_back(rEntry.GetId()); }
            return ids;
        }

    }  // namespace

    TEST_CASE("BackgroundGridTestSuite - Type Checks")
    {


        auto grid = CreateMixedGrid();
        const auto& r_const_grid = grid;

        static_assert(std::same_as<decltype(r_const_grid.GetElementView(1)), std::optional<ElementViewType>>);
        static_assert(std::same_as<decltype(grid.pGetElement<ElementFilter::untrimmed>(1)), UntrimmedElementType*>);
        static_assert(std::same_as<decltype(grid.pGetElement<ElementFilter::trimmed>(2)), TrimmedElementType*>);
        static_assert(
            std::same_as<decltype(grid.GetElements<ElementFilter::untrimmed>()), std::span<UntrimmedElementType>>
        );
        static_assert(
            std::same_as<decltype(grid.GetElements<ElementFilter::trimmed>()), std::span<TrimmedElementType>>
        );
        static_assert(std::same_as<
                      decltype(r_const_grid.GetElements<ElementFilter::untrimmed>()),
                      std::span<const UntrimmedElementType>>);
        static_assert(std::same_as<
                      decltype(r_const_grid.GetElements<ElementFilter::trimmed>()),
                      std::span<const TrimmedElementType>>);
        static_assert(std::same_as<
                      decltype(r_const_grid.GetNextElementView(1, Direction::x_forward)),
                      BackgroundGridType::NextElementViewResult>);
        static_assert(std::same_as<
                      decltype(grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::x_forward)),
                      BackgroundGridType::NextElementResult<ElementFilter::untrimmed>>);
        static_assert(std::same_as<
                      decltype(grid.GetNextElement<ElementFilter::trimmed>(1, Direction::x_forward)),
                      BackgroundGridType::NextElementResult<ElementFilter::trimmed>>);
    }

    TEST_CASE("BackgroundGridTestSuite - Access Requires Lock")
    {


        if constexpr (!NOTDEBUG) {
            auto grid = CreateMixedGrid(false);
            CHECK(!grid.ElementsAreLocked());

            CHECK_THROWS_AS((void)grid.GetElementView(1), queso::Exception);
            CHECK_THROWS_AS((void)grid.pGetElement<ElementFilter::untrimmed>(1), queso::Exception);
            CHECK_THROWS_AS((void)grid.GetElementViews(), queso::Exception);
            CHECK_THROWS_AS((void)grid.GetElements<ElementFilter::untrimmed>(), queso::Exception);
            const auto& r_const_grid = grid;
            CHECK_THROWS_AS((void)r_const_grid.GetElements<ElementFilter::trimmed>(), queso::Exception);
            grid.LockElements();
            CHECK(grid.ElementsAreLocked());
        }

        auto locked_grid = CreateMixedGrid();
        CHECK_THROWS_AS(locked_grid.ReserveElements(100, 100), queso::Exception);
        NoOpBuilder builder;
        if constexpr (!NOTDEBUG) { CHECK_THROWS_AS((void)locked_grid.MakeElement(builder, 3), queso::Exception); }
    }

    TEST_CASE("BackgroundGridTestSuite - Make Element Validates Cell Index And Builder Id")
    {
        auto settings = MakeSettings({ 2, 2, 2 });
        BackgroundGridType grid(settings);
        NoOpBuilder valid_builder;
        if constexpr (!NOTDEBUG) {
            CHECK_THROWS_AS(
                (void)grid.MakeElement(valid_builder, grid.GetGridIndexer().NumberOfElements()), queso::Exception
            );
        }

        if constexpr (!NOTDEBUG) {
            WrongIdBuilder wrong_id_builder;
            CHECK_THROWS_AS((void)grid.MakeElement(wrong_id_builder, 0), queso::Exception);
        }
    }

    TEST_CASE("BackgroundGridTestSuite - Make Element Uses Authoritative Cell Bounds")
    {
        auto settings = MakeSettings({ 2, 2, 2 });
        BackgroundGridType grid(settings);
        CapturingBuilder builder;
        constexpr IndexType cell_index = 3;
        REQUIRE(grid.MakeElement(builder, cell_index));
        REQUIRE(builder.bounds.has_value());
        const auto& r_indexer = grid.GetGridIndexer();
        const auto expected_xyz = r_indexer.GetBoundingBoxXYZFromIndex(cell_index);
        const auto expected_uvw = r_indexer.GetBoundingBoxUVWFromIndex(cell_index);
        CHECK((builder.bounds->global.lower) == (expected_xyz.lower));
        CHECK((builder.bounds->global.upper) == (expected_xyz.upper));
        CHECK((builder.bounds->parametric.lower) == (expected_uvw.lower));
        CHECK((builder.bounds->parametric.upper) == (expected_uvw.upper));
    }

    TEST_CASE("BackgroundGridTestSuite - Counts For Mixed Grid")
    {


        const auto grid = CreateMixedGrid();

        CHECK((grid.NumberOfActiveElements()) == (6UL));
        CHECK((grid.GetGridIndexer().NumberOfElements()) == (8UL));
        CHECK((grid.NumberOfUntrimmedElements()) == (3UL));
        CHECK((grid.NumberOfTrimmedElements()) == (3UL));
    }

    TEST_CASE("BackgroundGridTestSuite - Get Element View Filters")
    {


        const auto grid = CreateMixedGrid();

        const auto view_all_untrimmed = grid.GetElementView(1);
        CHECK_FALSE(!view_all_untrimmed.has_value());
        CHECK((view_all_untrimmed->GetId()) == (1UL));
        CHECK_FALSE(view_all_untrimmed->IsTrimmed());

        const auto view_all_trimmed = grid.GetElementView(2);
        CHECK_FALSE(!view_all_trimmed.has_value());
        CHECK((view_all_trimmed->GetId()) == (2UL));
        CHECK_FALSE(!view_all_trimmed->IsTrimmed());

        CHECK_FALSE(grid.GetElementView(4).has_value());

        const auto untrimmed_hit = grid.GetElementView<ElementFilter::untrimmed>(1);
        CHECK_FALSE(!untrimmed_hit.has_value());
        CHECK((untrimmed_hit->GetId()) == (1UL));
        CHECK_FALSE(untrimmed_hit->IsTrimmed());
        CHECK_FALSE(grid.GetElementView<ElementFilter::untrimmed>(2).has_value());
        CHECK_FALSE(grid.GetElementView<ElementFilter::untrimmed>(4).has_value());

        const auto trimmed_hit = grid.GetElementView<ElementFilter::trimmed>(2);
        CHECK_FALSE(!trimmed_hit.has_value());
        CHECK((trimmed_hit->GetId()) == (2UL));
        CHECK_FALSE(!trimmed_hit->IsTrimmed());
        CHECK_FALSE(grid.GetElementView<ElementFilter::trimmed>(1).has_value());
        CHECK_FALSE(grid.GetElementView<ElementFilter::trimmed>(4).has_value());
    }

    TEST_CASE("BackgroundGridTestSuite - p Get Element Filters")
    {


        auto grid = CreateMixedGrid();

        auto* p_untrimmed = grid.pGetElement<ElementFilter::untrimmed>(1);
        CHECK_FALSE(p_untrimmed == nullptr);
        CHECK((p_untrimmed->GetId()) == (1UL));
        CHECK_FALSE(p_untrimmed->IsTrimmed());

        auto* p_trimmed = grid.pGetElement<ElementFilter::trimmed>(2);
        CHECK_FALSE(p_trimmed == nullptr);
        CHECK((p_trimmed->GetId()) == (2UL));
        CHECK_FALSE(!p_trimmed->IsTrimmed());

        CHECK(grid.pGetElement<ElementFilter::untrimmed>(2) == nullptr);
        CHECK(grid.pGetElement<ElementFilter::trimmed>(1) == nullptr);
        CHECK(grid.pGetElement<ElementFilter::untrimmed>(4) == nullptr);
        CHECK(grid.pGetElement<ElementFilter::trimmed>(4) == nullptr);
    }

    TEST_CASE("BackgroundGridTestSuite - Get Element Views All And Filtered")
    {


        const auto grid = CreateMixedGrid();

        const auto ids_all = CollectIds(grid.GetElementViews());
        const auto ids_untrimmed = CollectIds(grid.GetElementViews<ElementFilter::untrimmed>());
        const auto ids_trimmed = CollectIds(grid.GetElementViews<ElementFilter::trimmed>());

        REQUIRE((ids_all.size()) == (6UL));
        REQUIRE((ids_untrimmed.size()) == (3UL));
        REQUIRE((ids_trimmed.size()) == (3UL));

        CHECK((ids_all[0]) == (1UL));
        CHECK((ids_all[1]) == (3UL));
        CHECK((ids_all[2]) == (6UL));
        CHECK((ids_all[3]) == (2UL));
        CHECK((ids_all[4]) == (5UL));
        CHECK((ids_all[5]) == (8UL));

        CHECK((ids_untrimmed[0]) == (1UL));
        CHECK((ids_untrimmed[1]) == (3UL));
        CHECK((ids_untrimmed[2]) == (6UL));

        CHECK((ids_trimmed[0]) == (2UL));
        CHECK((ids_trimmed[1]) == (5UL));
        CHECK((ids_trimmed[2]) == (8UL));

        for (const auto& rView : grid.GetElementViews<ElementFilter::untrimmed>()) { CHECK_FALSE(rView.IsTrimmed()); }
        for (const auto& rView : grid.GetElementViews<ElementFilter::trimmed>()) { CHECK_FALSE(!rView.IsTrimmed()); }
    }

    TEST_CASE("BackgroundGridTestSuite - Get Elements All Types")
    {


        auto grid = CreateMixedGrid();
        const auto& r_const_grid = grid;

        const auto untrimmed = grid.GetElements<ElementFilter::untrimmed>();
        const auto trimmed = grid.GetElements<ElementFilter::trimmed>();
        const auto const_untrimmed = r_const_grid.GetElements<ElementFilter::untrimmed>();
        const auto const_trimmed = r_const_grid.GetElements<ElementFilter::trimmed>();

        REQUIRE((untrimmed.size()) == (3UL));
        REQUIRE((trimmed.size()) == (3UL));
        REQUIRE((const_untrimmed.size()) == (3UL));
        REQUIRE((const_trimmed.size()) == (3UL));

        CHECK((untrimmed[0].GetId()) == (1UL));
        CHECK((untrimmed[1].GetId()) == (3UL));
        CHECK((untrimmed[2].GetId()) == (6UL));
        for (const auto& rElement : untrimmed) { CHECK_FALSE(rElement.IsTrimmed()); }

        CHECK((trimmed[0].GetId()) == (2UL));
        CHECK((trimmed[1].GetId()) == (5UL));
        CHECK((trimmed[2].GetId()) == (8UL));
        for (const auto& rElement : trimmed) { CHECK_FALSE(!rElement.IsTrimmed()); }
    }

    TEST_CASE("BackgroundGridTestSuite - Get Next Element Filtered X")
    {


        auto grid = CreateMixedGrid();

        const auto next_all = grid.GetNextElementView(1, Direction::x_forward);
        CHECK((next_all.next_id) == (2UL));
        CHECK_FALSE(next_all.is_end);
        CHECK_FALSE(!next_all.element.has_value());
        CHECK((next_all.element->GetId()) == (2UL));
        CHECK_FALSE(!next_all.element->IsTrimmed());

        const auto next_untrimmed_view = grid.GetNextElementView<ElementFilter::untrimmed>(1, Direction::x_forward);
        CHECK((next_untrimmed_view.next_id) == (2UL));
        CHECK(next_untrimmed_view.is_end);
        CHECK_FALSE(next_untrimmed_view.element.has_value());

        const auto next_trimmed_view = grid.GetNextElementView<ElementFilter::trimmed>(1, Direction::x_forward);
        CHECK((next_trimmed_view.next_id) == (2UL));
        CHECK_FALSE(next_trimmed_view.is_end);
        CHECK_FALSE(!next_trimmed_view.element.has_value());
        CHECK((next_trimmed_view.element->GetId()) == (2UL));

        const auto next_untrimmed = grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::x_forward);
        CHECK((next_untrimmed.next_id) == (2UL));
        CHECK(next_untrimmed.is_end);
        REQUIRE(next_untrimmed.p_element == nullptr);

        const auto next_trimmed = grid.GetNextElement<ElementFilter::trimmed>(1, Direction::x_forward);
        CHECK((next_trimmed.next_id) == (2UL));
        CHECK_FALSE(next_trimmed.is_end);
        REQUIRE(next_trimmed.p_element != nullptr);
        CHECK((next_trimmed.p_element->GetId()) == (2UL));
        CHECK_FALSE(!next_trimmed.p_element->IsTrimmed());
    }

    TEST_CASE("BackgroundGridTestSuite - Get Next Element Filtered Y")
    {


        auto grid = CreateMixedGrid();

        const auto next_all = grid.GetNextElementView(1, Direction::y_forward);
        CHECK((next_all.next_id) == (3UL));
        CHECK_FALSE(next_all.is_end);
        CHECK_FALSE(!next_all.element.has_value());
        CHECK((next_all.element->GetId()) == (3UL));
        CHECK_FALSE(next_all.element->IsTrimmed());

        const auto next_untrimmed_view = grid.GetNextElementView<ElementFilter::untrimmed>(1, Direction::y_forward);
        CHECK((next_untrimmed_view.next_id) == (3UL));
        CHECK_FALSE(next_untrimmed_view.is_end);
        CHECK_FALSE(!next_untrimmed_view.element.has_value());
        CHECK((next_untrimmed_view.element->GetId()) == (3UL));

        const auto next_trimmed_view = grid.GetNextElementView<ElementFilter::trimmed>(1, Direction::y_forward);
        CHECK((next_trimmed_view.next_id) == (3UL));
        CHECK(next_trimmed_view.is_end);
        CHECK_FALSE(next_trimmed_view.element.has_value());

        const auto next_untrimmed = grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::y_forward);
        CHECK((next_untrimmed.next_id) == (3UL));
        CHECK_FALSE(next_untrimmed.is_end);
        REQUIRE(next_untrimmed.p_element != nullptr);
        CHECK((next_untrimmed.p_element->GetId()) == (3UL));
        CHECK_FALSE(next_untrimmed.p_element->IsTrimmed());

        const auto next_trimmed = grid.GetNextElement<ElementFilter::trimmed>(1, Direction::y_forward);
        CHECK((next_trimmed.next_id) == (3UL));
        CHECK(next_trimmed.is_end);
        REQUIRE(next_trimmed.p_element == nullptr);
    }

    TEST_CASE("BackgroundGridTestSuite - Get Next Element Filtered Z")
    {


        auto grid = CreateMixedGrid();

        const auto next_all = grid.GetNextElementView(1, Direction::z_forward);
        CHECK((next_all.next_id) == (5UL));
        CHECK_FALSE(next_all.is_end);
        CHECK_FALSE(!next_all.element.has_value());
        CHECK((next_all.element->GetId()) == (5UL));
        CHECK_FALSE(!next_all.element->IsTrimmed());

        const auto next_untrimmed_view = grid.GetNextElementView<ElementFilter::untrimmed>(1, Direction::z_forward);
        CHECK((next_untrimmed_view.next_id) == (5UL));
        CHECK(next_untrimmed_view.is_end);
        CHECK_FALSE(next_untrimmed_view.element.has_value());

        const auto next_trimmed_view = grid.GetNextElementView<ElementFilter::trimmed>(1, Direction::z_forward);
        CHECK((next_trimmed_view.next_id) == (5UL));
        CHECK_FALSE(next_trimmed_view.is_end);
        CHECK_FALSE(!next_trimmed_view.element.has_value());
        CHECK((next_trimmed_view.element->GetId()) == (5UL));

        const auto next_untrimmed = grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::z_forward);
        CHECK((next_untrimmed.next_id) == (5UL));
        CHECK(next_untrimmed.is_end);
        REQUIRE(next_untrimmed.p_element == nullptr);

        const auto next_trimmed = grid.GetNextElement<ElementFilter::trimmed>(1, Direction::z_forward);
        CHECK((next_trimmed.next_id) == (5UL));
        CHECK_FALSE(next_trimmed.is_end);
        REQUIRE(next_trimmed.p_element != nullptr);
        CHECK((next_trimmed.p_element->GetId()) == (5UL));
        CHECK_FALSE(!next_trimmed.p_element->IsTrimmed());
    }

    TEST_CASE("BackgroundGridTestSuite - Get Next Element View All X")
    {


        constexpr Vector3i number_of_elements{ 3, 4, 2 };
        const auto grid = CreateUntrimmedTraversalGrid(number_of_elements);

        IndexType current_id = 1;
        CHECK((grid.NumberOfActiveElements()) == (23UL));
        IndexType active_element_counter = 1;
        for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
            const auto next_result = grid.GetNextElementView(current_id, Direction::x_forward);
            const auto element = next_result.element;
            bool found = false;
            if (element.has_value()) {
                found = true;
                const auto reverse_result = grid.GetNextElementView(next_result.next_id, Direction::x_backward);
                if (next_result.is_end) { CHECK(reverse_result.is_end); }
                CHECK((current_id) == (reverse_result.next_id));
                active_element_counter++;
            }
            CHECK((next_result.next_id) == (i + 1));
            if (next_result.next_id == 17) {
                CHECK(next_result.is_end);
                CHECK_FALSE(found);
            } else {
                CHECK((element->GetId()) == (i + 1));
                CHECK(found);

                if (current_id % 3 == 0) {
                    CHECK(next_result.is_end);
                } else {
                    CHECK_FALSE(next_result.is_end);
                }
            }
            current_id = next_result.next_id;
        }
        CHECK((active_element_counter) == (23UL));
    }

    TEST_CASE("BackgroundGridTestSuite - Get Next Element View All Y")
    {


        constexpr Vector3i number_of_elements{ 3, 4, 2 };
        const auto grid = CreateUntrimmedTraversalGrid(number_of_elements);

        IndexType current_id = 1;
        const std::vector<int> test_next_ids{ 1,  4,  7,  10, 2,  5,  8,  11, 3,  6,  9,  12,
                                              13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24 };
        const std::vector<int> test_local_ends{ 10, 11, 12, 22, 23, 24 };

        CHECK((grid.NumberOfActiveElements()) == (23UL));
        IndexType active_element_counter = 1;
        for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
            const auto next_result = grid.GetNextElementView(current_id, Direction::y_forward);
            const auto element = next_result.element;
            bool found = false;
            if (element.has_value()) {
                found = true;
                const auto reverse_result = grid.GetNextElementView(next_result.next_id, Direction::y_backward);
                if (next_result.is_end) { CHECK(reverse_result.is_end); }
                CHECK((current_id) == (reverse_result.next_id));
                active_element_counter++;
            }
            CHECK((static_cast<int>(next_result.next_id)) == (test_next_ids[i]));

            if (next_result.next_id == 17) {
                CHECK(next_result.is_end);
                CHECK_FALSE(found);
            } else {
                CHECK(found);
                CHECK((static_cast<int>(element->GetId())) == (test_next_ids[i]));

                if (Contains(test_local_ends, static_cast<int>(current_id))) {
                    CHECK(next_result.is_end);
                } else {
                    CHECK_FALSE(next_result.is_end);
                }
            }
            current_id = next_result.next_id;
        }
        CHECK((active_element_counter) == (23UL));
    }

    TEST_CASE("BackgroundGridTestSuite - Get Next Element View All Z")
    {


        constexpr Vector3i number_of_elements{ 3, 4, 2 };
        const auto grid = CreateUntrimmedTraversalGrid(number_of_elements);

        IndexType current_id = 1;
        const std::vector<int> test_next_ids{ 1, 13, 2, 14, 3, 15, 4,  16, 5,  17, 6,  18,
                                              7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24 };
        const std::vector<int> test_local_ends{ 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };

        CHECK((grid.NumberOfActiveElements()) == (23UL));
        IndexType active_element_counter = 1;
        for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
            const auto next_result = grid.GetNextElementView(current_id, Direction::z_forward);
            const auto element = next_result.element;
            bool found = false;
            if (element.has_value()) {
                found = true;
                const auto reverse_result = grid.GetNextElementView(next_result.next_id, Direction::z_backward);
                if (next_result.is_end) { CHECK(reverse_result.is_end); }
                CHECK((current_id) == (reverse_result.next_id));
                active_element_counter++;
            }
            CHECK((static_cast<int>(next_result.next_id)) == (test_next_ids[i]));

            if (next_result.next_id == 17) {
                CHECK(next_result.is_end);
                CHECK_FALSE(found);
            } else {
                CHECK(found);
                CHECK((static_cast<int>(element->GetId())) == (test_next_ids[i]));

                if (Contains(test_local_ends, static_cast<int>(current_id))) {
                    CHECK(next_result.is_end);
                } else {
                    CHECK_FALSE(next_result.is_end);
                }
            }
            current_id = next_result.next_id;
        }
        CHECK((active_element_counter) == (23UL));
    }

}  // TEST_SUITE

}  // namespace queso::testing
