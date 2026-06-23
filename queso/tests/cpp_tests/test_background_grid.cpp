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

//// STL includes
#include <algorithm>
#include <span>

//// Project includes
#include "queso/containers/background_grid.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/clipped_triangle_mesh.hpp"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso::Testing {

BOOST_AUTO_TEST_SUITE(BackgroundGridTestSuite)

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
        r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ -24.0, -43.0, 5.0 });
        r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 85.0, 46.0, 115.0 });
        r_background_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ -1.0, -1.0, 1.0 });
        r_background_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 3.0 });
        r_background_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_background_grid_settings.CheckRequired();

        return std::move(r_settings);
    }

	struct NoOpBuilder
	{
		static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;
		[[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType Id, const ElementBounds& rBounds)
		{ return UntrimmedElementType(Id, rBounds); }
	};

    [[nodiscard]] BackgroundGridType CreateUntrimmedTraversalGrid(const Vector3i& rNumberOfElements)
    {
        auto settings = MakeSettings(rNumberOfElements);
        auto& r_settings = settings;

        BackgroundGridType grid(r_settings);
        GridIndexer grid_indexer(r_settings);
        const IndexType number_of_elements = rNumberOfElements[0] * rNumberOfElements[1] * rNumberOfElements[2];


        NoOpBuilder builder{};
        for (IndexType element_id = 1; element_id <= number_of_elements; ++element_id) {
            if (element_id == 17) continue;
            const auto bounds_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(element_id - 1);
            const auto bounds_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(element_id - 1);
            grid.MakeElement(builder, element_id, ElementBounds{ bounds_xyz, bounds_uvw });
        }

        grid.LockElements();
        return grid;
    }

	struct UntrimmedBuilder
	{
		static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;
		[[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType Id, const ElementBounds& rBounds)
		{ return UntrimmedElementType(Id, rBounds); }
	};

	struct TrimmedBuilder
	{
		static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::trimmed;
		[[nodiscard]] std::optional<TrimmedElementType> Build(IndexType Id, const ElementBounds& rBounds)
		{
			ClippedTriangleMesh clipped_mesh{};
			return TrimmedElementType(
					Id,
					rBounds,
					TrimmedDomain(std::move(clipped_mesh), rBounds.global.lower, rBounds.global.upper, nullptr, 0)
					);
		}
	};

    [[nodiscard]] BackgroundGridType CreateMixedGrid(bool Lock = true)
    {
        constexpr Vector3i number_of_elements{ 2, 2, 2 };
        auto settings = MakeSettings(number_of_elements);
        auto& r_settings = settings;

        BackgroundGridType grid(r_settings);
        GridIndexer grid_indexer(r_settings);


        const std::vector<IndexType> untrimmed_ids{ 1, 3, 6 };
        const std::vector<IndexType> trimmed_ids{ 2, 5, 8 };

        UntrimmedBuilder untrimmed_builder{};
        TrimmedBuilder trimmed_builder{};

        for (IndexType element_id = 1; element_id <= 8; ++element_id) {
            if (element_id == 4 || element_id == 7) continue;
            const auto bounds_xyz = grid_indexer.GetBoundingBoxXYZFromIndex(element_id - 1);
            const auto bounds_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(element_id - 1);
            const ElementBounds bounds{ bounds_xyz, bounds_uvw };

            if (std::ranges::find(untrimmed_ids, element_id) != untrimmed_ids.end()) {
                grid.MakeElement(untrimmed_builder, element_id, bounds);
            } else {
                grid.MakeElement(trimmed_builder, element_id, bounds);
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

}// namespace

BOOST_AUTO_TEST_CASE(TypeChecks)
{
    QuESo_INFO << "Testing :: Background Grid :: TypeChecks" << std::endl;

    auto grid = CreateMixedGrid();
    const auto& r_const_grid = grid;

    static_assert(std::same_as<decltype(r_const_grid.GetElementView(1)), std::optional<ElementViewType>>);
    static_assert(std::same_as<decltype(grid.pGetElement<ElementFilter::untrimmed>(1)), UntrimmedElementType*>);
    static_assert(std::same_as<decltype(grid.pGetElement<ElementFilter::trimmed>(2)), TrimmedElementType*>);
    static_assert(
        std::same_as<decltype(grid.GetElements<ElementFilter::untrimmed>()), std::span<UntrimmedElementType>>
    );
    static_assert(std::same_as<decltype(grid.GetElements<ElementFilter::trimmed>()), std::span<TrimmedElementType>>);
    static_assert(std::same_as<
                  decltype(r_const_grid.GetElements<ElementFilter::untrimmed>()),
                  std::span<const UntrimmedElementType>>);
    static_assert(
        std::same_as<decltype(r_const_grid.GetElements<ElementFilter::trimmed>()), std::span<const TrimmedElementType>>
    );
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

BOOST_AUTO_TEST_CASE(AccessRequiresLock)
{
    QuESo_INFO << "Testing :: Background Grid :: AccessRequiresLock" << std::endl;

    if constexpr (!NOTDEBUG) {
        auto grid = CreateMixedGrid(false);

        BOOST_CHECK_THROW((void)grid.GetElementView(1), queso::Exception);
        BOOST_CHECK_THROW((void)grid.pGetElement<ElementFilter::untrimmed>(1), queso::Exception);
    }
}

BOOST_AUTO_TEST_CASE(CountsForMixedGrid)
{
    QuESo_INFO << "Testing :: Background Grid :: CountsForMixedGrid" << std::endl;

    const auto grid = CreateMixedGrid();

    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 6UL);
    QuESo_CHECK_EQUAL(grid.NumberOfUntrimmedElements(), 3UL);
    QuESo_CHECK_EQUAL(grid.NumberOfTrimmedElements(), 3UL);
}

BOOST_AUTO_TEST_CASE(GetElementViewFilters)
{
    QuESo_INFO << "Testing :: Background Grid :: GetElementViewFilters" << std::endl;

    const auto grid = CreateMixedGrid();

    const auto view_all_untrimmed = grid.GetElementView(1);
    QuESo_CHECK_IS_FALSE(!view_all_untrimmed.has_value());
    QuESo_CHECK_EQUAL(view_all_untrimmed->GetId(), 1UL);
    QuESo_CHECK_IS_FALSE(view_all_untrimmed->IsTrimmed());

    const auto view_all_trimmed = grid.GetElementView(2);
    QuESo_CHECK_IS_FALSE(!view_all_trimmed.has_value());
    QuESo_CHECK_EQUAL(view_all_trimmed->GetId(), 2UL);
    QuESo_CHECK_IS_FALSE(!view_all_trimmed->IsTrimmed());

    QuESo_CHECK_IS_FALSE(grid.GetElementView(4).has_value());

    const auto untrimmed_hit = grid.GetElementView<ElementFilter::untrimmed>(1);
    QuESo_CHECK_IS_FALSE(!untrimmed_hit.has_value());
    QuESo_CHECK_EQUAL(untrimmed_hit->GetId(), 1UL);
    QuESo_CHECK_IS_FALSE(untrimmed_hit->IsTrimmed());
    QuESo_CHECK_IS_FALSE(grid.GetElementView<ElementFilter::untrimmed>(2).has_value());
    QuESo_CHECK_IS_FALSE(grid.GetElementView<ElementFilter::untrimmed>(4).has_value());

    const auto trimmed_hit = grid.GetElementView<ElementFilter::trimmed>(2);
    QuESo_CHECK_IS_FALSE(!trimmed_hit.has_value());
    QuESo_CHECK_EQUAL(trimmed_hit->GetId(), 2UL);
    QuESo_CHECK_IS_FALSE(!trimmed_hit->IsTrimmed());
    QuESo_CHECK_IS_FALSE(grid.GetElementView<ElementFilter::trimmed>(1).has_value());
    QuESo_CHECK_IS_FALSE(grid.GetElementView<ElementFilter::trimmed>(4).has_value());
}

BOOST_AUTO_TEST_CASE(pGetElementFilters)
{
    QuESo_INFO << "Testing :: Background Grid :: pGetElementFilters" << std::endl;

    auto grid = CreateMixedGrid();

    auto* p_untrimmed = grid.pGetElement<ElementFilter::untrimmed>(1);
    QuESo_CHECK_IS_FALSE(p_untrimmed == nullptr);
    QuESo_CHECK_EQUAL(p_untrimmed->GetId(), 1UL);
    QuESo_CHECK_IS_FALSE(p_untrimmed->IsTrimmed());

    auto* p_trimmed = grid.pGetElement<ElementFilter::trimmed>(2);
    QuESo_CHECK_IS_FALSE(p_trimmed == nullptr);
    QuESo_CHECK_EQUAL(p_trimmed->GetId(), 2UL);
    QuESo_CHECK_IS_FALSE(!p_trimmed->IsTrimmed());

    QuESo_CHECK(grid.pGetElement<ElementFilter::untrimmed>(2) == nullptr);
    QuESo_CHECK(grid.pGetElement<ElementFilter::trimmed>(1) == nullptr);
    QuESo_CHECK(grid.pGetElement<ElementFilter::untrimmed>(4) == nullptr);
    QuESo_CHECK(grid.pGetElement<ElementFilter::trimmed>(4) == nullptr);
}

BOOST_AUTO_TEST_CASE(GetElementViewsAllAndFiltered)
{
    QuESo_INFO << "Testing :: Background Grid :: GetElementViewsAllAndFiltered" << std::endl;

    const auto grid = CreateMixedGrid();

    const auto ids_all = CollectIds(grid.GetElementViews());
    const auto ids_untrimmed = CollectIds(grid.GetElementViews<ElementFilter::untrimmed>());
    const auto ids_trimmed = CollectIds(grid.GetElementViews<ElementFilter::trimmed>());

    QuESo_CHECK_EQUAL(ids_all.size(), 6UL);
    QuESo_CHECK_EQUAL(ids_untrimmed.size(), 3UL);
    QuESo_CHECK_EQUAL(ids_trimmed.size(), 3UL);

    QuESo_CHECK_EQUAL(ids_all[0], 1UL);
    QuESo_CHECK_EQUAL(ids_all[1], 3UL);
    QuESo_CHECK_EQUAL(ids_all[2], 6UL);
    QuESo_CHECK_EQUAL(ids_all[3], 2UL);
    QuESo_CHECK_EQUAL(ids_all[4], 5UL);
    QuESo_CHECK_EQUAL(ids_all[5], 8UL);

    QuESo_CHECK_EQUAL(ids_untrimmed[0], 1UL);
    QuESo_CHECK_EQUAL(ids_untrimmed[1], 3UL);
    QuESo_CHECK_EQUAL(ids_untrimmed[2], 6UL);

    QuESo_CHECK_EQUAL(ids_trimmed[0], 2UL);
    QuESo_CHECK_EQUAL(ids_trimmed[1], 5UL);
    QuESo_CHECK_EQUAL(ids_trimmed[2], 8UL);

    for (const auto& rView : grid.GetElementViews<ElementFilter::untrimmed>()) {
        QuESo_CHECK_IS_FALSE(rView.IsTrimmed());
    }
    for (const auto& rView : grid.GetElementViews<ElementFilter::trimmed>()) {
        QuESo_CHECK_IS_FALSE(!rView.IsTrimmed());
    }
}

BOOST_AUTO_TEST_CASE(GetElementsAllTypes)
{
    QuESo_INFO << "Testing :: Background Grid :: GetElementsAllTypes" << std::endl;

    auto grid = CreateMixedGrid();
    const auto& r_const_grid = grid;

    const auto untrimmed = grid.GetElements<ElementFilter::untrimmed>();
    const auto trimmed = grid.GetElements<ElementFilter::trimmed>();
    const auto const_untrimmed = r_const_grid.GetElements<ElementFilter::untrimmed>();
    const auto const_trimmed = r_const_grid.GetElements<ElementFilter::trimmed>();

    QuESo_CHECK_EQUAL(untrimmed.size(), 3UL);
    QuESo_CHECK_EQUAL(trimmed.size(), 3UL);
    QuESo_CHECK_EQUAL(const_untrimmed.size(), 3UL);
    QuESo_CHECK_EQUAL(const_trimmed.size(), 3UL);

    QuESo_CHECK_EQUAL(untrimmed[0].GetId(), 1UL);
    QuESo_CHECK_EQUAL(untrimmed[1].GetId(), 3UL);
    QuESo_CHECK_EQUAL(untrimmed[2].GetId(), 6UL);
    for (const auto& rElement : untrimmed) { QuESo_CHECK_IS_FALSE(rElement.IsTrimmed()); }

    QuESo_CHECK_EQUAL(trimmed[0].GetId(), 2UL);
    QuESo_CHECK_EQUAL(trimmed[1].GetId(), 5UL);
    QuESo_CHECK_EQUAL(trimmed[2].GetId(), 8UL);
    for (const auto& rElement : trimmed) { QuESo_CHECK_IS_FALSE(!rElement.IsTrimmed()); }
}

BOOST_AUTO_TEST_CASE(GetNextElementFilteredX)
{
    QuESo_INFO << "Testing :: Background Grid :: GetNextElementFilteredX" << std::endl;

    auto grid = CreateMixedGrid();

    const auto next_all = grid.GetNextElementView(1, Direction::x_forward);
    QuESo_CHECK_EQUAL(next_all.next_id, 2UL);
    QuESo_CHECK_IS_FALSE(next_all.is_end);
    QuESo_CHECK_IS_FALSE(!next_all.element.has_value());
    QuESo_CHECK_EQUAL(next_all.element->GetId(), 2UL);
    QuESo_CHECK_IS_FALSE(!next_all.element->IsTrimmed());

    const auto next_untrimmed_view = grid.GetNextElementView<ElementFilter::untrimmed>(1, Direction::x_forward);
    QuESo_CHECK_EQUAL(next_untrimmed_view.next_id, 2UL);
    QuESo_CHECK(next_untrimmed_view.is_end);
    QuESo_CHECK_IS_FALSE(next_untrimmed_view.element.has_value());

    const auto next_trimmed_view = grid.GetNextElementView<ElementFilter::trimmed>(1, Direction::x_forward);
    QuESo_CHECK_EQUAL(next_trimmed_view.next_id, 2UL);
    QuESo_CHECK_IS_FALSE(next_trimmed_view.is_end);
    QuESo_CHECK_IS_FALSE(!next_trimmed_view.element.has_value());
    QuESo_CHECK_EQUAL(next_trimmed_view.element->GetId(), 2UL);

    const auto next_untrimmed = grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::x_forward);
    QuESo_CHECK_EQUAL(next_untrimmed.next_id, 2UL);
    QuESo_CHECK(next_untrimmed.is_end);
    QuESo_CHECK(next_untrimmed.p_element == nullptr);

    const auto next_trimmed = grid.GetNextElement<ElementFilter::trimmed>(1, Direction::x_forward);
    QuESo_CHECK_EQUAL(next_trimmed.next_id, 2UL);
    QuESo_CHECK_IS_FALSE(next_trimmed.is_end);
    QuESo_CHECK(next_trimmed.p_element != nullptr);
    QuESo_CHECK_EQUAL(next_trimmed.p_element->GetId(), 2UL);
    QuESo_CHECK_IS_FALSE(!next_trimmed.p_element->IsTrimmed());
}

BOOST_AUTO_TEST_CASE(GetNextElementFilteredY)
{
    QuESo_INFO << "Testing :: Background Grid :: GetNextElementFilteredY" << std::endl;

    auto grid = CreateMixedGrid();

    const auto next_all = grid.GetNextElementView(1, Direction::y_forward);
    QuESo_CHECK_EQUAL(next_all.next_id, 3UL);
    QuESo_CHECK_IS_FALSE(next_all.is_end);
    QuESo_CHECK_IS_FALSE(!next_all.element.has_value());
    QuESo_CHECK_EQUAL(next_all.element->GetId(), 3UL);
    QuESo_CHECK_IS_FALSE(next_all.element->IsTrimmed());

    const auto next_untrimmed_view = grid.GetNextElementView<ElementFilter::untrimmed>(1, Direction::y_forward);
    QuESo_CHECK_EQUAL(next_untrimmed_view.next_id, 3UL);
    QuESo_CHECK_IS_FALSE(next_untrimmed_view.is_end);
    QuESo_CHECK_IS_FALSE(!next_untrimmed_view.element.has_value());
    QuESo_CHECK_EQUAL(next_untrimmed_view.element->GetId(), 3UL);

    const auto next_trimmed_view = grid.GetNextElementView<ElementFilter::trimmed>(1, Direction::y_forward);
    QuESo_CHECK_EQUAL(next_trimmed_view.next_id, 3UL);
    QuESo_CHECK(next_trimmed_view.is_end);
    QuESo_CHECK_IS_FALSE(next_trimmed_view.element.has_value());

    const auto next_untrimmed = grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::y_forward);
    QuESo_CHECK_EQUAL(next_untrimmed.next_id, 3UL);
    QuESo_CHECK_IS_FALSE(next_untrimmed.is_end);
    QuESo_CHECK(next_untrimmed.p_element != nullptr);
    QuESo_CHECK_EQUAL(next_untrimmed.p_element->GetId(), 3UL);
    QuESo_CHECK_IS_FALSE(next_untrimmed.p_element->IsTrimmed());

    const auto next_trimmed = grid.GetNextElement<ElementFilter::trimmed>(1, Direction::y_forward);
    QuESo_CHECK_EQUAL(next_trimmed.next_id, 3UL);
    QuESo_CHECK(next_trimmed.is_end);
    QuESo_CHECK(next_trimmed.p_element == nullptr);
}

BOOST_AUTO_TEST_CASE(GetNextElementFilteredZ)
{
    QuESo_INFO << "Testing :: Background Grid :: GetNextElementFilteredZ" << std::endl;

    auto grid = CreateMixedGrid();

    const auto next_all = grid.GetNextElementView(1, Direction::z_forward);
    QuESo_CHECK_EQUAL(next_all.next_id, 5UL);
    QuESo_CHECK_IS_FALSE(next_all.is_end);
    QuESo_CHECK_IS_FALSE(!next_all.element.has_value());
    QuESo_CHECK_EQUAL(next_all.element->GetId(), 5UL);
    QuESo_CHECK_IS_FALSE(!next_all.element->IsTrimmed());

    const auto next_untrimmed_view = grid.GetNextElementView<ElementFilter::untrimmed>(1, Direction::z_forward);
    QuESo_CHECK_EQUAL(next_untrimmed_view.next_id, 5UL);
    QuESo_CHECK(next_untrimmed_view.is_end);
    QuESo_CHECK_IS_FALSE(next_untrimmed_view.element.has_value());

    const auto next_trimmed_view = grid.GetNextElementView<ElementFilter::trimmed>(1, Direction::z_forward);
    QuESo_CHECK_EQUAL(next_trimmed_view.next_id, 5UL);
    QuESo_CHECK_IS_FALSE(next_trimmed_view.is_end);
    QuESo_CHECK_IS_FALSE(!next_trimmed_view.element.has_value());
    QuESo_CHECK_EQUAL(next_trimmed_view.element->GetId(), 5UL);

    const auto next_untrimmed = grid.GetNextElement<ElementFilter::untrimmed>(1, Direction::z_forward);
    QuESo_CHECK_EQUAL(next_untrimmed.next_id, 5UL);
    QuESo_CHECK(next_untrimmed.is_end);
    QuESo_CHECK(next_untrimmed.p_element == nullptr);

    const auto next_trimmed = grid.GetNextElement<ElementFilter::trimmed>(1, Direction::z_forward);
    QuESo_CHECK_EQUAL(next_trimmed.next_id, 5UL);
    QuESo_CHECK_IS_FALSE(next_trimmed.is_end);
    QuESo_CHECK(next_trimmed.p_element != nullptr);
    QuESo_CHECK_EQUAL(next_trimmed.p_element->GetId(), 5UL);
    QuESo_CHECK_IS_FALSE(!next_trimmed.p_element->IsTrimmed());
}

BOOST_AUTO_TEST_CASE(GetNextElementViewAllX)
{
    QuESo_INFO << "Testing :: Background Grid :: GetNextElementViewAllX" << std::endl;

    constexpr Vector3i number_of_elements{ 3, 4, 2 };
    const auto grid = CreateUntrimmedTraversalGrid(number_of_elements);

    IndexType current_id = 1;
    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 23UL);
    IndexType active_element_counter = 1;
    for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
        const auto next_result = grid.GetNextElementView(current_id, Direction::x_forward);
        const auto element = next_result.element;
        bool found = false;
        if (element.has_value()) {
            found = true;
            const auto reverse_result = grid.GetNextElementView(next_result.next_id, Direction::x_backward);
            if (next_result.is_end) { QuESo_CHECK(reverse_result.is_end); }
            QuESo_CHECK_EQUAL(current_id, reverse_result.next_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(next_result.next_id, i + 1);
        if (next_result.next_id == 17) {
            QuESo_CHECK(next_result.is_end);
            QuESo_CHECK_IS_FALSE(found);
        } else {
            QuESo_CHECK_EQUAL(element->GetId(), i + 1);
            QuESo_CHECK(found);

            if (current_id % 3 == 0) {
                QuESo_CHECK(next_result.is_end);
            } else {
                QuESo_CHECK_IS_FALSE(next_result.is_end);
            }
        }
        current_id = next_result.next_id;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23UL);
}

BOOST_AUTO_TEST_CASE(GetNextElementViewAllY)
{
    QuESo_INFO << "Testing :: Background Grid :: GetNextElementViewAllY" << std::endl;

    constexpr Vector3i number_of_elements{ 3, 4, 2 };
    const auto grid = CreateUntrimmedTraversalGrid(number_of_elements);

    IndexType current_id = 1;
    const std::vector<int> test_next_ids{ 1,  4,  7,  10, 2,  5,  8,  11, 3,  6,  9,  12,
                                          13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24 };
    const std::vector<int> test_local_ends{ 10, 11, 12, 22, 23, 24 };

    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 23UL);
    IndexType active_element_counter = 1;
    for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
        const auto next_result = grid.GetNextElementView(current_id, Direction::y_forward);
        const auto element = next_result.element;
        bool found = false;
        if (element.has_value()) {
            found = true;
            const auto reverse_result = grid.GetNextElementView(next_result.next_id, Direction::y_backward);
            if (next_result.is_end) { QuESo_CHECK(reverse_result.is_end); }
            QuESo_CHECK_EQUAL(current_id, reverse_result.next_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_result.next_id), test_next_ids[i]);

        if (next_result.next_id == 17) {
            QuESo_CHECK(next_result.is_end);
            QuESo_CHECK_IS_FALSE(found);
        } else {
            QuESo_CHECK(found);
            QuESo_CHECK_EQUAL(static_cast<int>(element->GetId()), test_next_ids[i]);

            if (Contains(test_local_ends, static_cast<int>(current_id))) {
                QuESo_CHECK(next_result.is_end);
            } else {
                QuESo_CHECK_IS_FALSE(next_result.is_end);
            }
        }
        current_id = next_result.next_id;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23UL);
}

BOOST_AUTO_TEST_CASE(GetNextElementViewAllZ)
{
    QuESo_INFO << "Testing :: Background Grid :: GetNextElementViewAllZ" << std::endl;

    constexpr Vector3i number_of_elements{ 3, 4, 2 };
    const auto grid = CreateUntrimmedTraversalGrid(number_of_elements);

    IndexType current_id = 1;
    const std::vector<int> test_next_ids{ 1, 13, 2, 14, 3, 15, 4,  16, 5,  17, 6,  18,
                                          7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24 };
    const std::vector<int> test_local_ends{ 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };

    QuESo_CHECK_EQUAL(grid.NumberOfActiveElements(), 23UL);
    IndexType active_element_counter = 1;
    for (IndexType i = 1; i <= grid.NumberOfActiveElements(); ++i) {
        const auto next_result = grid.GetNextElementView(current_id, Direction::z_forward);
        const auto element = next_result.element;
        bool found = false;
        if (element.has_value()) {
            found = true;
            const auto reverse_result = grid.GetNextElementView(next_result.next_id, Direction::z_backward);
            if (next_result.is_end) { QuESo_CHECK(reverse_result.is_end); }
            QuESo_CHECK_EQUAL(current_id, reverse_result.next_id);
            active_element_counter++;
        }
        QuESo_CHECK_EQUAL(static_cast<int>(next_result.next_id), test_next_ids[i]);

        if (next_result.next_id == 17) {
            QuESo_CHECK(next_result.is_end);
            QuESo_CHECK_IS_FALSE(found);
        } else {
            QuESo_CHECK(found);
            QuESo_CHECK_EQUAL(static_cast<int>(element->GetId()), test_next_ids[i]);

            if (Contains(test_local_ends, static_cast<int>(current_id))) {
                QuESo_CHECK(next_result.is_end);
            } else {
                QuESo_CHECK_IS_FALSE(next_result.is_end);
            }
        }
        current_id = next_result.next_id;
    }
    QuESo_CHECK_EQUAL(active_element_counter, 23UL);
}

BOOST_AUTO_TEST_SUITE_END()

}// namespace queso::Testing
