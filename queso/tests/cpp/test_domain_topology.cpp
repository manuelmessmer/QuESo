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

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/embedding/domain_topology.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso::testing {
namespace {

    [[nodiscard]] Unique<Dictionary<queso::key::MainValuesTypeTag>> MakeGridSettings()
    {
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_grid_settings = (*p_settings)[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ 0.0, 0.0, 0.0 });
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 2.0, 2.0, 2.0 });
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ 0.0, 0.0, 0.0 });
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 2.0, 2.0, 2.0 });
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{ 2, 2, 2 });
        r_grid_settings.CheckRequired();
        return p_settings;
    }

}  // namespace

TEST_SUITE("DomainTopologyTestSuite")
{

    TEST_CASE("DomainTopologyTestSuite - Domain Topology Initializes Cell Face And Vote Storage")
    {
        const auto p_settings = MakeGridSettings();
        const GridIndexer grid_indexer(*p_settings);
        embedding::detail::DomainTopology topology(grid_indexer);
        using Direction = GridIndexer::Direction;
        using State = embedding::detail::GridFaceState;

        for (IndexType cell_index = 0; cell_index < grid_indexer.NumberOfElements(); ++cell_index) {
            const auto& r_cell = topology.cell_topology[cell_index];
            CHECK(!r_cell.is_trimmed);
            REQUIRE(!r_cell.surface_section.has_value());
            CHECK(!r_cell.input_consumed);
            for (const Direction direction : EnumRange<Direction>()) {
                const State expected = grid_indexer.IsEnd(cell_index, direction) ? State::blocked : State::traversable;
                CHECK(topology.face_topology.Get(cell_index, direction) == expected);
                CHECK((topology.face_votes.Get(cell_index, direction)) == (0));
            }
        }
    }

    TEST_CASE("DomainTopologyTestSuite - Canonical Faces Resolve Opposite Directions To One State")
    {
        const auto p_settings = MakeGridSettings();
        const GridIndexer grid_indexer(*p_settings);
        embedding::detail::GridFaceStates states(grid_indexer);
        using Direction = GridIndexer::Direction;
        using State = embedding::detail::GridFaceState;

        CHECK(states.Get(0, Direction::x_forward) == State::traversable);
        states.Get(0, Direction::x_forward) = State::blocked;
        CHECK(states.Get(1, Direction::x_backward) == State::blocked);

        states.Get(0, Direction::y_forward) = State::blocked;
        CHECK(states.Get(2, Direction::y_backward) == State::blocked);

        states.Get(0, Direction::z_forward) = State::blocked;
        CHECK(states.Get(4, Direction::z_backward) == State::blocked);
    }

    TEST_CASE("DomainTopologyTestSuite - Exterior Faces Are Blocked And Votes Remain Directional")
    {
        const auto p_settings = MakeGridSettings();
        const GridIndexer grid_indexer(*p_settings);
        const embedding::detail::GridFaceStates states(grid_indexer);
        embedding::detail::GridFaceVotes votes(grid_indexer);
        using Direction = GridIndexer::Direction;
        using State = embedding::detail::GridFaceState;

        CHECK(states.Get(0, Direction::x_backward) == State::blocked);
        CHECK(states.Get(0, Direction::y_backward) == State::blocked);
        CHECK(states.Get(0, Direction::z_backward) == State::blocked);
        CHECK(states.Get(7, Direction::x_forward) == State::blocked);
        CHECK(states.Get(7, Direction::y_forward) == State::blocked);
        CHECK(states.Get(7, Direction::z_forward) == State::blocked);

        votes.Get(0, Direction::x_forward) = 1;
        votes.Get(1, Direction::x_backward) = -1;
        const GridFaceId face = grid_indexer.GetFace(0, Direction::x_forward);
        CHECK((votes.Get(0, Direction::x_forward)) == (1));
        CHECK((votes.Get(1, Direction::x_backward)) == (-1));
        CHECK((votes.Get(face, GridFaceSide::negative)) == (1));
        CHECK((votes.Get(face, GridFaceSide::positive)) == (-1));
    }

}  // TEST_SUITE

}  // namespace queso::testing
