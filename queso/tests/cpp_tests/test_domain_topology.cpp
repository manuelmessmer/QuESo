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

//// Project includes
#include "queso/embedding/domain_topology.hpp"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso::Testing {
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

BOOST_AUTO_TEST_SUITE(DomainTopologyTestSuite)

BOOST_AUTO_TEST_CASE(DomainTopologyInitializesCellFaceAndVoteStorage)
{
    const auto p_settings = MakeGridSettings();
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    using Direction = GridIndexer::Direction;
    using State = embedding::detail::GridFaceState;

    for (IndexType cell_index = 0; cell_index < grid_indexer.NumberOfElements(); ++cell_index) {
        const auto& r_cell = topology.cell_topology[cell_index];
        QuESo_CHECK(!r_cell.is_trimmed);
        QuESo_CHECK(!r_cell.surface_section.has_value());
        QuESo_CHECK(!r_cell.input_consumed);
        for (const Direction direction : EnumRange<Direction>()) {
            const State expected = grid_indexer.IsEnd(cell_index, direction) ? State::blocked : State::traversable;
            QuESo_CHECK(topology.face_topology.Get(cell_index, direction) == expected);
            QuESo_CHECK_EQUAL(topology.face_votes.Get(cell_index, direction), 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(CanonicalFacesResolveOppositeDirectionsToOneState)
{
    const auto p_settings = MakeGridSettings();
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::GridFaceStates states(grid_indexer);
    using Direction = GridIndexer::Direction;
    using State = embedding::detail::GridFaceState;

    QuESo_CHECK(states.Get(0, Direction::x_forward) == State::traversable);
    states.Get(0, Direction::x_forward) = State::blocked;
    QuESo_CHECK(states.Get(1, Direction::x_backward) == State::blocked);

    states.Get(0, Direction::y_forward) = State::blocked;
    QuESo_CHECK(states.Get(2, Direction::y_backward) == State::blocked);

    states.Get(0, Direction::z_forward) = State::blocked;
    QuESo_CHECK(states.Get(4, Direction::z_backward) == State::blocked);
}

BOOST_AUTO_TEST_CASE(ExteriorFacesAreBlockedAndVotesRemainDirectional)
{
    const auto p_settings = MakeGridSettings();
    const GridIndexer grid_indexer(*p_settings);
    const embedding::detail::GridFaceStates states(grid_indexer);
    embedding::detail::GridFaceVotes votes(grid_indexer);
    using Direction = GridIndexer::Direction;
    using State = embedding::detail::GridFaceState;

    QuESo_CHECK(states.Get(0, Direction::x_backward) == State::blocked);
    QuESo_CHECK(states.Get(0, Direction::y_backward) == State::blocked);
    QuESo_CHECK(states.Get(0, Direction::z_backward) == State::blocked);
    QuESo_CHECK(states.Get(7, Direction::x_forward) == State::blocked);
    QuESo_CHECK(states.Get(7, Direction::y_forward) == State::blocked);
    QuESo_CHECK(states.Get(7, Direction::z_forward) == State::blocked);

    votes.Get(0, Direction::x_forward) = 1;
    votes.Get(1, Direction::x_backward) = -1;
    const GridFaceId face = grid_indexer.GetFace(0, Direction::x_forward);
    QuESo_CHECK_EQUAL(votes.Get(0, Direction::x_forward), 1);
    QuESo_CHECK_EQUAL(votes.Get(1, Direction::x_backward), -1);
    QuESo_CHECK_EQUAL(votes.Get(face, GridFaceSide::negative), 1);
    QuESo_CHECK_EQUAL(votes.Get(face, GridFaceSide::positive), -1);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
