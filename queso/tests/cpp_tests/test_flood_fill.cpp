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
#include "queso/embedding/flood_fill.h"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso::Testing {
namespace {

    [[nodiscard]] Unique<Dictionary<queso::key::MainValuesTypeTag>> MakeGridSettings(const Vector3i& rNumberOfElements)
    {
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_grid_settings = (*p_settings)[MainSettings::background_grid_settings];
        const PointType upper{ static_cast<double>(rNumberOfElements[0]),
                               static_cast<double>(rNumberOfElements[1]),
                               static_cast<double>(rNumberOfElements[2]) };
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ 0.0, 0.0, 0.0 });
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ 0.0, 0.0, 0.0 });
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
        r_grid_settings.CheckRequired();
        return p_settings;
    }

    void BlockAllFaces(embedding::detail::DomainTopology& rTopology, IndexType CellIndex)
    {
        for (const auto direction : EnumRange<GridIndexer::Direction>()) {
            rTopology.face_topology.Get(CellIndex, direction) = embedding::detail::GridFaceState::blocked;
        }
    }

    void CheckEqualStates(const embedding::ElementStates& rFirst, const embedding::ElementStates& rSecond)
    {
        BOOST_REQUIRE_EQUAL(rFirst.size(), rSecond.size());
        for (IndexType i = 0; i < rFirst.size(); ++i) { QuESo_CHECK(rFirst[i] == rSecond[i]); }
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(FloodFillTestSuite)

BOOST_AUTO_TEST_CASE(AllClearGridIsExteriorOutside)
{
    const auto p_settings = MakeGridSettings({ 2, 2, 2 });
    const GridIndexer grid_indexer(*p_settings);
    const embedding::detail::DomainTopology topology(grid_indexer);
    const auto states = embedding::detail::FloodFill(topology, { .partition_count = 1 });

    for (const auto state : states) { QuESo_CHECK(state == IntersectionState::outside); }
}

BOOST_AUTO_TEST_CASE(OnlyMarkedCellsAreTrimmed)
{
    const auto p_settings = MakeGridSettings({ 3, 1, 1 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    topology.cell_topology[1].is_trimmed = true;
    topology.face_topology.Get(0, GridIndexer::Direction::x_forward) = embedding::detail::GridFaceState::blocked;
    topology.face_topology.Get(1, GridIndexer::Direction::x_forward) = embedding::detail::GridFaceState::blocked;

    const auto states = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    QuESo_CHECK(states[0] == IntersectionState::outside);
    QuESo_CHECK(states[1] == IntersectionState::trimmed);
    QuESo_CHECK(states[2] == IntersectionState::outside);
}

BOOST_AUTO_TEST_CASE(BlockedFillableCellUsesComponentVoteAndIsNotTrimmed)
{
    const auto p_settings = MakeGridSettings({ 3, 3, 3 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    const IndexType center = grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1);
    BlockAllFaces(topology, center);
    for (const auto direction : EnumRange<GridIndexer::Direction>()) { topology.face_votes.Get(center, direction) = 1; }

    const auto states = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    QuESo_CHECK(states[center] == IntersectionState::inside);
}

BOOST_AUTO_TEST_CASE(ZeroAndNegativeVotesClassifyOutside)
{
    const auto p_settings = MakeGridSettings({ 3, 3, 3 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    const IndexType center = grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1);
    BlockAllFaces(topology, center);

    auto states = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    QuESo_CHECK(states[center] == IntersectionState::outside);

    topology.face_votes.Get(center, GridIndexer::Direction::x_forward) = -1;
    states = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    QuESo_CHECK(states[center] == IntersectionState::outside);
}

BOOST_AUTO_TEST_CASE(ExteriorConnectivityOverridesPositiveVotes)
{
    const auto p_settings = MakeGridSettings({ 2, 1, 1 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    topology.face_topology.Get(0, GridIndexer::Direction::x_forward) = embedding::detail::GridFaceState::blocked;
    topology.face_votes.Get(0, GridIndexer::Direction::x_forward) = 100;

    const auto states = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    QuESo_CHECK(states[0] == IntersectionState::outside);
}

BOOST_AUTO_TEST_CASE(VotesOnTraversableFacesAreIgnored)
{
    const auto p_settings = MakeGridSettings({ 4, 3, 3 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    const IndexType left = grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1);
    const IndexType right = grid_indexer.GetVectorIndexFromMatrixIndices(2, 1, 1);
    BlockAllFaces(topology, left);
    BlockAllFaces(topology, right);
    topology.face_topology.Get(left, GridIndexer::Direction::x_forward) = embedding::detail::GridFaceState::traversable;
    topology.face_votes.Get(left, GridIndexer::Direction::x_forward) = 100;
    topology.face_votes.Get(right, GridIndexer::Direction::x_backward) = 100;

    const auto states = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    QuESo_CHECK(states[left] == IntersectionState::outside);
    QuESo_CHECK(states[right] == IntersectionState::outside);
}

BOOST_AUTO_TEST_CASE(PartitionInterfacesCollectBothDirectedVotesDeterministically)
{
    const auto p_settings = MakeGridSettings({ 4, 3, 3 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    const IndexType left = grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1);
    const IndexType right = grid_indexer.GetVectorIndexFromMatrixIndices(2, 1, 1);
    BlockAllFaces(topology, left);
    BlockAllFaces(topology, right);
    topology.face_votes.Get(left, GridIndexer::Direction::x_forward) = 1;
    topology.face_votes.Get(right, GridIndexer::Direction::x_backward) = 1;

    const auto serial = embedding::detail::FloodFill(topology, { .partition_count = 1 });
    const auto two_partitions = embedding::detail::FloodFill(topology, { .partition_count = 2 });
    const auto four_partitions = embedding::detail::FloodFill(topology, { .partition_count = 4 });
    CheckEqualStates(serial, two_partitions);
    CheckEqualStates(serial, four_partitions);
    QuESo_CHECK(serial[left] == IntersectionState::inside);
    QuESo_CHECK(serial[right] == IntersectionState::inside);
}

BOOST_AUTO_TEST_CASE(RejectsInvalidPartitionCounts)
{
    const auto p_settings = MakeGridSettings({ 3, 2, 1 });
    const GridIndexer grid_indexer(*p_settings);
    const embedding::detail::DomainTopology topology(grid_indexer);
    BOOST_CHECK_THROW((void)embedding::detail::FloodFill(topology, { .partition_count = 0 }), queso::Exception);
    BOOST_CHECK_THROW((void)embedding::detail::FloodFill(topology, { .partition_count = 4 }), queso::Exception);
}

BOOST_AUTO_TEST_CASE(RejectsMalformedTopology)
{
    const auto p_settings = MakeGridSettings({ 1, 1, 4 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    topology.cell_topology.pop_back();
    BOOST_CHECK_THROW((void)embedding::detail::FloodFill(topology), queso::Exception);
}

BOOST_AUTO_TEST_CASE(RejectsTraversableFacesAdjacentToTrimmedCells)
{
    const auto p_settings = MakeGridSettings({ 2, 1, 1 });
    const GridIndexer grid_indexer(*p_settings);
    embedding::detail::DomainTopology topology(grid_indexer);
    topology.cell_topology[0].is_trimmed = true;
    BOOST_CHECK_THROW((void)embedding::detail::FloodFill(topology), queso::Exception);
}

BOOST_AUTO_TEST_CASE(PartitioningAlongYAndZMatchesSerialFill)
{
    for (const Vector3i counts : { Vector3i{ 2, 5, 2 }, Vector3i{ 2, 2, 5 } }) {
        const auto p_settings = MakeGridSettings(counts);
        const GridIndexer grid_indexer(*p_settings);
        const embedding::detail::DomainTopology topology(grid_indexer);
        const auto serial = embedding::detail::FloodFill(topology, { .partition_count = 1 });
        const auto partitioned = embedding::detail::FloodFill(topology, { .partition_count = 5 });
        CheckEqualStates(serial, partitioned);
    }
}

BOOST_AUTO_TEST_CASE(OptionsSupportEquality)
{
    QuESo_CHECK(embedding::FloodFillOptions{} == embedding::FloodFillOptions{});
    QuESo_CHECK(
        embedding::FloodFillOptions{ .partition_count = 2 } != embedding::FloodFillOptions{ .partition_count = 1 }
    );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
