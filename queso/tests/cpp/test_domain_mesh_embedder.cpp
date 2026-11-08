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
#include <optional>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/embedding/domain_mesh_embedder.h"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/tests/cpp/helpers/test_data.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::testing {
namespace {

    [[nodiscard]] Unique<Dictionary<queso::key::MainValuesTypeTag>>
        MakeGridSettings(const BoundingBoxType& rBounds, const Vector3i& rNumberOfElements)
    {
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_grid_settings = (*p_settings)[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, rBounds.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, rBounds.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, rBounds.lower);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, rBounds.upper);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
        r_grid_settings.CheckRequired();
        return p_settings;
    }

    void CheckEqualStates(const embedding::ElementStates& rFirst, const embedding::ElementStates& rSecond)
    {
        REQUIRE((rFirst.size()) == (rSecond.size()));
        for (IndexType i = 0; i < rFirst.size(); ++i) { CHECK(rFirst[i] == rSecond[i]); }
    }

}  // namespace

TEST_SUITE("DomainMeshEmbedderTestSuite")
{

    TEST_CASE("DomainMeshEmbedderTestSuite - Classifies Interior And Trimmed Cells Deterministically")
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 3.0, 3.0, 3.0 }), { 3, 3, 3 });
        const GridIndexer grid_indexer(*p_settings);
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.75, 0.75, 0.75 }, { 2.25, 2.25, 2.25 });
        embedding::DomainMeshEmbedder embedder(mesh.View(), grid_indexer);

        const auto serial = embedder.Classify({ .partition_count = 1 });
        const auto serial_cached = embedder.Classify({ .partition_count = 3 });
        embedding::DomainMeshEmbedder partitioned_embedder(mesh.View(), grid_indexer);
        const auto partitioned = partitioned_embedder.Classify({ .partition_count = 3 });
        CheckEqualStates(serial, serial_cached);
        CheckEqualStates(serial, partitioned);

        const IndexType center = grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1);
        const IndexType corner = grid_indexer.GetVectorIndexFromMatrixIndices(0, 0, 0);
        CHECK(serial[center] == IntersectionState::inside);
        CHECK(serial[corner] == IntersectionState::trimmed);
    }

    TEST_CASE("DomainMeshEmbedderTestSuite - Grid Aligned Closed Cube Does Not Leak Or Create Trimmed Cells")
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 4.0, 4.0, 4.0 }), { 4, 4, 4 });
        const GridIndexer grid_indexer(*p_settings);
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 1.0, 1.0, 1.0 }, { 3.0, 3.0, 3.0 });
        embedding::DomainMeshEmbedder embedder(mesh.View(), grid_indexer);
        const auto states = embedder.Classify({ .partition_count = 2 });

        for (const auto state : states) { CHECK(state != IntersectionState::trimmed); }
        CHECK(states[grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1)] == IntersectionState::inside);
        CHECK(states[grid_indexer.GetVectorIndexFromMatrixIndices(2, 2, 2)] == IntersectionState::inside);
        CHECK(states[grid_indexer.GetVectorIndexFromMatrixIndices(0, 0, 0)] == IntersectionState::outside);
        CHECK(states[grid_indexer.GetVectorIndexFromMatrixIndices(3, 3, 3)] == IntersectionState::outside);
    }

    TEST_CASE("DomainMeshEmbedderTestSuite - Near Aligned Closed Cube Snaps To Closest Cell Faces")
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 4.0, 4.0, 4.0 }), { 4, 4, 4 });
        const GridIndexer grid_indexer(*p_settings);
        const BoundingBoxType cell_bounds = grid_indexer.GetBoundingBoxXYZFromIndex(0);
        const double offset = 0.5 * grid_indexer.GetGeometryTolerance().SnapDistance();
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox(
            { 1.0 + offset, 1.0 + offset, 1.0 + offset }, { 3.0 - offset, 3.0 - offset, 3.0 - offset }
        );
        embedding::DomainMeshEmbedder embedder(mesh.View(), grid_indexer);
        const auto states = embedder.Classify({ .partition_count = 2 });

        for (const auto state : states) { CHECK(state != IntersectionState::trimmed); }
        CHECK(states[grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1)] == IntersectionState::inside);
        CHECK(states[grid_indexer.GetVectorIndexFromMatrixIndices(0, 0, 0)] == IntersectionState::outside);
    }

    TEST_CASE("DomainMeshEmbedderTestSuite - Mixed Aligned And Trimmed Surfaces Classify Consistently")
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 4.0, 4.0, 4.0 }), { 4, 4, 4 });
        const GridIndexer grid_indexer(*p_settings);
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 1.0, 0.75, 0.75 }, { 3.0, 3.25, 3.25 });
        embedding::DomainMeshEmbedder serial_embedder(mesh.View(), grid_indexer);
        const auto serial_states = serial_embedder.Classify({ .partition_count = 1 });
        embedding::DomainMeshEmbedder partitioned_embedder(mesh.View(), grid_indexer);
        const auto partitioned_states = partitioned_embedder.Classify({ .partition_count = 4 });

        CheckEqualStates(serial_states, partitioned_states);
        CHECK(serial_states[grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1)] == IntersectionState::inside);
        CHECK(serial_states[grid_indexer.GetVectorIndexFromMatrixIndices(2, 2, 2)] == IntersectionState::inside);
        CHECK(serial_states[grid_indexer.GetVectorIndexFromMatrixIndices(0, 1, 1)] == IntersectionState::outside);
        CHECK(serial_states[grid_indexer.GetVectorIndexFromMatrixIndices(3, 1, 1)] == IntersectionState::outside);
        CHECK(serial_states[grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 1)] == IntersectionState::trimmed);
        CHECK(serial_states[grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 0)] == IntersectionState::trimmed);
    }

    TEST_CASE("DomainMeshEmbedderTestSuite - Trimmed Domain Consumes Section And Remains Self Contained")
    {
        auto p_settings = MakeGridSettings(MakeBox({ 0.0, 0.0, 0.0 }, { 3.0, 3.0, 3.0 }), { 3, 3, 3 });
        const GridIndexer grid_indexer(*p_settings);
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.75, 0.75, 0.75 }, { 2.25, 2.25, 2.25 });
        std::optional<TrimmedDomain> domain;

        {
            embedding::DomainMeshEmbedder embedder(mesh.View(), grid_indexer);
            const IndexType cell_index = grid_indexer.GetVectorIndexFromMatrixIndices(0, 0, 0);
            if constexpr (!NOTDEBUG) {
                CHECK_THROWS_AS((void)embedder.MakeTrimmedDomain(grid_indexer.NumberOfElements(), 0), queso::Exception);
                CHECK_THROWS_AS((void)embedder.MakeTrimmedDomain(cell_index, 0), queso::Exception);
            }
            static_cast<void>(embedder.Classify());
            if constexpr (!NOTDEBUG) {
                CHECK_THROWS_AS(
                    (void)embedder.MakeTrimmedDomain(grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1), 0),
                    queso::Exception
                );
            }
            domain.emplace(embedder.MakeTrimmedDomain(cell_index, 0));
            if constexpr (!NOTDEBUG) {
                CHECK_THROWS_AS((void)embedder.MakeTrimmedDomain(cell_index, 0), queso::Exception);
            }
        }

        CHECK(domain->IsInside(PointType{ 0.9, 0.9, 0.9 }));
        CHECK(!domain->IsInside(PointType{ 0.1, 0.1, 0.1 }));
        CHECK((std::abs(MeshUtilities::Volume(domain->GetBoundaryMesh()) - 0.25 * 0.25 * 0.25)) < (1e-12));
    }

    TEST_CASE("DomainMeshEmbedderTestSuite - Steering Knuckle Face Topology Matches Direct Cell Centers")
    {
        auto p_settings = MakeGridSettings(MakeBox({ -130.0, -110.0, -110.0 }, { 20.0, 190.0, 190.0 }), { 10, 25, 25 });
        const GridIndexer grid_indexer(*p_settings);
        TriangleMesh mesh;
        IO::ReadMeshFromSTL(mesh, TestDataDirectory() + "/steering_knuckle/steering_knuckle.stl");
        embedding::DomainMeshEmbedder embedder(mesh.View(), grid_indexer);
        const auto states = embedder.Classify({ .partition_count = 1 });
        const embedding::MeshOperator mesh_operator(mesh.View(), grid_indexer.GetGeometryTolerance());
        for (IndexType cell_index = 0; cell_index < states.size(); ++cell_index) {
            if (states[cell_index] == IntersectionState::trimmed) { continue; }
            const BoundingBoxType bounds = grid_indexer.GetBoundingBoxXYZFromIndex(cell_index);
            const bool is_inside = mesh_operator.IsInside(0.5 * (bounds.lower + bounds.upper));
            CHECK((states[cell_index] == IntersectionState::inside) == is_inside);
        }
    }

}  // TEST_SUITE

}  // namespace queso::testing
