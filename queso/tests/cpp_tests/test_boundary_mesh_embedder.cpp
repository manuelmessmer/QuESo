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

//// STL includes
#include <array>
#include <cmath>
#include <optional>
#include <span>

//// External includes
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/embedding/boundary_mesh_embedder.h"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {
namespace {

    using BackgroundGridType = embedding::BoundaryMeshEmbedder::BackgroundGridType;
    using UntrimmedElementType = BackgroundGridType::UntrimmedElementType;

    struct UntrimmedBuilder
    {
        static constexpr BackgroundGridType::ElementFilter Builds = BackgroundGridType::ElementFilter::untrimmed;

        [[nodiscard]] std::optional<UntrimmedElementType> Build(IndexType CellIndex, const ElementBounds& rBounds)
        { return UntrimmedElementType(CellIndex + 1, rBounds); }
    };

    [[nodiscard]] Unique<Dictionary<queso::key::MainValuesTypeTag>> MakeSettings(const Vector3i& rNumberOfElements)
    {
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_grid_settings = (*p_settings)[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ 0.0, 0.0, 0.0 });
        r_grid_settings.SetValue(
            BackgroundGridSettings::upper_bound_xyz,
            PointType{ static_cast<double>(rNumberOfElements[0]),
                       static_cast<double>(rNumberOfElements[1]),
                       static_cast<double>(rNumberOfElements[2]) }
        );
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ 0.0, 0.0, 0.0 });
        r_grid_settings.SetValue(
            BackgroundGridSettings::upper_bound_uvw,
            PointType{ static_cast<double>(rNumberOfElements[0]),
                       static_cast<double>(rNumberOfElements[1]),
                       static_cast<double>(rNumberOfElements[2]) }
        );
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid_settings.CheckRequired();
        return p_settings;
    }

    [[nodiscard]] BackgroundGridType
        MakeGrid(const Vector3i& rNumberOfElements, std::span<const IndexType> rActiveCells)
    {
        const auto p_settings = MakeSettings(rNumberOfElements);
        BackgroundGridType grid(*p_settings);
        UntrimmedBuilder builder;
        for (const IndexType cell_index : rActiveCells) {
            QuESo_CHECK(grid.MakeElement(builder, cell_index));
        }
        grid.LockElements();
        return grid;
    }

    void AddTriangle(TriangleMesh& rMesh, const PointType& rP1, const PointType& rP2, const PointType& rP3)
    {
        const IndexType first = rMesh.AddVertex(rP1);
        const IndexType second = rMesh.AddVertex(rP2);
        const IndexType third = rMesh.AddVertex(rP3);
        rMesh.AddTriangle({ first, second, third }, { 1.0, 0.0, 0.0 });
    }

    [[nodiscard]] TriangleMesh MakeTriangleAtX(double X)
    {
        TriangleMesh mesh;
        AddTriangle(mesh, { X, 0.2, 0.2 }, { X, 0.8, 0.2 }, { X, 0.2, 0.8 });
        return mesh;
    }

    [[nodiscard]] TriangleMesh MakeTriangleOnCenterCellFace(SignedAxis Face)
    {
        const IndexType face_index = static_cast<IndexType>(Face);
        const IndexType axis = face_index / 2;
        const double coordinate = face_index % 2 == 0 ? 1.0 : 2.0;
        std::array<PointType, 3> points{ PointType{ 1.2, 1.2, 1.2 },
                                         PointType{ 1.2, 1.2, 1.2 },
                                         PointType{ 1.2, 1.2, 1.2 } };
        std::array<IndexType, 2> tangential_axes{};
        IndexType tangential_index = 0;
        for (IndexType candidate_axis = 0; candidate_axis < 3; ++candidate_axis) {
            if (candidate_axis != axis) { tangential_axes[tangential_index++] = candidate_axis; }
        }
        points[1][tangential_axes[0]] = 1.8;
        points[2][tangential_axes[1]] = 1.8;
        for (auto& r_point : points) { r_point[axis] = coordinate; }
        TriangleMesh mesh;
        AddTriangle(mesh, points[0], points[1], points[2]);
        return mesh;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(BoundaryMeshEmbedderTestSuite)

BOOST_AUTO_TEST_CASE(RequiresLockedBackgroundGrid)
{
    if constexpr (!NOTDEBUG) {
        const auto p_settings = MakeSettings({ 2, 1, 1 });
        const BackgroundGridType grid(*p_settings);
        const TriangleMesh mesh = MakeTriangleAtX(0.5);
        BOOST_CHECK_THROW((void)embedding::BoundaryMeshEmbedder(mesh.View(), grid), queso::Exception);
    }
}

BOOST_AUTO_TEST_CASE(UsesActiveCanonicalOwnerOnce)
{
    constexpr std::array<IndexType, 1> active_cells{ 1 };
    const auto grid = MakeGrid({ 2, 1, 1 }, active_cells);
    const TriangleMesh mesh = MakeTriangleAtX(1.0);
    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);

    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 1);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices().front(), 1);
    const auto section = embedder.TakeBoundarySection(1);
    QuESo_CHECK_EQUAL(section.NumOfTriangles(), 1);
    if constexpr (!NOTDEBUG) { BOOST_CHECK_THROW((void)embedder.TakeBoundarySection(1), queso::Exception); }
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices().front(), 1);
}

BOOST_AUTO_TEST_CASE(FallsBackToIncidentActiveCell)
{
    constexpr std::array<IndexType, 1> active_cells{ 0 };
    const auto grid = MakeGrid({ 2, 1, 1 }, active_cells);
    const TriangleMesh mesh = MakeTriangleAtX(1.0);
    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);

    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 1);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices().front(), 0);
    const auto section = embedder.TakeBoundarySection(0);
    QuESo_CHECK_EQUAL(section.NumOfTriangles(), 1);
    QuESo_CHECK_POINT_NEAR(section.Normal(0), mesh.Normal(0), 1e-12);
}

BOOST_AUTO_TEST_CASE(RetainsCanonicalParentWhenNeitherAdjacentCellIsActive)
{
    constexpr std::array<IndexType, 1> active_cells{ 2 };
    const auto grid = MakeGrid({ 4, 1, 1 }, active_cells);
    const TriangleMesh mesh = MakeTriangleAtX(1.0);
    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);

    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 1);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices().front(), 1);
    QuESo_CHECK_EQUAL(embedder.TakeBoundarySection(1).NumOfTriangles(), IndexType{ 1 });
    if constexpr (!NOTDEBUG) { BOOST_CHECK_THROW((void)embedder.TakeBoundarySection(2), queso::Exception); }
}

BOOST_AUTO_TEST_CASE(RetainsInGridGeometryOutsideActiveCellHalo)
{
    constexpr std::array<IndexType, 1> active_cells{ 0 };
    const auto grid = MakeGrid({ 4, 1, 1 }, active_cells);
    const TriangleMesh mesh = MakeTriangleAtX(3.5);
    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);

    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 1);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices().front(), 3);
    QuESo_CHECK_EQUAL(embedder.TakeBoundarySection(3).NumOfTriangles(), IndexType{ 1 });
}

BOOST_AUTO_TEST_CASE(ReturnsSortedUniqueParents)
{
    constexpr std::array<IndexType, 2> active_cells{ 0, 2 };
    const auto grid = MakeGrid({ 3, 1, 1 }, active_cells);
    TriangleMesh mesh = MakeTriangleAtX(0.5);
    AddTriangle(mesh, { 2.5, 0.2, 0.2 }, { 2.5, 0.8, 0.2 }, { 2.5, 0.2, 0.8 });
    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);

    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 2);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[0], 0);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[1], 2);
    if constexpr (!NOTDEBUG) { BOOST_CHECK_THROW((void)embedder.TakeBoundarySection(1), queso::Exception); }
}

BOOST_AUTO_TEST_CASE(AssignsAllSixFacesToCanonicalOrFallbackParent)
{
    const auto p_settings = MakeSettings({ 3, 3, 3 });
    const GridIndexer indexer(*p_settings);
    for (const SignedAxis face : EnumRange<SignedAxis>()) {
        const IndexType face_index = static_cast<IndexType>(face);
        const IndexType axis = face_index / 2;
        const bool is_upper = face_index % 2 == 1;
        Vector3i canonical_indices{ 1, 1, 1 };
        if (is_upper) { ++canonical_indices[axis]; }
        const IndexType canonical = indexer.GetVectorIndexFromMatrixIndices(canonical_indices);
        const TriangleMesh mesh = MakeTriangleOnCenterCellFace(face);

        const std::array canonical_active{ canonical };
        const auto canonical_grid = MakeGrid({ 3, 3, 3 }, canonical_active);
        embedding::BoundaryMeshEmbedder canonical_embedder(mesh.View(), canonical_grid);
        BOOST_REQUIRE_EQUAL(canonical_embedder.GetParentCellIndices().size(), 1);
        QuESo_CHECK_EQUAL(canonical_embedder.GetParentCellIndices().front(), canonical);

        Vector3i fallback_indices = canonical_indices;
        --fallback_indices[axis];
        const IndexType fallback = indexer.GetVectorIndexFromMatrixIndices(fallback_indices);
        const std::array fallback_active{ fallback };
        const auto fallback_grid = MakeGrid({ 3, 3, 3 }, fallback_active);
        embedding::BoundaryMeshEmbedder fallback_embedder(mesh.View(), fallback_grid);
        BOOST_REQUIRE_EQUAL(fallback_embedder.GetParentCellIndices().size(), 1);
        QuESo_CHECK_EQUAL(fallback_embedder.GetParentCellIndices().front(), fallback);
    }
}

BOOST_AUTO_TEST_CASE(AssignsMixedInteriorAndFaceSectionsIndependently)
{
    constexpr std::array<IndexType, 2> active_cells{ 0, 1 };
    const auto grid = MakeGrid({ 2, 1, 1 }, active_cells);
    TriangleMesh mesh = MakeTriangleAtX(1.0);
    AddTriangle(mesh, { 0.2, 0.2, 0.2 }, { 0.8, 0.3, 0.2 }, { 0.3, 0.8, 0.7 });
    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);

    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 2);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[0], 0);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[1], 1);
    const auto interior = embedder.TakeBoundarySection(0);
    const auto face = embedder.TakeBoundarySection(1);
    QuESo_CHECK_EQUAL(interior.NumOfTriangles(), 1);
    QuESo_CHECK_EQUAL(face.NumOfTriangles(), 1);
}

BOOST_AUTO_TEST_CASE(CombinesInteriorAndFaceSectionsForOneParent)
{
    constexpr std::array<IndexType, 2> active_cells{ 0, 1 };
    const auto grid = MakeGrid({ 2, 1, 1 }, active_cells);
    TriangleMesh mesh = MakeTriangleAtX(1.0);
    AddTriangle(mesh, { 1.2, 0.2, 0.2 }, { 1.8, 0.3, 0.2 }, { 1.3, 0.8, 0.7 });

    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);
    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 1);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices().front(), 1);
    const auto section = embedder.TakeBoundarySection(1);
    QuESo_CHECK_EQUAL(section.NumOfTriangles(), IndexType{ 2 });
    QuESo_CHECK_LT(std::abs(MeshUtilities::Area(section.View()) - MeshUtilities::Area(mesh.View())), 1e-12);
}

BOOST_AUTO_TEST_CASE(AssignsExteriorFaceSectionsToOnlyAdjacentCells)
{
    constexpr std::array<IndexType, 0> active_cells{};
    const auto grid = MakeGrid({ 2, 1, 1 }, active_cells);
    TriangleMesh mesh = MakeTriangleAtX(0.0);
    AddTriangle(mesh, { 2.0, 0.2, 0.2 }, { 2.0, 0.8, 0.2 }, { 2.0, 0.2, 0.8 });

    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);
    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 2);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[0], 0);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[1], 1);
    double area{};
    for (const IndexType cell_index : embedder.GetParentCellIndices()) {
        const auto section = embedder.TakeBoundarySection(cell_index);
        QuESo_CHECK_EQUAL(section.NumOfTriangles(), IndexType{ 1 });
        area += MeshUtilities::Area(section.View());
    }
    QuESo_CHECK_LT(std::abs(area - MeshUtilities::Area(mesh.View())), 1e-12);
}

BOOST_AUTO_TEST_CASE(AssignsAsymmetricSplitPiecesExactlyOnce)
{
    constexpr std::array<IndexType, 0> active_cells{};
    const auto grid = MakeGrid({ 2, 1, 1 }, active_cells);
    TriangleMesh mesh;
    AddTriangle(mesh, { 0.25, 0.2, 0.5 }, { 1.75, 0.2, 0.5 }, { 0.5, 0.8, 0.5 });

    embedding::BoundaryMeshEmbedder embedder(mesh.View(), grid);
    BOOST_REQUIRE_EQUAL(embedder.GetParentCellIndices().size(), 2);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[0], 0);
    QuESo_CHECK_EQUAL(embedder.GetParentCellIndices()[1], 1);
    const auto left = embedder.TakeBoundarySection(0);
    const auto right = embedder.TakeBoundarySection(1);
    QuESo_CHECK_LT(std::abs(MeshUtilities::Area(left.View()) - 0.315), 1e-12);
    QuESo_CHECK_LT(std::abs(MeshUtilities::Area(right.View()) - 0.135), 1e-12);
    QuESo_CHECK_LT(
        std::abs(
            MeshUtilities::Area(left.View()) + MeshUtilities::Area(right.View()) - MeshUtilities::Area(mesh.View())
        ),
        1e-12
    );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
