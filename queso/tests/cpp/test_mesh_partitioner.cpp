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

//// External includes
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/embedding/convex_polygon.h"
#include "queso/embedding/mesh_partitioner.h"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {
namespace {

    using DomainMeshPartitioner = embedding::MeshPartitioner<embedding::CellProduct::Domain>;
    using SurfaceMeshPartitioner = embedding::MeshPartitioner<embedding::CellProduct::Surface>;

    [[nodiscard]] GridIndexer MakeGridIndexer(const BoundingBoxType& rBounds, const Vector3i& rNumberOfElements)
    {
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_grid = (*p_settings)[MainSettings::background_grid_settings];
        r_grid.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        r_grid.SetValue(BackgroundGridSettings::lower_bound_xyz, rBounds.lower);
        r_grid.SetValue(BackgroundGridSettings::upper_bound_xyz, rBounds.upper);
        r_grid.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ 0.0, 0.0, 0.0 });
        r_grid.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 1.0, 1.0 });
        r_grid.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 2, 2, 2 });
        r_grid.SetValue(BackgroundGridSettings::number_of_elements, rNumberOfElements);
        r_grid.CheckRequired();
        return GridIndexer(*p_settings);
    }

    [[nodiscard]] GridIndexer MakeGridIndexer(const Vector3i& rNumberOfElements = { 2, 1, 1 })
    {
        return MakeGridIndexer(
            MakeBox(
                { 0.0, 0.0, 0.0 },
                { static_cast<double>(rNumberOfElements[0]),
                  static_cast<double>(rNumberOfElements[1]),
                  static_cast<double>(rNumberOfElements[2]) }
            ),
            rNumberOfElements
        );
    }

    void AddTriangle(TriangleMesh& rMesh, PointView rFirst, PointView rSecond, PointView rThird)
    {
        const IndexType first = rMesh.AddVertex({ rFirst[0], rFirst[1], rFirst[2] });
        const IndexType second = rMesh.AddVertex({ rSecond[0], rSecond[1], rSecond[2] });
        const IndexType third = rMesh.AddVertex({ rThird[0], rThird[1], rThird[2] });
        const Vector3d normal = Math::Cross(rSecond - rFirst, rThird - rFirst);
        rMesh.AddTriangle({ first, second, third }, normal / Math::Norm(normal));
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(MeshPartitionerTestSuite)

BOOST_AUTO_TEST_CASE(ConvexPolygonSplitConservesArea)
{
    embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
    polygon.AddVertex(PointType{ 0.0, 0.0, 0.5 });
    polygon.AddVertex(PointType{ 2.0, 0.0, 0.5 });
    polygon.AddVertex(PointType{ 1.0, 1.0, 0.5 });

    const auto split = embedding::detail::SplitByPlane(
        polygon, 0, 1.0, GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
    );
    BOOST_REQUIRE(split.negative.has_value());
    BOOST_REQUIRE(split.positive.has_value());
    QuESo_CHECK_LT(std::abs(split.negative->Area() + split.positive->Area() - polygon.Area()), 1e-12);
    for (const auto& r_piece : { *split.negative, *split.positive }) {
        QuESo_CHECK_EQUAL(r_piece.NumberOfVertices(), IndexType{ 3 });
        QuESo_CHECK_POINT_NEAR(r_piece.Normal(), polygon.Normal(), 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(ConvexPolygonOnPlaneProducesBothSides)
{
    embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
    polygon.AddVertex(PointType{ 1.0, 0.0, 0.0 });
    polygon.AddVertex(PointType{ 1.0, 1.0, 0.0 });
    polygon.AddVertex(PointType{ 1.0, 0.0, 1.0 });

    const auto split = embedding::detail::SplitByPlane(
        polygon, 0, 1.0, GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
    );
    BOOST_REQUIRE(split.negative.has_value());
    BOOST_REQUIRE(split.positive.has_value());
    for (const auto& r_piece : { *split.negative, *split.positive }) {
        QuESo_CHECK_LT(std::abs(r_piece.Area() - polygon.Area()), 1e-12);
        QuESo_CHECK_EQUAL(r_piece.NumberOfVertices(), polygon.NumberOfVertices());
        QuESo_CHECK_POINT_NEAR(r_piece.Normal(), polygon.Normal(), 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(ConvexPolygonSnapsAndDeduplicatesVertices)
{
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1e11, .coordinate_scale = 1.0 });
    const double snap_distance = tolerance.SnapDistance();
    embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
    polygon.AddVertex(PointType{ -snap_distance, 0.0, 0.0 });
    polygon.AddVertex(PointType{ snap_distance, 0.0, 0.0 });
    polygon.AddVertex(PointType{ 1.0, 1.0, 0.0 });
    polygon.AddVertex(PointType{ 1.0, 0.0, 0.0 });

    const auto split = embedding::detail::SplitByPlane(polygon, 0, 0.0, tolerance);
    QuESo_CHECK(!split.negative.has_value());
    BOOST_REQUIRE(split.positive.has_value());
    QuESo_CHECK_EQUAL(split.positive->NumberOfVertices(), IndexType{ 3 });
    bool found_snapped_vertex = false;
    for (PointView r_vertex : split.positive->Vertices()) {
        if (r_vertex[0] == 0.0) {
            QuESo_CHECK_EQUAL(r_vertex[1], 0.0);
            QuESo_CHECK_EQUAL(r_vertex[2], 0.0);
            found_snapped_vertex = true;
        }
    }
    QuESo_CHECK(found_snapped_vertex);
    for (IndexType i = 0; i < split.positive->NumberOfVertices(); ++i) {
        const IndexType next = (i + 1) % split.positive->NumberOfVertices();
        QuESo_CHECK(
            Math::SquaredNorm(split.positive->Vertex(i) - split.positive->Vertex(next)) > snap_distance * snap_distance
        );
    }
}

BOOST_AUTO_TEST_CASE(ConvexPolygonDoesNotEmitSideWithoutStrictVertex)
{
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1e11, .coordinate_scale = 1.0 });
    embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 1.0, 0.0 });
    polygon.AddVertex(PointType{ 0.0, 0.05, 0.0 });
    polygon.AddVertex(PointType{ 1.0, 0.06, 0.0 });
    polygon.AddVertex(PointType{ 1.0, 0.2, 14.0 });
    polygon.AddVertex(PointType{ 0.0, 0.1, 5.0 });

    const auto split = embedding::detail::SplitByPlane(polygon, 1, 0.0, tolerance);
    QuESo_CHECK(!split.negative.has_value());
    BOOST_REQUIRE(split.positive.has_value());
    for (PointView r_vertex : split.positive->Vertices()) { QuESo_CHECK(r_vertex[1] >= 0.0); }
}

BOOST_AUTO_TEST_CASE(ConvexPolygonRetainsThinGeometryAtLegacyScale)
{
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
    polygon.AddVertex(PointType{ 0.0, 0.0, 0.0 });
    polygon.AddVertex(PointType{ 1e-7, 0.0, 0.0 });
    polygon.AddVertex(PointType{ 0.0, 1e-7, 0.0 });

    const auto split = embedding::detail::SplitByPlane(polygon, 0, 0.0, tolerance);
    QuESo_CHECK(!split.negative.has_value());
    QuESo_CHECK(split.positive.has_value());
}

BOOST_AUTO_TEST_CASE(LongTriangleIsPartitionedWithoutAreaLoss)
{
    const GridIndexer grid_indexer = MakeGridIndexer();
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 0.25, 0.2, 0.5 }, PointType{ 1.75, 0.2, 0.5 }, PointType{ 1.0, 0.8, 0.5 });
    const double reference_area = MeshUtilities::Area(mesh.View());

    DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
    QuESo_CHECK_EQUAL(partitioner.GetCellIndices().size(), IndexType{ 2 });
    QuESo_CHECK(partitioner.GetFaceIds().empty());
    double area{};
    for (const IndexType cell : partitioner.GetCellIndices()) {
        area += MeshUtilities::Area(partitioner.TakeCellProduct(cell).SurfaceView());
    }
    QuESo_CHECK_LT(std::abs(area - reference_area), 1e-12);
}

BOOST_AUTO_TEST_CASE(GridAlignedTriangleProducesOneFaceSection)
{
    const GridIndexer grid_indexer = MakeGridIndexer();
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 1.0, 0.2, 0.2 }, PointType{ 1.0, 0.8, 0.2 }, PointType{ 1.0, 0.2, 0.8 });

    DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
    QuESo_CHECK(partitioner.GetCellIndices().empty());
    BOOST_REQUIRE_EQUAL(partitioner.GetFaceIds().size(), IndexType{ 1 });
    const GridFaceId expected = grid_indexer.GetFace(0, GridIndexer::Direction::x_forward);
    QuESo_CHECK(partitioner.GetFaceIds().front() == expected);
    const auto section = partitioner.TakeGridFaceSurface(expected);
    QuESo_CHECK_EQUAL(section.NumOfTriangles(), IndexType{ 1 });
}

BOOST_AUTO_TEST_CASE(PartitioningDependsOnCellScaleNotWholeGridExtent)
{
    const GridIndexer two_cells = MakeGridIndexer({ 2, 1, 1 });
    const GridIndexer ten_cells = MakeGridIndexer({ 10, 1, 1 });
    const double offset = 5.0 * two_cells.GetGeometryTolerance().SnapDistance();
    TriangleMesh mesh;
    AddTriangle(
        mesh,
        PointType{ 1.0 + offset, 0.2, 0.2 },
        PointType{ 1.0 + offset, 0.8, 0.2 },
        PointType{ 1.0 + offset, 0.2, 0.8 }
    );

    DomainMeshPartitioner short_partitioner(mesh.View(), two_cells);
    DomainMeshPartitioner long_partitioner(mesh.View(), ten_cells);
    QuESo_CHECK(short_partitioner.GetFaceIds().empty());
    QuESo_CHECK(long_partitioner.GetFaceIds().empty());
    BOOST_REQUIRE_EQUAL(short_partitioner.GetCellIndices().size(), IndexType{ 1 });
    BOOST_REQUIRE_EQUAL(long_partitioner.GetCellIndices().size(), IndexType{ 1 });
    QuESo_CHECK_EQUAL(short_partitioner.GetCellIndices().front(), IndexType{ 1 });
    QuESo_CHECK_EQUAL(long_partitioner.GetCellIndices().front(), IndexType{ 1 });
    const double short_area = MeshUtilities::Area(short_partitioner.TakeCellProduct(1).SurfaceView());
    const double long_area = MeshUtilities::Area(long_partitioner.TakeCellProduct(1).SurfaceView());
    QuESo_CHECK_NEAR(short_area, long_area, EPS4);
}

BOOST_AUTO_TEST_CASE(PartitioningUsesCoordinateResolutionAfterTranslation)
{
    const BoundingBoxType bounds = MakeBox({ 1e12, 1e12, 1e12 }, { 1e12 + 2.0, 1e12 + 1.0, 1e12 + 1.0 });
    const GridIndexer grid_indexer = MakeGridIndexer(bounds, { 2, 1, 1 });
    const double plane = 1e12 + 1.0;
    const double offset = 0.5 * grid_indexer.GetGeometryTolerance().SnapDistance();
    TriangleMesh mesh;
    AddTriangle(
        mesh,
        PointType{ plane + offset, 1e12 + 0.2, 1e12 + 0.2 },
        PointType{ plane + offset, 1e12 + 0.8, 1e12 + 0.2 },
        PointType{ plane + offset, 1e12 + 0.2, 1e12 + 0.8 }
    );

    DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
    QuESo_CHECK(partitioner.GetCellIndices().empty());
    BOOST_REQUIRE_EQUAL(partitioner.GetFaceIds().size(), IndexType{ 1 });
    const GridFaceId expected = grid_indexer.GetFace(0, GridIndexer::Direction::x_forward);
    QuESo_CHECK(partitioner.GetFaceIds().front() == expected);
}

BOOST_AUTO_TEST_CASE(ExplicitSegmentsUseOutwardFaceOrientation)
{
    const GridIndexer grid_indexer = MakeGridIndexer({ 1, 1, 1 });
    for (IndexType face_index = 0; face_index < 6; ++face_index) {
        const embedding::CellFace face = embedding::CellFace::FromIndex(face_index);
        const IndexType first_in_plane = (face.axis + 1) % 3;
        const IndexType second_in_plane = (face.axis + 2) % 3;
        const double face_coordinate = face.is_upper ? 1.0 : 0.0;
        const double outside_coordinate = face.is_upper ? 1.5 : -0.5;
        TriangleMesh mesh;
        PointType outside{ 0.2, 0.2, 0.2 };
        PointType inside_first{ 0.2, 0.2, 0.2 };
        PointType inside_second{ 0.2, 0.2, 0.2 };
        outside[face.axis] = outside_coordinate;
        inside_first[face.axis] = 0.5;
        inside_second[face.axis] = 0.5;
        inside_second[first_in_plane] = 0.8;
        outside[second_in_plane] = 0.3;
        inside_first[second_in_plane] = 0.3;
        inside_second[second_in_plane] = 0.3;
        AddTriangle(mesh, outside, inside_first, inside_second);

        DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
        BOOST_REQUIRE_EQUAL(partitioner.GetCellIndices().size(), IndexType{ 1 });
        const auto section = partitioner.TakeCellProduct(0);
        const auto segments = section.FaceSegments(face);
        BOOST_REQUIRE_EQUAL(segments.size(), IndexType{ 1 });
        const auto& r_segment = segments.front();
        QuESo_CHECK_EQUAL(r_segment.first[face.axis], face_coordinate);
        QuESo_CHECK_EQUAL(r_segment.second[face.axis], face_coordinate);
        QuESo_CHECK_EQUAL(r_segment.source_triangle, IndexType{ 0 });

        PointType face_normal{};
        face_normal[face.axis] = face.is_upper ? 1.0 : -1.0;
        const Vector3d tangent = r_segment.second - r_segment.first;
        const Vector3d left = Math::Cross(face_normal, tangent);
        const Vector3d projected_source_normal =
            r_segment.source_normal - Math::Dot(r_segment.source_normal, face_normal) * face_normal;
        QuESo_CHECK_GT(Math::Dot(left, -1.0 * projected_source_normal), 0.0);
    }
}

BOOST_AUTO_TEST_CASE(CellEdgeTraceIsRecordedForBothIncidentFaces)
{
    const GridIndexer grid_indexer = MakeGridIndexer({ 1, 1, 1 });
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 0.0, 0.0, 0.2 }, PointType{ 0.0, 0.0, 0.8 }, PointType{ 0.5, 0.5, 0.5 });

    DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
    const auto section = partitioner.TakeCellProduct(0);
    const auto first_face = section.FaceSegments({ 0, false });
    const auto second_face = section.FaceSegments({ 1, false });
    BOOST_REQUIRE_EQUAL(first_face.size(), IndexType{ 1 });
    BOOST_REQUIRE_EQUAL(second_face.size(), IndexType{ 1 });
    const bool same_direction =
        section.GetGeometryTolerance().PointsAreSame(first_face[0].first, second_face[0].first)
        && section.GetGeometryTolerance().PointsAreSame(first_face[0].second, second_face[0].second);
    const bool opposite_direction =
        section.GetGeometryTolerance().PointsAreSame(first_face[0].first, second_face[0].second)
        && section.GetGeometryTolerance().PointsAreSame(first_face[0].second, second_face[0].first);
    QuESo_CHECK(same_direction || opposite_direction);
}

BOOST_AUTO_TEST_CASE(ExplicitSegmentGeometryIsIndependentOfSourceOrder)
{
    const GridIndexer grid_indexer = MakeGridIndexer({ 1, 1, 1 });
    const std::array first_triangle{ PointType{ 0.0, 0.1, 0.2 },
                                     PointType{ 0.0, 0.3, 0.2 },
                                     PointType{ 0.5, 0.2, 0.5 } };
    const std::array second_triangle{ PointType{ 0.0, 0.6, 0.2 },
                                      PointType{ 0.0, 0.8, 0.2 },
                                      PointType{ 0.5, 0.7, 0.5 } };
    TriangleMesh forward_mesh;
    TriangleMesh reverse_mesh;
    AddTriangle(forward_mesh, first_triangle[0], first_triangle[1], first_triangle[2]);
    AddTriangle(forward_mesh, second_triangle[0], second_triangle[1], second_triangle[2]);
    AddTriangle(reverse_mesh, second_triangle[0], second_triangle[1], second_triangle[2]);
    AddTriangle(reverse_mesh, first_triangle[0], first_triangle[1], first_triangle[2]);

    DomainMeshPartitioner forward_partitioner(forward_mesh.View(), grid_indexer);
    DomainMeshPartitioner reverse_partitioner(reverse_mesh.View(), grid_indexer);
    const auto forward_section = forward_partitioner.TakeCellProduct(0);
    const auto reverse_section = reverse_partitioner.TakeCellProduct(0);
    const auto forward = forward_section.FaceSegments({ 0, false });
    const auto reverse = reverse_section.FaceSegments({ 0, false });
    BOOST_REQUIRE_EQUAL(forward.size(), reverse.size());
    for (IndexType i = 0; i < forward.size(); ++i) {
        QuESo_CHECK_POINT_NEAR(forward[i].first, reverse[i].first, 0.0);
        QuESo_CHECK_POINT_NEAR(forward[i].second, reverse[i].second, 0.0);
        QuESo_CHECK_POINT_NEAR(forward[i].source_normal, reverse[i].source_normal, 0.0);
        QuESo_CHECK_EQUAL(forward[i].source_triangle, IndexType{ 1 } - reverse[i].source_triangle);
    }
}

BOOST_AUTO_TEST_CASE(DomainAndSurfaceProductsMatchSurfaceGeometry)
{
    const GridIndexer grid_indexer = MakeGridIndexer();
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 0.25, 0.2, 0.5 }, PointType{ 1.75, 0.2, 0.5 }, PointType{ 1.0, 0.8, 0.5 });

    DomainMeshPartitioner domain_partitioner(mesh.View(), grid_indexer);
    SurfaceMeshPartitioner surface_partitioner(mesh.View(), grid_indexer);
    for (const IndexType cell : domain_partitioner.GetCellIndices()) {
        const auto explicit_section = domain_partitioner.TakeCellProduct(cell);
        const auto surface = surface_partitioner.TakeCellProduct(cell);
        QuESo_CHECK_NEAR(
            MeshUtilities::Area(explicit_section.SurfaceView()), MeshUtilities::Area(surface.View()), EPS4
        );
    }
}

BOOST_AUTO_TEST_CASE(GridAlignedTriangleProducesExteriorFaceSection)
{
    const GridIndexer grid_indexer = MakeGridIndexer();
    for (const bool IsUpper : { false, true }) {
        const double coordinate = IsUpper ? 2.0 : 0.0;
        TriangleMesh mesh;
        AddTriangle(
            mesh,
            PointType{ coordinate, 0.2, 0.2 },
            PointType{ coordinate, 0.8, 0.2 },
            PointType{ coordinate, 0.2, 0.8 }
        );

        DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
        QuESo_CHECK(partitioner.GetCellIndices().empty());
        BOOST_REQUIRE_EQUAL(partitioner.GetFaceIds().size(), IndexType{ 1 });
        const GridIndexer::Direction direction =
            IsUpper ? GridIndexer::Direction::x_forward : GridIndexer::Direction::x_backward;
        const IndexType cell = IsUpper ? 1 : 0;
        const GridFaceId expected = grid_indexer.GetFace(cell, direction);
        QuESo_CHECK(partitioner.GetFaceIds().front() == expected);
        QuESo_CHECK_LT(
            std::abs(
                MeshUtilities::Area(partitioner.TakeGridFaceSurface(expected).View()) - MeshUtilities::Area(mesh.View())
            ),
            1e-12
        );
    }
}

BOOST_AUTO_TEST_CASE(GridAlignedTriangleProducesFacePatchesWithoutAreaLoss)
{
    const GridIndexer grid_indexer = MakeGridIndexer(Vector3i{ 3, 3, 3 });
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 1.0, 0.2, 0.2 }, PointType{ 1.0, 2.8, 0.2 }, PointType{ 1.0, 0.2, 2.8 });

    DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
    QuESo_CHECK(partitioner.GetCellIndices().empty());
    BOOST_REQUIRE_EQUAL(partitioner.GetFaceIds().size(), IndexType{ 6 });
    const std::array<GridFaceId, 6> expected_faces{
        grid_indexer.GetFace(grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 0), GridIndexer::Direction::x_backward),
        grid_indexer.GetFace(grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 1), GridIndexer::Direction::x_backward),
        grid_indexer.GetFace(grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 2), GridIndexer::Direction::x_backward),
        grid_indexer.GetFace(grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 0), GridIndexer::Direction::x_backward),
        grid_indexer.GetFace(grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1), GridIndexer::Direction::x_backward),
        grid_indexer.GetFace(grid_indexer.GetVectorIndexFromMatrixIndices(1, 2, 0), GridIndexer::Direction::x_backward)
    };
    for (const GridFaceId expected : expected_faces) {
        bool is_present = false;
        for (const GridFaceId face : partitioner.GetFaceIds()) { is_present = is_present || face == expected; }
        QuESo_CHECK(is_present);
    }
    double partitioned_area{};
    for (const GridFaceId face : partitioner.GetFaceIds()) {
        QuESo_CHECK_EQUAL(grid_indexer.GetFaceAxis(face), IndexType{ 0 });
        partitioned_area += MeshUtilities::Area(partitioner.TakeGridFaceSurface(face).View());
    }
    QuESo_CHECK_LT(std::abs(partitioned_area - MeshUtilities::Area(mesh.View())), 1e-12);
}

BOOST_AUTO_TEST_CASE(ClipsToExactGridOverlapAndIgnoresExternalGeometry)
{
    const GridIndexer grid_indexer = MakeGridIndexer();
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ -1.0, 0.2, 0.5 }, PointType{ 1.0, 0.2, 0.5 }, PointType{ 1.0, 0.8, 0.5 });
    AddTriangle(mesh, PointType{ -2.0, 0.2, 0.5 }, PointType{ -1.5, 0.2, 0.5 }, PointType{ -1.5, 0.8, 0.5 });

    DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
    double retained_area{};
    for (const IndexType cell : partitioner.GetCellIndices()) {
        retained_area += MeshUtilities::Area(partitioner.TakeCellProduct(cell).SurfaceView());
    }
    QuESo_CHECK_LT(std::abs(retained_area - 0.45), 1e-12);
    QuESo_CHECK(partitioner.GetFaceIds().empty());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
