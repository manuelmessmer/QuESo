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
#include <doctest/doctest.h>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/embedding/convex_polygon.h"
#include "queso/embedding/mesh_partitioner.h"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::testing {
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

TEST_SUITE("MeshPartitionerTestSuite")
{

    TEST_CASE("MeshPartitionerTestSuite - Product Consumption Preconditions Are Checked In Debug")
    {
        if constexpr (NOTDEBUG) { return; }
        const auto grid = MakeGridIndexer({ 2, 1, 1 });
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.2, 0.2, 0.5 }, PointType{ 0.8, 0.2, 0.5 }, PointType{ 0.2, 0.8, 0.5 });
        SurfaceMeshPartitioner partitioner(mesh.View(), grid);
        CHECK_THROWS_AS((void)partitioner.TakeCellProduct(2), Exception);
        CHECK_THROWS_AS((void)partitioner.TakeCellProduct(1), Exception);
        CHECK_NOTHROW((void)partitioner.TakeCellProduct(0));
        CHECK_THROWS_AS((void)partitioner.TakeCellProduct(0), Exception);
        CHECK_THROWS_AS((void)partitioner.TakeGridFaceSurface(GridFaceId{ grid.NumberOfFaces() }), Exception);
        CHECK_THROWS_AS(
            (void)partitioner.TakeGridFaceSurface(grid.GetFace(0, GridIndexer::Direction::x_forward)), Exception
        );
    }

    TEST_CASE("MeshPartitionerTestSuite - Convex Polygon Split Conserves Area")
    {
        embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
        polygon.AddVertex(PointType{ 0.0, 0.0, 0.5 });
        polygon.AddVertex(PointType{ 2.0, 0.0, 0.5 });
        polygon.AddVertex(PointType{ 1.0, 1.0, 0.5 });

        const auto split = embedding::detail::SplitByPlane(
            polygon, 0, 1.0, GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
        );
        REQUIRE(split.negative.has_value());
        REQUIRE(split.positive.has_value());
        CHECK((std::abs(split.negative->Area() + split.positive->Area() - polygon.Area())) < (1e-12));
        for (const auto& r_piece : { *split.negative, *split.positive }) {
            CHECK((r_piece.NumberOfVertices()) == (IndexType{ 3 }));
            CheckPointNear(r_piece.Normal(), polygon.Normal(), 1e-12);
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Convex Polygon On Plane Produces Both Sides")
    {
        embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
        polygon.AddVertex(PointType{ 1.0, 0.0, 0.0 });
        polygon.AddVertex(PointType{ 1.0, 1.0, 0.0 });
        polygon.AddVertex(PointType{ 1.0, 0.0, 1.0 });

        const auto split = embedding::detail::SplitByPlane(
            polygon, 0, 1.0, GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
        );
        REQUIRE(split.negative.has_value());
        REQUIRE(split.positive.has_value());
        for (const auto& r_piece : { *split.negative, *split.positive }) {
            CHECK((std::abs(r_piece.Area() - polygon.Area())) < (1e-12));
            CHECK((r_piece.NumberOfVertices()) == (polygon.NumberOfVertices()));
            CheckPointNear(r_piece.Normal(), polygon.Normal(), 1e-12);
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Convex Polygon Snaps And Deduplicates Vertices")
    {
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1e11, .coordinate_scale = 1.0 });
        const double snap_distance = tolerance.SnapDistance();
        embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
        polygon.AddVertex(PointType{ -snap_distance, 0.0, 0.0 });
        polygon.AddVertex(PointType{ snap_distance, 0.0, 0.0 });
        polygon.AddVertex(PointType{ 1.0, 1.0, 0.0 });
        polygon.AddVertex(PointType{ 1.0, 0.0, 0.0 });

        const auto split = embedding::detail::SplitByPlane(polygon, 0, 0.0, tolerance);
        REQUIRE(!split.negative.has_value());
        REQUIRE(split.positive.has_value());
        CHECK((split.positive->NumberOfVertices()) == (IndexType{ 3 }));
        bool found_snapped_vertex = false;
        for (PointView r_vertex : split.positive->Vertices()) {
            if (r_vertex[0] == 0.0) {
                CHECK((r_vertex[1]) == (0.0));
                CHECK((r_vertex[2]) == (0.0));
                found_snapped_vertex = true;
            }
        }
        CHECK(found_snapped_vertex);
        for (IndexType i = 0; i < split.positive->NumberOfVertices(); ++i) {
            const IndexType next = (i + 1) % split.positive->NumberOfVertices();
            CHECK(
                Math::SquaredNorm(split.positive->Vertex(i) - split.positive->Vertex(next))
                > snap_distance * snap_distance
            );
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Convex Polygon Does Not Emit Side Without Strict Vertex")
    {
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1e11, .coordinate_scale = 1.0 });
        embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 1.0, 0.0 });
        polygon.AddVertex(PointType{ 0.0, 0.05, 0.0 });
        polygon.AddVertex(PointType{ 1.0, 0.06, 0.0 });
        polygon.AddVertex(PointType{ 1.0, 0.2, 14.0 });
        polygon.AddVertex(PointType{ 0.0, 0.1, 5.0 });

        const auto split = embedding::detail::SplitByPlane(polygon, 1, 0.0, tolerance);
        REQUIRE(!split.negative.has_value());
        REQUIRE(split.positive.has_value());
        for (PointView r_vertex : split.positive->Vertices()) { CHECK(r_vertex[1] >= 0.0); }
    }

    TEST_CASE("MeshPartitionerTestSuite - Convex Polygon Retains Thin Geometry At Legacy Scale")
    {
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        embedding::detail::ConvexPolygon polygon(PointType{ 0.0, 0.0, 1.0 });
        polygon.AddVertex(PointType{ 0.0, 0.0, 0.0 });
        polygon.AddVertex(PointType{ 1e-7, 0.0, 0.0 });
        polygon.AddVertex(PointType{ 0.0, 1e-7, 0.0 });

        const auto split = embedding::detail::SplitByPlane(polygon, 0, 0.0, tolerance);
        REQUIRE(!split.negative.has_value());
        REQUIRE(split.positive.has_value());
    }

    TEST_CASE("MeshPartitionerTestSuite - Long Triangle Is Partitioned Without Area Loss")
    {
        const GridIndexer grid_indexer = MakeGridIndexer();
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.25, 0.2, 0.5 }, PointType{ 1.75, 0.2, 0.5 }, PointType{ 1.0, 0.8, 0.5 });
        const double reference_area = MeshUtilities::Area(mesh.View());

        DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
        REQUIRE((partitioner.GetCellIndices().size()) == (IndexType{ 2 }));
        REQUIRE(partitioner.GetFaceIds().empty());
        double area{};
        for (const IndexType cell : partitioner.GetCellIndices()) {
            area += MeshUtilities::Area(partitioner.TakeCellProduct(cell).SurfaceView());
        }
        CHECK((std::abs(area - reference_area)) < (1e-12));
    }

    TEST_CASE("MeshPartitionerTestSuite - Grid Aligned Triangle Produces One Face Section")
    {
        const GridIndexer grid_indexer = MakeGridIndexer();
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 1.0, 0.2, 0.2 }, PointType{ 1.0, 0.8, 0.2 }, PointType{ 1.0, 0.2, 0.8 });

        DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
        REQUIRE(partitioner.GetCellIndices().empty());
        REQUIRE((partitioner.GetFaceIds().size()) == (IndexType{ 1 }));
        const GridFaceId expected = grid_indexer.GetFace(0, GridIndexer::Direction::x_forward);
        CHECK(partitioner.GetFaceIds().front() == expected);
        const auto section = partitioner.TakeGridFaceSurface(expected);
        CHECK((section.NumOfTriangles()) == (IndexType{ 1 }));
    }

    TEST_CASE("MeshPartitionerTestSuite - Partitioning Depends On Cell Scale Not Whole Grid Extent")
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
        REQUIRE(short_partitioner.GetFaceIds().empty());
        REQUIRE(long_partitioner.GetFaceIds().empty());
        REQUIRE((short_partitioner.GetCellIndices().size()) == (IndexType{ 1 }));
        REQUIRE((long_partitioner.GetCellIndices().size()) == (IndexType{ 1 }));
        CHECK((short_partitioner.GetCellIndices().front()) == (IndexType{ 1 }));
        CHECK((long_partitioner.GetCellIndices().front()) == (IndexType{ 1 }));
        const double short_area = MeshUtilities::Area(short_partitioner.TakeCellProduct(1).SurfaceView());
        const double long_area = MeshUtilities::Area(long_partitioner.TakeCellProduct(1).SurfaceView());
        CheckNear(short_area, long_area, EPS4);
    }

    TEST_CASE("MeshPartitionerTestSuite - Partitioning Uses Coordinate Resolution After Translation")
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
        REQUIRE(partitioner.GetCellIndices().empty());
        REQUIRE((partitioner.GetFaceIds().size()) == (IndexType{ 1 }));
        const GridFaceId expected = grid_indexer.GetFace(0, GridIndexer::Direction::x_forward);
        CHECK(partitioner.GetFaceIds().front() == expected);
    }

    TEST_CASE("MeshPartitionerTestSuite - Explicit Segments Use Outward Face Orientation")
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
            REQUIRE((partitioner.GetCellIndices().size()) == (IndexType{ 1 }));
            const auto section = partitioner.TakeCellProduct(0);
            const auto segments = section.FaceSegments(face);
            REQUIRE((segments.size()) == (IndexType{ 1 }));
            const auto& r_segment = segments.front();
            CHECK((r_segment.first[face.axis]) == (face_coordinate));
            CHECK((r_segment.second[face.axis]) == (face_coordinate));
            CHECK((r_segment.source_triangle) == (IndexType{ 0 }));

            PointType face_normal{};
            face_normal[face.axis] = face.is_upper ? 1.0 : -1.0;
            const Vector3d tangent = r_segment.second - r_segment.first;
            const Vector3d left = Math::Cross(face_normal, tangent);
            const Vector3d projected_source_normal =
                r_segment.source_normal - Math::Dot(r_segment.source_normal, face_normal) * face_normal;
            CHECK((Math::Dot(left, -1.0 * projected_source_normal)) > (0.0));
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Cell Edge Trace Is Recorded For Both Incident Faces")
    {
        const GridIndexer grid_indexer = MakeGridIndexer({ 1, 1, 1 });
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.0, 0.0, 0.2 }, PointType{ 0.0, 0.0, 0.8 }, PointType{ 0.5, 0.5, 0.5 });

        DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
        const auto section = partitioner.TakeCellProduct(0);
        const auto first_face = section.FaceSegments({ 0, false });
        const auto second_face = section.FaceSegments({ 1, false });
        REQUIRE((first_face.size()) == (IndexType{ 1 }));
        REQUIRE((second_face.size()) == (IndexType{ 1 }));
        const bool same_direction =
            section.GetGeometryTolerance().PointsAreSame(first_face[0].first, second_face[0].first)
            && section.GetGeometryTolerance().PointsAreSame(first_face[0].second, second_face[0].second);
        const bool opposite_direction =
            section.GetGeometryTolerance().PointsAreSame(first_face[0].first, second_face[0].second)
            && section.GetGeometryTolerance().PointsAreSame(first_face[0].second, second_face[0].first);
        CHECK((same_direction || opposite_direction));
    }

    TEST_CASE("MeshPartitionerTestSuite - Explicit Segment Geometry Is Independent Of Source Order")
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
        REQUIRE((forward.size()) == (reverse.size()));
        for (IndexType i = 0; i < forward.size(); ++i) {
            CheckPointNear(forward[i].first, reverse[i].first, 0.0);
            CheckPointNear(forward[i].second, reverse[i].second, 0.0);
            CheckPointNear(forward[i].source_normal, reverse[i].source_normal, 0.0);
            CHECK((forward[i].source_triangle) == (IndexType{ 1 } - reverse[i].source_triangle));
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Domain And Surface Products Match Surface Geometry")
    {
        const GridIndexer grid_indexer = MakeGridIndexer();
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.25, 0.2, 0.5 }, PointType{ 1.75, 0.2, 0.5 }, PointType{ 1.0, 0.8, 0.5 });

        DomainMeshPartitioner domain_partitioner(mesh.View(), grid_indexer);
        SurfaceMeshPartitioner surface_partitioner(mesh.View(), grid_indexer);
        for (const IndexType cell : domain_partitioner.GetCellIndices()) {
            const auto explicit_section = domain_partitioner.TakeCellProduct(cell);
            const auto surface = surface_partitioner.TakeCellProduct(cell);
        CheckNear(MeshUtilities::Area(explicit_section.SurfaceView()), MeshUtilities::Area(surface.View()), EPS4);
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Grid Aligned Triangle Produces Exterior Face Section")
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
            REQUIRE(partitioner.GetCellIndices().empty());
            REQUIRE((partitioner.GetFaceIds().size()) == (IndexType{ 1 }));
            const GridIndexer::Direction direction =
                IsUpper ? GridIndexer::Direction::x_forward : GridIndexer::Direction::x_backward;
            const IndexType cell = IsUpper ? 1 : 0;
            const GridFaceId expected = grid_indexer.GetFace(cell, direction);
            CHECK(partitioner.GetFaceIds().front() == expected);
            CHECK(
                (std::abs(
                    MeshUtilities::Area(partitioner.TakeGridFaceSurface(expected).View())
                    - MeshUtilities::Area(mesh.View())
                ))
                < (1e-12)
            );
        }
    }

    TEST_CASE("MeshPartitionerTestSuite - Grid Aligned Triangle Produces Face Patches Without Area Loss")
    {
        const GridIndexer grid_indexer = MakeGridIndexer(Vector3i{ 3, 3, 3 });
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 1.0, 0.2, 0.2 }, PointType{ 1.0, 2.8, 0.2 }, PointType{ 1.0, 0.2, 2.8 });

        DomainMeshPartitioner partitioner(mesh.View(), grid_indexer);
        REQUIRE(partitioner.GetCellIndices().empty());
        REQUIRE((partitioner.GetFaceIds().size()) == (IndexType{ 6 }));
        const std::array<GridFaceId, 6> expected_faces{
            grid_indexer.GetFace(
                grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 0), GridIndexer::Direction::x_backward
            ),
            grid_indexer.GetFace(
                grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 1), GridIndexer::Direction::x_backward
            ),
            grid_indexer.GetFace(
                grid_indexer.GetVectorIndexFromMatrixIndices(1, 0, 2), GridIndexer::Direction::x_backward
            ),
            grid_indexer.GetFace(
                grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 0), GridIndexer::Direction::x_backward
            ),
            grid_indexer.GetFace(
                grid_indexer.GetVectorIndexFromMatrixIndices(1, 1, 1), GridIndexer::Direction::x_backward
            ),
            grid_indexer.GetFace(
                grid_indexer.GetVectorIndexFromMatrixIndices(1, 2, 0), GridIndexer::Direction::x_backward
            )
        };
        for (const GridFaceId expected : expected_faces) {
            bool is_present = false;
            for (const GridFaceId face : partitioner.GetFaceIds()) { is_present = is_present || face == expected; }
            CHECK(is_present);
        }
        double partitioned_area{};
        for (const GridFaceId face : partitioner.GetFaceIds()) {
            CHECK((grid_indexer.GetFaceAxis(face)) == (IndexType{ 0 }));
            partitioned_area += MeshUtilities::Area(partitioner.TakeGridFaceSurface(face).View());
        }
        CHECK((std::abs(partitioned_area - MeshUtilities::Area(mesh.View()))) < (1e-12));
    }

    TEST_CASE("MeshPartitionerTestSuite - Clips To Exact Grid Overlap And Ignores External Geometry")
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
        CHECK((std::abs(retained_area - 0.45)) < (1e-12));
        REQUIRE(partitioner.GetFaceIds().empty());
    }

}  // TEST_SUITE

}  // namespace queso::testing
