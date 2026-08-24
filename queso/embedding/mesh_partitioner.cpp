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
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>

//// Project includes
#include "queso/embedding/convex_polygon.h"
#include "queso/embedding/mesh_partitioner.h"

namespace queso::embedding {
namespace {

    using detail::ConvexPolygon;

    // Stores one clipped triangle piece with stable ordering keys. Per-source-triangle storage avoids synchronization
    // while source triangles are processed in parallel.
    struct PolygonRecord
    {
        // Cell index or GridFaceId value of the destination section.
        IndexType target;
        // Index of the input triangle that produced this polygon.
        IndexType source_triangle;
        // Sequential index of this polygon among pieces from the same source triangle.
        IndexType piece;
        ConvexPolygon polygon;
    };

    struct PendingCellSurfaceSection
    {
        TriangleMesh surface;
        CellFaceContours contours;
    };

    // Converts a source triangle to the mutable convex representation used by the clipping operations below.
    [[nodiscard]] ConvexPolygon MakePolygon(TriangleProxy<WithNormals> rTriangle)
    {
        ConvexPolygon polygon(rTriangle.Normal);
        polygon.AddVertex(rTriangle.P1);
        polygon.AddVertex(rTriangle.P2);
        polygon.AddVertex(rTriangle.P3);
        return polygon;
    }

    void SortFaceSegments(CellFaceContours& rContours)
    {
        for (auto& r_segments : rContours.segments) {
            std::sort(r_segments.begin(), r_segments.end(), [](const auto& rLeft, const auto& rRight) {
                return std::tie(rLeft.first, rLeft.second, rLeft.source_normal, rLeft.source_triangle)
                       < std::tie(rRight.first, rRight.second, rRight.source_normal, rRight.source_triangle);
            });
        }
    }

    // Retains the intersection with the closed grid AABB by clipping against its lower and upper plane on each axis.
    [[nodiscard]] std::optional<ConvexPolygon>
        ClipToGrid(ConvexPolygon Polygon, const BoundingBoxType& rBounds, const GeometryTolerance& rTolerance)
    {
        for (IndexType axis = 0; axis < 3; ++axis) {
            auto lower = detail::SplitByPlane(Polygon, axis, rBounds.lower[axis], rTolerance);
            if (!lower.positive) { return std::nullopt; }
            auto upper = detail::SplitByPlane(*lower.positive, axis, rBounds.upper[axis], rTolerance);
            if (!upper.negative) { return std::nullopt; }
            Polygon = std::move(*upper.negative);
        }
        return Polygon;
    }

    // Splits retained polygons at internal grid planes intersecting their extent along Axis.
    void SplitAtInternalPlanes(
        std::vector<ConvexPolygon>& rPolygons,
        const GridIndexer& rGridIndexer,
        IndexType Axis,
        const GeometryTolerance& rTolerance
    )
    {
        const double snap_distance = rTolerance.SnapDistance();
        const auto& r_bounds = rGridIndexer.BoundsXYZ();
        const Vector3i counts = rGridIndexer.ElementCounts();
        const double delta = (r_bounds.upper[Axis] - r_bounds.lower[Axis]) / static_cast<double>(counts[Axis]);
        double minimum = MAXD;
        double maximum = LOWESTD;
        for (const auto& r_polygon : rPolygons) {
            for (PointView r_vertex : r_polygon.Vertices()) {
                minimum = std::min(minimum, r_vertex[Axis]);
                maximum = std::max(maximum, r_vertex[Axis]);
            }
        }
        // Plane positions may fall outside the grid before clamping to the valid internal-plane range.
        const auto first_plane =
            static_cast<std::ptrdiff_t>(std::ceil((minimum - snap_distance - r_bounds.lower[Axis]) / delta));
        const auto last_plane =
            static_cast<std::ptrdiff_t>(std::floor((maximum + snap_distance - r_bounds.lower[Axis]) / delta));
        const auto first_internal_plane = std::max<std::ptrdiff_t>(1, first_plane);
        const auto last_internal_plane =
            std::min<std::ptrdiff_t>(static_cast<std::ptrdiff_t>(counts[Axis]) - 1, last_plane);
        for (std::ptrdiff_t raw_plane = first_internal_plane; raw_plane <= last_internal_plane; ++raw_plane) {
            const IndexType plane = static_cast<IndexType>(raw_plane);
            const double coordinate = r_bounds.lower[Axis] + delta * static_cast<double>(plane);
            std::vector<ConvexPolygon> split_polygons;
            split_polygons.reserve(rPolygons.size() + 1);
            for (const auto& r_polygon : rPolygons) {
                double polygon_minimum = MAXD;
                double polygon_maximum = LOWESTD;
                for (PointView r_vertex : r_polygon.Vertices()) {
                    polygon_minimum = std::min(polygon_minimum, r_vertex[Axis]);
                    polygon_maximum = std::max(polygon_maximum, r_vertex[Axis]);
                }
                if (polygon_minimum >= coordinate - snap_distance && polygon_maximum <= coordinate + snap_distance) {
                    // Keep one canonical copy of a polygon lying on a shared grid plane.
                    auto on_plane = detail::SplitByPlane(r_polygon, Axis, coordinate, rTolerance);
                    QuESo_ASSERT(on_plane.positive.has_value(), "An on-plane polygon must remain non-empty.");
                    split_polygons.push_back(std::move(*on_plane.positive));
                    continue;
                }
                if (polygon_maximum < coordinate - snap_distance || polygon_minimum > coordinate + snap_distance) {
                    // This plane does not intersect the polygon, so avoid an unnecessary split.
                    split_polygons.push_back(r_polygon);
                    continue;
                }
                // Retain each non-degenerate side of a polygon crossing the plane.
                auto split = detail::SplitByPlane(r_polygon, Axis, coordinate, rTolerance);
                if (split.negative) { split_polygons.push_back(std::move(*split.negative)); }
                if (split.positive) { split_polygons.push_back(std::move(*split.positive)); }
            }
            rPolygons = std::move(split_polygons);
        }
    }

    // Returns the canonical grid face for a polygon lying completely on a grid plane, if any.
    [[nodiscard]] std::optional<GridFaceId> GetFaceTarget(
        const ConvexPolygon& rPolygon,
        const GridIndexer& rGridIndexer,
        const GeometryTolerance& rTolerance
    )
    {
        PointType center{};
        for (PointView r_vertex : rPolygon.Vertices()) { center += r_vertex; }
        center /= static_cast<double>(rPolygon.NumberOfVertices());
        const auto bounds = rGridIndexer.BoundsXYZ();
        const Vector3i counts = rGridIndexer.ElementCounts();
        const IndexType cell = rGridIndexer.GetContainingCell(center);
        const auto cell_indices = rGridIndexer.GetMatrixIndicesFromVectorIndex(cell);
        for (IndexType axis = 0; axis < 3; ++axis) {
            const double delta = (bounds.upper[axis] - bounds.lower[axis]) / static_cast<double>(counts[axis]);
            const double grid_position = (center[axis] - bounds.lower[axis]) / delta;
            const IndexType plane = static_cast<IndexType>(std::llround(grid_position));
            const double coordinate = bounds.lower[axis] + delta * static_cast<double>(plane);
            bool is_on_plane = plane <= counts[axis];
            for (PointView r_vertex : rPolygon.Vertices()) {
                is_on_plane = is_on_plane && rTolerance.CoordinatesAreSame(r_vertex[axis], coordinate);
            }
            if (!is_on_plane) { continue; }

            GridIndexer::Direction direction{};
            // Lower cell boundaries belong to the preceding cell's forward face; upper boundaries belong to this cell.
            if (plane == 0 || plane == cell_indices[axis]) {
                direction = static_cast<GridIndexer::Direction>(2 * axis + 1);
            } else {
                direction = static_cast<GridIndexer::Direction>(2 * axis);
            }
            return rGridIndexer.GetFace(cell, direction);
        }
        return std::nullopt;
    }

}  // namespace

template<CellProduct TProduct>
MeshPartitioner<TProduct>::MeshPartitioner(const TriangleMeshView& rMesh, const GridIndexer& rGridIndexer)
    : mCellProducts(rGridIndexer.NumberOfElements()), mGridFaceSurfaces(rGridIndexer.NumberOfFaces())
{
    const GeometryTolerance& r_tolerance = rGridIndexer.GetGeometryTolerance();
    const auto number_of_triangles = static_cast<std::int64_t>(rMesh.NumOfTriangles());
    std::vector<std::vector<PolygonRecord>> cell_records(static_cast<IndexType>(number_of_triangles));
    std::vector<std::vector<PolygonRecord>> face_records(static_cast<IndexType>(number_of_triangles));

    // Clip each source triangle to the grid, split it at internal planes, and classify its pieces by target section.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::int64_t raw_triangle_id = 0; raw_triangle_id < number_of_triangles; ++raw_triangle_id) {
        const IndexType triangle_id = static_cast<IndexType>(raw_triangle_id);
        auto clipped_polygon =
            ClipToGrid(MakePolygon(rMesh.Triangle<WithNormals>(triangle_id)), rGridIndexer.BoundsXYZ(), r_tolerance);
        if (!clipped_polygon) { continue; }

        std::vector<ConvexPolygon> polygons{ std::move(*clipped_polygon) };
        for (IndexType axis = 0; axis < 3 && !polygons.empty(); ++axis) {
            SplitAtInternalPlanes(polygons, rGridIndexer, axis, r_tolerance);
        }
        IndexType piece = 0;
        for (auto& r_polygon : polygons) {
            if (const auto target = GetFaceTarget(r_polygon, rGridIndexer, r_tolerance)) {
                face_records[triangle_id].push_back({ target->value, triangle_id, piece++, std::move(r_polygon) });
                continue;
            }

            PointType center{};
            for (PointView r_vertex : r_polygon.Vertices()) { center += r_vertex; }
            center /= static_cast<double>(r_polygon.NumberOfVertices());
            const IndexType cell_index = rGridIndexer.GetContainingCell(center);
            cell_records[triangle_id].push_back({ cell_index, triangle_id, piece++, std::move(r_polygon) });
        }
    }

    // Merge per-triangle results in source order, then sort by target and piece order for deterministic output.
    std::vector<PolygonRecord> cells;
    std::vector<PolygonRecord> faces;
    for (IndexType triangle_id = 0; triangle_id < static_cast<IndexType>(number_of_triangles); ++triangle_id) {
        std::ranges::move(cell_records[triangle_id], std::back_inserter(cells));
        std::ranges::move(face_records[triangle_id], std::back_inserter(faces));
    }
    const auto Less = [](const PolygonRecord& rLeft, const PolygonRecord& rRight) {
        return std::tie(rLeft.target, rLeft.source_triangle, rLeft.piece)
               < std::tie(rRight.target, rRight.source_triangle, rRight.piece);
    };
    std::sort(cells.begin(), cells.end(), Less);
    std::sort(faces.begin(), faces.end(), Less);

    if constexpr (TProduct == CellProduct::Domain) {
        std::vector<std::optional<PendingCellSurfaceSection>> pending_sections(rGridIndexer.NumberOfElements());
        for (const auto& r_record : cells) {
            auto& r_pending = pending_sections[r_record.target];
            if (!r_pending) {
                r_pending.emplace();
                mCellIndices.push_back(r_record.target);
            }
            detail::AppendCellSurfacePolygon(
                r_record.polygon,
                r_record.source_triangle,
                rGridIndexer.GetBoundingBoxXYZFromIndex(r_record.target),
                r_tolerance,
                r_pending->surface,
                r_pending->contours
            );
        }
        for (const IndexType cell_index : mCellIndices) {
            auto& r_pending = pending_sections[cell_index];
            QuESo_ASSERT(r_pending.has_value(), "Closure input was not assembled.");
            SortFaceSegments(r_pending->contours);
            mCellProducts[cell_index].emplace(
                std::move(r_pending->surface),
                std::move(r_pending->contours),
                rGridIndexer.GetBoundingBoxXYZFromIndex(cell_index),
                r_tolerance
            );
        }
    } else {
        for (const auto& r_record : cells) {
            auto& r_surface = mCellProducts[r_record.target];
            if (!r_surface) {
                r_surface.emplace();
                mCellIndices.push_back(r_record.target);
            }
            detail::Triangulate(r_record.polygon, *r_surface);
        }
    }
    for (const auto& r_record : faces) {
        auto& r_section = mGridFaceSurfaces[r_record.target];
        if (!r_section) {
            r_section.emplace();
            mFaceIds.push_back({ r_record.target });
        }
        detail::Triangulate(r_record.polygon, *r_section);
    }
}

template<CellProduct TProduct>
typename MeshPartitioner<TProduct>::CellProductType
    MeshPartitioner<TProduct>::TakeCellProduct(IndexType CellIndex) noexcept(NOTDEBUG)
{
    QuESo_ASSERT(CellIndex < mCellProducts.size(), "Cell index is out-of-bounds.");
    QuESo_ASSERT(mCellProducts[CellIndex].has_value(), "Cell product is missing or already consumed.");
    CellProductType result = std::move(*mCellProducts[CellIndex]);
    mCellProducts[CellIndex].reset();
    return result;
}

template<CellProduct TProduct>
TriangleMesh MeshPartitioner<TProduct>::TakeGridFaceSurface(GridFaceId Face) noexcept(NOTDEBUG)
{
    QuESo_ASSERT(Face.value < mGridFaceSurfaces.size(), "Grid-face ID is out-of-bounds.");
    QuESo_ASSERT(mGridFaceSurfaces[Face.value].has_value(), "Grid-face surface is missing or already consumed.");
    TriangleMesh result = std::move(*mGridFaceSurfaces[Face.value]);
    mGridFaceSurfaces[Face.value].reset();
    return result;
}

template class MeshPartitioner<CellProduct::Surface>;
template class MeshPartitioner<CellProduct::Domain>;

}  // namespace queso::embedding
