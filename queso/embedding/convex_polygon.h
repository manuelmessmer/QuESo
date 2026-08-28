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

#pragma once

//// STL includes
#include <array>
#include <optional>
#include <span>

//// Project includes
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/embedding/cell_surface_section.h"
#include "queso/includes/define.hpp"

namespace queso::embedding::detail {

/// @brief Planar convex polygon with compact inline vertex storage.
class ConvexPolygon
{
public:
    /// @brief Maximum number of vertices produced by clipping an input triangle against six cell planes.
    /// @details Each clipping plane can add at most one vertex, so three input vertices require at most nine slots.
    static constexpr IndexType max_vertices = 9;

    /// @brief Constructs an empty polygon with the source surface normal.
    /// @param rNormal Source surface normal preserved by splitting and triangulation.
    explicit ConvexPolygon(PointView rNormal) noexcept : mNormal{ rNormal[0], rNormal[1], rNormal[2] }
    {}

    /// @brief Returns the number of polygon vertices.
    /// @return Number of ordered vertices.
    [[nodiscard]] IndexType NumberOfVertices() const noexcept
    { return mNumberOfVertices; }

    /// @brief Returns one ordered polygon vertex.
    /// @param Index Vertex index.
    /// @return Read-only point view.
    [[nodiscard]] PointView Vertex(IndexType Index) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(Index < mNumberOfVertices, "Polygon vertex index is out-of-bounds.");
        return mVertices[Index];
    }

    /// @brief Returns all ordered polygon vertices.
    /// @return Read-only vertex span.
    [[nodiscard]] std::span<const PointType> Vertices() const noexcept
    { return { mVertices.data(), mNumberOfVertices }; }

    /// @brief Returns the source surface normal.
    /// @return Read-only normal view.
    [[nodiscard]] Vector3dView Normal() const noexcept
    { return mNormal; }

    /// @brief Computes polygon area.
    /// @return Non-negative polygon area.
    [[nodiscard]] double Area() const noexcept;

    /// @brief Appends one ordered vertex.
    /// @param rPoint Vertex coordinates.
    void AddVertex(PointView rPoint) noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(mNumberOfVertices < max_vertices, "ConvexPolygon vertex capacity exceeded.");
        mVertices[mNumberOfVertices++] = PointType{ rPoint[0], rPoint[1], rPoint[2] };
    }

    /// @brief Removes the final ordered vertex.
    void RemoveLastVertex() noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(mNumberOfVertices > 0, "Cannot remove a vertex from an empty polygon.");
        --mNumberOfVertices;
    }

private:
    std::array<PointType, max_vertices> mVertices{};
    IndexType mNumberOfVertices{};
    PointType mNormal{};
};

/// @brief Result of splitting a polygon by an axis-aligned plane.
struct PolygonSplit
{
    std::optional<ConvexPolygon> negative;
    std::optional<ConvexPolygon> positive;
};

/// @brief Splits a convex polygon by an axis-aligned plane using a resolved geometry tolerance.
/// @details Coordinates use the snap distance while degenerate polygons use the dimensionally resolved zero area.
/// @param rPolygon Polygon to split.
/// @param Axis Plane-normal axis in `[0, 3)`.
/// @param Coordinate Exact plane coordinate.
/// @param rTolerance Resolved geometry tolerance.
/// @return Optional negative-side and positive-side polygons.
[[nodiscard]] PolygonSplit
    SplitByPlane(const ConvexPolygon& rPolygon, IndexType Axis, double Coordinate, const GeometryTolerance& rTolerance);

/// @brief Triangulates a convex polygon into a plain triangle mesh.
/// @param rPolygon Polygon whose source orientation is preserved.
/// @param rMesh Target mesh.
void Triangulate(const ConvexPolygon& rPolygon, TriangleMesh& rMesh);

/// @brief Appends a polygon to an explicit cell surface and records its oriented perimeter traces before triangulation.
/// @param rPolygon Polygon whose source orientation is preserved.
/// @param SourceTriangle Source-triangle provenance ID.
/// @param rBounds Exact target-cell bounds.
/// @param rTolerance Authoritative cell geometry tolerance.
/// @param rSurface Target plain source-surface mesh.
/// @param rContours Target oriented face-contour storage.
void AppendCellSurfacePolygon(
    const ConvexPolygon& rPolygon,
    IndexType SourceTriangle,
    const BoundingBoxType& rBounds,
    const GeometryTolerance& rTolerance,
    TriangleMesh& rSurface,
    CellFaceContours& rContours
);

}  // namespace queso::embedding::detail
