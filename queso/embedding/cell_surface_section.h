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
#include <span>
#include <utility>
#include <vector>

//// Project includes
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/containers/triangle_mesh_view.hpp"

namespace queso {

class TrimmedDomain;

namespace embedding {

    /// @brief Identifies one axis-aligned face of a background-grid cell.
    struct CellFace
    {
        IndexType axis;
        bool is_upper;

        /// @brief Returns the dense storage index in `[0, 6)`.
        /// @return Dense face index with lower/upper faces adjacent for each axis.
        [[nodiscard]] IndexType Index() const noexcept(NOTDEBUG)
        {
            QuESo_ASSERT(axis < 3, "Cell-face axis is out-of-bounds.");
            return 2 * axis + static_cast<IndexType>(is_upper);
        }

        /// @brief Constructs a cell face from its dense storage index.
        /// @param Index Dense face index in `[0, 6)`.
        /// @return Cell face identity.
        [[nodiscard]] static CellFace FromIndex(IndexType Index) noexcept(NOTDEBUG)
        {
            QuESo_ASSERT(Index < 6, "Cell-face index is out-of-bounds.");
            return { Index / 2, Index % 2 == 1 };
        }

        friend bool operator==(const CellFace&, const CellFace&) = default;
    };

    /// @brief Oriented source-surface trace on one cell face.
    /// @details In the outward-oriented face frame, the filled region lies to the left of `first -> second`.
    struct CellFaceSegment
    {
        PointType first;
        PointType second;
        Vector3d source_normal;
        IndexType source_triangle;
    };

    /// @brief Segment storage for all six faces of one cell.
    struct CellFaceContours
    {
        std::array<std::vector<CellFaceSegment>, 6> segments;
    };

    /// @brief Owns source-surface geometry and face-contour provenance for one cell.
    class CellSurfaceSection
    {
    public:
        /// @brief Constructs a cell surface section with checked mesh connectivity.
        /// @details Debug builds additionally validate bounds, clipped geometry, and face-segment invariants. Release
        ///          callers must provide finite bounds and geometry generated under the supplied tolerance.
        /// @param rSurface Clipped source-surface geometry.
        /// @param rContours Oriented traces grouped by cell face.
        /// @param rCellBounds Exact original cell bounds.
        /// @param rTolerance Geometry tolerance used to generate the section.
        CellSurfaceSection(
            TriangleMesh&& rSurface,
            CellFaceContours&& rContours,
            const BoundingBoxType& rCellBounds,
            const GeometryTolerance& rTolerance
        );

        CellSurfaceSection(const CellSurfaceSection&) = delete;
        CellSurfaceSection& operator=(const CellSurfaceSection&) = delete;
        CellSurfaceSection(CellSurfaceSection&&) noexcept = default;
        CellSurfaceSection& operator=(CellSurfaceSection&&) noexcept = default;

        /// @brief Returns a non-owning view of the clipped source surface.
        /// @return Surface view valid for this section's lifetime or until it is moved.
        [[nodiscard]] TriangleMeshView SurfaceView() const noexcept
        { return mSurface.View(); }

        /// @brief Returns oriented segments for one cell face.
        /// @param Face Requested face.
        /// @return Read-only segment span.
        [[nodiscard]] std::span<const CellFaceSegment> FaceSegments(CellFace Face) const noexcept(NOTDEBUG)
        { return mFaceContours.segments[Face.Index()]; }

        /// @brief Returns the exact original cell bounds.
        /// @return Read-only cell bounds.
        [[nodiscard]] const BoundingBoxType& CellBounds() const noexcept
        { return mCellBounds; }

        /// @brief Returns the geometry tolerance snapshot used to generate this section.
        /// @return Read-only geometry tolerance.
        [[nodiscard]] const GeometryTolerance& GetGeometryTolerance() const noexcept
        { return mGeometryTolerance; }

    private:
        /// TrimmedDomain may take mSurface.
        friend class ::queso::TrimmedDomain;

        [[nodiscard]] TriangleMesh TakeSurface() noexcept
        { return std::move(mSurface); }

        TriangleMesh mSurface;
        CellFaceContours mFaceContours;
        BoundingBoxType mCellBounds;
        GeometryTolerance mGeometryTolerance;
    };

}  // namespace embedding
}  // namespace queso
