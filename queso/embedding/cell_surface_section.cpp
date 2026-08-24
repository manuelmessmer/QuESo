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
#include <cmath>

//// Project includes
#include "queso/embedding/cell_surface_section.h"

namespace queso::embedding {

#ifndef NDEBUG
namespace {

    [[nodiscard]] bool IsFinite(PointView rPoint) noexcept
    { return std::isfinite(rPoint[0]) && std::isfinite(rPoint[1]) && std::isfinite(rPoint[2]); }

}  // namespace
#endif

CellSurfaceSection::CellSurfaceSection(
    TriangleMesh&& rSurface,
    CellFaceContours&& rContours,
    const BoundingBoxType& rCellBounds,
    const GeometryTolerance& rTolerance
)
    : mSurface(std::move(rSurface)), mFaceContours(std::move(rContours)), mCellBounds(rCellBounds),
      mGeometryTolerance(rTolerance)
{
    mSurface.Check();

#ifndef NDEBUG
    // TODO: Add a QuESo_DEBUG environment-controlled diagnostic mode that runs these validations in release builds.
    const double snap_distance = mGeometryTolerance.SnapDistance();
    for (IndexType axis = 0; axis < 3; ++axis) {
        QuESo_ASSERT(
            std::isfinite(mCellBounds.lower[axis]) && std::isfinite(mCellBounds.upper[axis])
                && mCellBounds.lower[axis] < mCellBounds.upper[axis],
            "Cell-surface-section bounds must be finite and strictly increasing."
        );
        QuESo_ASSERT(
            mCellBounds.upper[axis] - mCellBounds.lower[axis] > 2.0 * snap_distance,
            "Cell-surface-section extent must exceed twice its snap distance."
        );
    }

    for (PointView r_vertex : mSurface.Vertices()) {
        QuESo_ASSERT(IsFinite(r_vertex), "Cell-surface-section vertices must be finite.");
    }
    for (Vector3dView r_normal : mSurface.Normals()) {
        QuESo_ASSERT(IsFinite(r_normal), "Cell-surface-section normals must be finite.");
    }

    for (IndexType face_index = 0; face_index < 6; ++face_index) {
        const CellFace face = CellFace::FromIndex(face_index);
        const double face_coordinate = face.is_upper ? mCellBounds.upper[face.axis] : mCellBounds.lower[face.axis];
        for (const CellFaceSegment& r_segment : mFaceContours.segments[face_index]) {
            QuESo_ASSERT(
                IsFinite(r_segment.first) && IsFinite(r_segment.second) && IsFinite(r_segment.source_normal),
                "Cell-face segment coordinates and normals must be finite."
            );
            QuESo_ASSERT(
                r_segment.first[face.axis] == face_coordinate && r_segment.second[face.axis] == face_coordinate,
                "Cell-face segment endpoints must lie exactly on their face plane."
            );
            for (IndexType axis = 0; axis < 3; ++axis) {
                if (axis == face.axis) { continue; }
                QuESo_ASSERT(
                    r_segment.first[axis] >= mCellBounds.lower[axis] && r_segment.first[axis] <= mCellBounds.upper[axis]
                        && r_segment.second[axis] >= mCellBounds.lower[axis]
                        && r_segment.second[axis] <= mCellBounds.upper[axis],
                    "Cell-face segment endpoints must lie within the exact face rectangle."
                );
            }
        }
    }
#endif
}

}  // namespace queso::embedding
