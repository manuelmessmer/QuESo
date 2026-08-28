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
#include <bit>
#include <cmath>
#include <cstdint>

//// Project includes
#include "queso/embedding/clipper.h"
#include "queso/embedding/mesh_operator.h"
#include "queso/embedding/ray.h"

namespace queso::embedding {
namespace {

    [[nodiscard]] constexpr std::uint64_t Mix(std::uint64_t Value) noexcept
    {
        Value += 0x9e3779b97f4a7c15ULL;
        Value = (Value ^ (Value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        Value = (Value ^ (Value >> 27U)) * 0x94d049bb133111ebULL;
        return Value ^ (Value >> 31U);
    }

    [[nodiscard]] std::uint64_t CoordinateBits(double Value) noexcept
    {
        if (Value == 0.0) { Value = 0.0; }  // Canonicalize negative zero.
        return std::bit_cast<std::uint64_t>(Value);
    }

    [[nodiscard]] double UnitInterval(std::uint64_t Value) noexcept
    { return static_cast<double>(Value >> 11U) * 0x1.0p-53; }

    [[nodiscard]] Vector3d RayDirection(PointView rPoint, IndexType Attempt) noexcept
    {
        std::uint64_t state = Mix(CoordinateBits(rPoint[0]));
        state = Mix(state ^ CoordinateBits(rPoint[1]));
        state = Mix(state ^ CoordinateBits(rPoint[2]));
        state = Mix(state ^ static_cast<std::uint64_t>(Attempt));
        Vector3d direction{};
        for (IndexType axis = 0; axis < 3; ++axis) {
            state = Mix(state);
            direction[axis] = 0.5 + UnitInterval(state);
        }
        return direction;
    }

}  // namespace

bool MeshOperator::IsInside(PointView rPoint) const
{
    if (!std::isfinite(rPoint[0]) || !std::isfinite(rPoint[1]) || !std::isfinite(rPoint[2])) { return false; }
    if (!mQuery.IsWithinBoundingBox(rPoint)) { return false; }
    const PointType point{ rPoint[0], rPoint[1], rPoint[2] };

    IndexType iteration = 0;
    IndexType success_count = 0;
    int inside_count = 0;
    while (success_count < 5) {
        if (iteration >= 100) { return false; }
        const Vector3d direction = RayDirection(rPoint, iteration++);
        const auto [is_inside, success] = mQuery.Classify(Ray(point, direction));
        if (success) {
            ++success_count;
            inside_count += is_inside ? 1 : -1;
        }
    }
    return inside_count > 0;
}

CellSurfaceSection
    MeshOperator::ClipCellSurfaceSection(std::span<const IndexType> rCandidateIds, const BoundingBoxType& rBounds) const
{
    const auto intersected_ids = mQuery.GetIntersectedTriangleIds(rCandidateIds, rBounds.lower, rBounds.upper);
    TriangleMesh surface;
    surface.Reserve(2 * intersected_ids.size());
    CellFaceContours contours;
    for (const IndexType source_triangle : intersected_ids) {
        auto polygon = detail::ClipTriangle(
            mQuery.MeshView().Triangle<WithNormals>(source_triangle), rBounds, mGeometryTolerance, true
        );
        if (polygon) {
            detail::AppendCellSurfacePolygon(*polygon, source_triangle, rBounds, mGeometryTolerance, surface, contours);
        }
    }
    for (auto& r_segments : contours.segments) {
        std::sort(r_segments.begin(), r_segments.end(), [](const auto& rLeft, const auto& rRight) {
            return std::tie(rLeft.first, rLeft.second, rLeft.source_normal, rLeft.source_triangle)
                   < std::tie(rRight.first, rRight.second, rRight.source_normal, rRight.source_triangle);
        });
    }
    return CellSurfaceSection(std::move(surface), std::move(contours), rBounds, mGeometryTolerance);
}

}  // namespace queso::embedding
