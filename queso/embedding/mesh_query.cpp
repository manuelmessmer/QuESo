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
#include "queso/embedding/aabb_triangle_intersector.h"
#include "queso/embedding/mesh_query.h"

namespace queso::embedding {

namespace {

#ifndef NDEBUG
    [[nodiscard]] bool BoundsAreValid(PointView rLowerBound, PointView rUpperBound) noexcept
    {
        for (IndexType axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(rLowerBound[axis]) || !std::isfinite(rUpperBound[axis])
                || rLowerBound[axis] > rUpperBound[axis]) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] bool BoundsSupportInset(PointView rLowerBound, PointView rUpperBound, double Inset) noexcept
    {
        if (!std::isfinite(Inset) || Inset < 0.0) { return false; }
        for (IndexType axis = 0; axis < 3; ++axis) {
            if (2.0 * Inset > rUpperBound[axis] - rLowerBound[axis]) { return false; }
        }
        return true;
    }
#endif

}  // namespace

bool MeshQuery::IsWithinBoundingBox(PointView rPoint) const noexcept
{
    if (!std::isfinite(rPoint[0]) || !std::isfinite(rPoint[1]) || !std::isfinite(rPoint[2])) { return false; }
    return mTriangleMesh.NumOfTriangles() > 0 && mTree.IsWithinBoundingBox(rPoint);
}

std::optional<BoundingBoxType> MeshQuery::BoundingBox() const noexcept
{
    if (mTriangleMesh.NumOfTriangles() == 0) { return std::nullopt; }
    return mTree.BoundingBox();
}

std::vector<IndexType> MeshQuery::GetAabbCandidates(PointView rLowerBound, PointView rUpperBound) const
{
#ifndef NDEBUG
    // TODO: Add a QuESo_DEBUG environment-controlled diagnostic mode that validates AABB contracts in release builds.
    QuESo_ASSERT(BoundsAreValid(rLowerBound, rUpperBound), "AABB bounds must be finite and non-inverted.");
#endif
    return mTree.Query(rLowerBound, rUpperBound);
}

bool MeshQuery::IntersectsAabb(
    std::span<const IndexType> rCandidateIds,
    PointView rLowerBound,
    PointView rUpperBound,
    AabbIntersectionPolicy Policy
) const
{
    const double inset = Policy == AabbIntersectionPolicy::SnapEroded ? mGeometryTolerance.SnapDistance() : 0.0;
#ifndef NDEBUG
    // TODO: Add a QuESo_DEBUG environment-controlled diagnostic mode that validates AABB contracts in release builds.
    QuESo_ASSERT(BoundsAreValid(rLowerBound, rUpperBound), "AABB bounds must be finite and non-inverted.");
    QuESo_ASSERT(
        Policy != AabbIntersectionPolicy::SnapEroded || BoundsSupportInset(rLowerBound, rUpperBound, inset),
        "Snap-eroded AABB must retain positive extent."
    );
#endif
    if constexpr (!NOTDEBUG) {
        for ([[maybe_unused]] const IndexType triangle_id : rCandidateIds) {
            QuESo_ASSERT(triangle_id < mTriangleMesh.NumOfTriangles(), "MeshQuery candidate ID is out-of-bounds.");
        }
    }
    const detail::AabbTriangleIntersector aabb(rLowerBound, rUpperBound);
    bool intersects = false;
    mTriangleMesh.VisitEachTriangle<WithoutNormals>(rCandidateIds, [&](const auto& rTriangle) {
        if (aabb.IntersectsTriangle(rTriangle, inset)) {
            intersects = true;
            return TriangleMeshView::VisitToken::stop_loop;
        }
        return TriangleMeshView::VisitToken::continue_loop;
    });
    return intersects;
}

std::vector<IndexType> MeshQuery::GetIntersectedTriangleIds(
    std::span<const IndexType> rCandidateIds,
    PointView rLowerBound,
    PointView rUpperBound
) const
{
#ifndef NDEBUG
    // TODO: Add a QuESo_DEBUG environment-controlled diagnostic mode that validates AABB contracts in release builds.
    QuESo_ASSERT(BoundsAreValid(rLowerBound, rUpperBound), "AABB bounds must be finite and non-inverted.");
#endif
    if constexpr (!NOTDEBUG) {
        for ([[maybe_unused]] const IndexType triangle_id : rCandidateIds) {
            QuESo_ASSERT(triangle_id < mTriangleMesh.NumOfTriangles(), "MeshQuery candidate ID is out-of-bounds.");
        }
    }
    const detail::AabbTriangleIntersector aabb(rLowerBound, rUpperBound);
    std::vector<IndexType> intersected_triangle_ids;
    intersected_triangle_ids.reserve(rCandidateIds.size());

    IndexType local_id = 0;
    mTriangleMesh.VisitEachTriangle<WithoutNormals>(rCandidateIds, [&](const auto& rTriangle) {
        const IndexType triangle_id = rCandidateIds[local_id++];
        if (aabb.IntersectsTriangle(rTriangle, 0.0)) { intersected_triangle_ids.push_back(triangle_id); }
    });
    return intersected_triangle_ids;
}

RayClassification MeshQuery::Classify(const Ray& rRay) const
{ return mMode == MeshQueryMode::Closed ? ClassifyClosed(rRay) : ClassifyOrientedSurface(rRay); }

RayClassification MeshQuery::ClassifyOrientedSurface(const Ray& rRay) const
{
    double min_distance = MAXD;
    bool is_inside = false;
    const auto potential_intersections = mTree.Query(rRay);
    std::optional<RayClassification> final_result;

    mTriangleMesh.VisitEachTriangle<WithoutNormals>(potential_intersections, [&](const auto& rTriangle) {
        const auto intersection = rRay.IntersectTriangle(rTriangle, mGeometryTolerance.ZeroLength());
        if (intersection.status != RayTriangleIntersectionStatus::Hit) {
            return TriangleMeshView::VisitToken::continue_loop;
        }
        const double sum_u_v = intersection.u + intersection.v;
        if (intersection.distance < mGeometryTolerance.ZeroLength()) {
            final_result = RayClassification{ false, true };
            return TriangleMeshView::VisitToken::stop_loop;
        }
        if (intersection.u < detail::ray_predicates::BarycentricBoundaryTolerance
            || intersection.v < detail::ray_predicates::BarycentricBoundaryTolerance
            || sum_u_v > 1.0 - detail::ray_predicates::BarycentricBoundaryTolerance) {
            final_result = RayClassification{ false, false };
            return TriangleMeshView::VisitToken::stop_loop;
        }
        if (intersection.distance < min_distance) {
            is_inside = intersection.is_back_facing;
            min_distance = intersection.distance;
        }
        return TriangleMeshView::VisitToken::continue_loop;
    });

    return final_result.value_or(RayClassification{ is_inside, true });
}

RayClassification MeshQuery::ClassifyClosed(const Ray& rRay) const
{
    const auto potential_intersections = mTree.Query(rRay);
    IndexType intersection_count = 0;
    std::optional<RayClassification> final_result;

    mTriangleMesh.VisitEachTriangle<WithoutNormals>(potential_intersections, [&](const auto& rTriangle) {
        const auto intersection = rRay.IntersectTriangle(rTriangle, mGeometryTolerance.ZeroLength());
        if (intersection.status == RayTriangleIntersectionStatus::Degenerate
            || intersection.status == RayTriangleIntersectionStatus::Miss) {
            return TriangleMeshView::VisitToken::continue_loop;
        }
        if (intersection.status == RayTriangleIntersectionStatus::Parallel) {
            final_result = RayClassification{ false, false };
            return TriangleMeshView::VisitToken::stop_loop;
        }
        ++intersection_count;
        const double sum_u_v = intersection.u + intersection.v;
        if (intersection.distance < mGeometryTolerance.ZeroLength()) {
            final_result = RayClassification{ false, true };
            return TriangleMeshView::VisitToken::stop_loop;
        }
        if (intersection.u < detail::ray_predicates::BarycentricBoundaryTolerance
            || intersection.v < detail::ray_predicates::BarycentricBoundaryTolerance
            || sum_u_v > 1.0 - detail::ray_predicates::BarycentricBoundaryTolerance) {
            final_result = RayClassification{ false, false };
            return TriangleMeshView::VisitToken::stop_loop;
        }
        return TriangleMeshView::VisitToken::continue_loop;
    });

    return final_result.value_or(RayClassification{ intersection_count % 2 == 1, true });
}

}  // namespace queso::embedding
