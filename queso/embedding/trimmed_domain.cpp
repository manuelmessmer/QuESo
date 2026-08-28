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

//// Project includes
#include "queso/embedding/trimmed_domain.h"
#include "queso/embedding/ray.h"
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso {

bool TrimmedDomain::IsInsideTrimmedDomain(const PointType& rPoint) const
{
    bool success = true;
    const bool val = IsInsideTrimmedDomain(rPoint, success);
    return success && val;
}

bool TrimmedDomain::IsInsideTrimmedDomain(const PointType& rPoint, bool& rSuccess) const
{

    const IndexType num_triangles = mClippedMesh.NumOfTriangles();
    if (num_triangles == 0) { return true; }
    rSuccess = true;
    bool success_local = false;
    bool is_inside = false;
    IndexType current_id = 0;
    while (!success_local) {
        // Return false if all triangles are tested, but non valid (all are parallel or on_boundary)
        if (current_id >= num_triangles) {
            rSuccess = false;
            return false;
        }

        // Get direction
        const auto triangle = mClippedMesh.Triangle<WithoutNormals>(current_id);
        const auto center_triangle = TriangleUtilities::Center(triangle);
        Vector3d direction = center_triangle - rPoint;

        if (Math::SquaredNorm(direction) <= mGeometryTolerance.ZeroLength() * mGeometryTolerance.ZeroLength()) {
            ++current_id;
            continue;
        }

        // Construct ray
        embedding::Ray ray(rPoint, direction);

        // Make sure the target triangle is not parallel or degenerate.
        if (!ray.IsParallel(triangle)) {
            const auto classification = mMeshQuery.Classify(ray);
            is_inside = classification.is_inside;
            success_local = classification.is_conclusive;
        }
        current_id++;
    }
    return is_inside;
}

IntersectionStateType TrimmedDomain::GetIntersectionState(
    const PointType& rLowerBound,
    const PointType& rUpperBound,
    AabbIntersectionPolicy Policy
) const
{
    const auto candidate_ids = mMeshQuery.GetAabbCandidates(rLowerBound, rUpperBound);
    if (mMeshQuery.IntersectsAabb(candidate_ids, rLowerBound, rUpperBound, Policy)) {
        return IntersectionState::trimmed;
    }

    // Test if center is inside or outside.
    const PointType center = (0.5 * (rLowerBound + rUpperBound));
    const auto status = IsInsideTrimmedDomain(center) ? IntersectionState::inside : IntersectionState::outside;

    // If triangle is not intersected, center location will determine if inside or outside.
    return status;
}

}  // End namespace queso
