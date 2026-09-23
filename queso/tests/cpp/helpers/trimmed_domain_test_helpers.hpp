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
#include <algorithm>
#include <cmath>

//// Project includes
#include "queso/embedding/mesh_operator.h"
#include "queso/embedding/trimmed_domain.h"

namespace queso::Testing::TrimmedDomainTestHelpers {

[[nodiscard]] inline GeometryTolerance MakeTolerance(const BoundingBoxType& rBounds)
{
    const PointType delta = rBounds.upper - rBounds.lower;
    const double length_scale = std::max({ delta[0], delta[1], delta[2] });
    const double coordinate_scale = std::max(
        { 1.0,
          std::abs(rBounds.lower[0]),
          std::abs(rBounds.lower[1]),
          std::abs(rBounds.lower[2]),
          std::abs(rBounds.upper[0]),
          std::abs(rBounds.upper[1]),
          std::abs(rBounds.upper[2]) }
    );
    return GeometryTolerance::FromScale({ .length_scale = length_scale, .coordinate_scale = coordinate_scale });
}

[[nodiscard]] inline Unique<TrimmedDomain> MakeTrimmedDomain(
    const embedding::MeshOperator& rMeshOperator,
    const BoundingBoxType& rBounds,
    IndexType MinNumberOfTriangles
)
{
    const GeometryTolerance tolerance = MakeTolerance(rBounds);
    BoundingBoxType candidate_bounds = rBounds;
    for (IndexType axis = 0; axis < 3; ++axis) {
        candidate_bounds.lower[axis] -= tolerance.SnapDistance();
        candidate_bounds.upper[axis] += tolerance.SnapDistance();
    }
    const auto candidates = rMeshOperator.Query().GetAabbCandidates(candidate_bounds.lower, candidate_bounds.upper);
    auto section = rMeshOperator.ClipCellSurfaceSection(candidates, rBounds);
    if (section.SurfaceView().NumOfTriangles() == 0) { return nullptr; }
    return MakeUnique<TrimmedDomain>(
        std::move(section),
        [&rMeshOperator](PointView rPoint) { return rMeshOperator.IsInside(rPoint); },
        MinNumberOfTriangles
    );
}

}  // namespace queso::Testing::TrimmedDomainTestHelpers
