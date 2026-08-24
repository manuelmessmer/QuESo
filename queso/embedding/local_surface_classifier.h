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

//// Project includes
#include "queso/containers/triangle_mesh_view.hpp"

namespace queso::embedding::detail {

/// @brief Tri-state classification against a nearby oriented local surface.
enum class LocalSurfaceClassification { Outside, Inside, Inconclusive };

/// @brief Classifies a point against a nearby clipped oriented surface.
/// @details Tests deterministic rays through section-triangle centers and applies a local majority rule without an
///          acceleration structure. Ambiguous rays and tied evidence remain inconclusive. Angular and barycentric
///          predicate thresholds remain dimensionless and separate from the supplied physical distance threshold.
/// @param rPoint Query point.
/// @param rBoundarySection Clipped oriented surface view.
/// @param ZeroLength Physical zero-distance threshold.
/// @return Tri-state local surface classification.
[[nodiscard]] LocalSurfaceClassification
    ClassifyOnBoundedSide(const PointType& rPoint, const TriangleMeshView& rBoundarySection, double ZeroLength);

}  // namespace queso::embedding::detail
