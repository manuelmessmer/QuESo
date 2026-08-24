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
//// Project includes
#include "queso/embedding/convex_polygon.h"

namespace queso::embedding::detail {

/// @brief Clips one candidate triangle against exact bounds using a resolved geometry tolerance.
/// @param rTriangle Triangle with its source orientation.
/// @param rBounds Exact clipping bounds.
/// @param rTolerance Authoritative geometry tolerance.
/// @param RejectAligned Whether triangles lying completely on a clipping plane are excluded.
/// @return Clipped polygon, or `std::nullopt` for rejected or disjoint geometry.
[[nodiscard]] std::optional<ConvexPolygon> ClipTriangle(
    TriangleProxy<WithNormals> rTriangle,
    const BoundingBoxType& rBounds,
    const GeometryTolerance& rTolerance,
    bool RejectAligned
);

}  // namespace queso::embedding::detail
