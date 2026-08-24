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
#include "queso/embedding/clipper.h"

namespace queso::embedding::detail {

std::optional<ConvexPolygon> ClipTriangle(
    TriangleProxy<WithNormals> rTriangle,
    const BoundingBoxType& rBounds,
    const GeometryTolerance& rTolerance,
    bool RejectAligned
)
{
    ConvexPolygon polygon(rTriangle.Normal);
    polygon.AddVertex(rTriangle.P1);
    polygon.AddVertex(rTriangle.P2);
    polygon.AddVertex(rTriangle.P3);

    for (IndexType axis = 0; axis < 3; ++axis) {
        const auto IsAligned = [&](double Coordinate) {
            return rTolerance.CoordinatesAreSame(rTriangle.P1[axis], Coordinate)
                   && rTolerance.CoordinatesAreSame(rTriangle.P2[axis], Coordinate)
                   && rTolerance.CoordinatesAreSame(rTriangle.P3[axis], Coordinate);
        };
        if (RejectAligned && (IsAligned(rBounds.lower[axis]) || IsAligned(rBounds.upper[axis]))) {
            return std::nullopt;
        }

        auto lower = SplitByPlane(polygon, axis, rBounds.lower[axis], rTolerance);
        if (!lower.positive) { return std::nullopt; }
        auto upper = SplitByPlane(*lower.positive, axis, rBounds.upper[axis], rTolerance);
        if (!upper.negative) { return std::nullopt; }
        polygon = std::move(*upper.negative);
    }
    return polygon;
}

}  // namespace queso::embedding::detail
