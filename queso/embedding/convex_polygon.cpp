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
#include "queso/embedding/convex_polygon.h"
#include "queso/utilities/math_utilities.hpp"

namespace queso::embedding::detail {
namespace {

    [[nodiscard]] PointType
        Intersection(PointView rFirst, PointView rSecond, IndexType Axis, double Coordinate) noexcept
    {
        const double parameter = (Coordinate - rFirst[Axis]) / (rSecond[Axis] - rFirst[Axis]);
        PointType result{ rFirst[0] + parameter * (rSecond[0] - rFirst[0]),
                          rFirst[1] + parameter * (rSecond[1] - rFirst[1]),
                          rFirst[2] + parameter * (rSecond[2] - rFirst[2]) };
        result[Axis] = Coordinate;
        return result;
    }

    void AddDistinctVertex(ConvexPolygon& rPolygon, PointView rPoint, double Tolerance)
    {
        if (rPolygon.NumberOfVertices() > 0) {
            const PointView last = rPolygon.Vertex(rPolygon.NumberOfVertices() - 1);
            if (Math::SquaredNorm(rPoint - last) <= Tolerance * Tolerance) { return; }
        }
        rPolygon.AddVertex(rPoint);
    }

    [[nodiscard]] std::optional<ConvexPolygon> Finalize(ConvexPolygon&& rPolygon, const GeometryTolerance& rTolerance)
    {
        if (rPolygon.NumberOfVertices() > 1
            && Math::SquaredNorm(rPolygon.Vertex(0) - rPolygon.Vertex(rPolygon.NumberOfVertices() - 1))
                   <= rTolerance.SnapDistance() * rTolerance.SnapDistance()) {
            rPolygon.RemoveLastVertex();
        }
        if (rPolygon.NumberOfVertices() < 3 || rPolygon.Area() <= rTolerance.ZeroArea()) { return std::nullopt; }
        return std::move(rPolygon);
    }

}  // namespace

double ConvexPolygon::Area() const noexcept
{
    if (mNumberOfVertices < 3) { return 0.0; }
    PointType area_vector{};
    const PointView origin = mVertices[0];
    for (IndexType i = 1; i + 1 < mNumberOfVertices; ++i) {
        area_vector += Math::Cross(mVertices[i] - origin, mVertices[i + 1] - origin);
    }
    return 0.5 * Math::Norm(area_vector);
}

PolygonSplit
    SplitByPlane(const ConvexPolygon& rPolygon, IndexType Axis, double Coordinate, const GeometryTolerance& rTolerance)
{
    QuESo_ASSERT(Axis < 3, "Plane axis is out-of-bounds.");
    const double snap_distance = rTolerance.SnapDistance();

    ConvexPolygon negative(rPolygon.Normal());
    ConvexPolygon positive(rPolygon.Normal());
    if (rPolygon.NumberOfVertices() == 0) { return {}; }

    bool has_negative_vertex = false;
    bool has_positive_vertex = false;

    const PointView last = rPolygon.Vertex(rPolygon.NumberOfVertices() - 1);
    PointType first{ last[0], last[1], last[2] };
    double first_distance = first[Axis] - Coordinate;
    if (std::abs(first_distance) <= snap_distance) {
        first[Axis] = Coordinate;
        first_distance = 0.0;
    }

    for (PointView r_vertex : rPolygon.Vertices()) {
        PointType second{ r_vertex[0], r_vertex[1], r_vertex[2] };
        double second_distance = second[Axis] - Coordinate;
        has_negative_vertex = has_negative_vertex || second_distance < -snap_distance;
        has_positive_vertex = has_positive_vertex || second_distance > snap_distance;
        if (std::abs(second_distance) <= snap_distance) {
            second[Axis] = Coordinate;
            second_distance = 0.0;
        }

        if ((first_distance < 0.0 && second_distance > 0.0) || (first_distance > 0.0 && second_distance < 0.0)) {
            const PointType intersection = Intersection(first, second, Axis, Coordinate);
            AddDistinctVertex(negative, intersection, snap_distance);
            AddDistinctVertex(positive, intersection, snap_distance);
        }
        if (second_distance <= 0.0) { AddDistinctVertex(negative, second, snap_distance); }
        if (second_distance >= 0.0) { AddDistinctVertex(positive, second, snap_distance); }

        first = second;
        first_distance = second_distance;
    }

    const bool is_on_plane = !has_negative_vertex && !has_positive_vertex;
    return { has_negative_vertex || is_on_plane ? Finalize(std::move(negative), rTolerance) : std::nullopt,
             has_positive_vertex || is_on_plane ? Finalize(std::move(positive), rTolerance) : std::nullopt };
}

void Triangulate(const ConvexPolygon& rPolygon, TriangleMesh& rMesh)
{
    if (rPolygon.NumberOfVertices() < 3) { return; }
    const auto MakePoint = [](PointView rPoint) { return PointType{ rPoint[0], rPoint[1], rPoint[2] }; };
    if (rPolygon.NumberOfVertices() == 3) {
        const IndexType first = rMesh.AddVertex(MakePoint(rPolygon.Vertex(0)));
        const IndexType second = rMesh.AddVertex(MakePoint(rPolygon.Vertex(1)));
        const IndexType third = rMesh.AddVertex(MakePoint(rPolygon.Vertex(2)));
        rMesh.AddTriangle({ first, second, third }, MakePoint(rPolygon.Normal()));
        return;
    }

    PointType center{};
    for (PointView r_vertex : rPolygon.Vertices()) { center += r_vertex; }
    center /= static_cast<double>(rPolygon.NumberOfVertices());
    for (IndexType i = 0; i < rPolygon.NumberOfVertices(); ++i) {
        const IndexType next = (i + 1) % rPolygon.NumberOfVertices();
        const IndexType first = rMesh.AddVertex(MakePoint(rPolygon.Vertex(i)));
        const IndexType second = rMesh.AddVertex(MakePoint(rPolygon.Vertex(next)));
        const IndexType centroid = rMesh.AddVertex(center);
        rMesh.AddTriangle({ first, second, centroid }, MakePoint(rPolygon.Normal()));
    }
}

void AppendCellSurfacePolygon(
    const ConvexPolygon& rPolygon,
    IndexType SourceTriangle,
    const BoundingBoxType& rBounds,
    const GeometryTolerance& rTolerance,
    TriangleMesh& rSurface,
    CellFaceContours& rContours
)
{
    const Vector3d source_normal{ rPolygon.Normal()[0], rPolygon.Normal()[1], rPolygon.Normal()[2] };
    for (IndexType edge = 0; edge < rPolygon.NumberOfVertices(); ++edge) {
        const IndexType next = (edge + 1) % rPolygon.NumberOfVertices();
        const PointView r_first = rPolygon.Vertex(edge);
        const PointView r_second = rPolygon.Vertex(next);
        for (IndexType face_index = 0; face_index < 6; ++face_index) {
            const CellFace face = CellFace::FromIndex(face_index);
            const double coordinate = face.is_upper ? rBounds.upper[face.axis] : rBounds.lower[face.axis];
            if (!rTolerance.CoordinatesAreSame(r_first[face.axis], coordinate)
                || !rTolerance.CoordinatesAreSame(r_second[face.axis], coordinate)) {
                continue;
            }

            const auto SnapToFaceRectangle = [&](PointView rPoint) noexcept(NOTDEBUG) {
                PointType result{ rPoint[0], rPoint[1], rPoint[2] };
                result[face.axis] = coordinate;
                for (IndexType axis = 0; axis < 3; ++axis) {
                    if (axis == face.axis) { continue; }
                    QuESo_ASSERT(
                        result[axis] >= rBounds.lower[axis] - rTolerance.SnapDistance()
                            && result[axis] <= rBounds.upper[axis] + rTolerance.SnapDistance(),
                        "Cell-face segment endpoint is outside the face rectangle."
                    );
                    if (rTolerance.CoordinatesAreSame(result[axis], rBounds.lower[axis])) {
                        result[axis] = rBounds.lower[axis];
                    } else if (rTolerance.CoordinatesAreSame(result[axis], rBounds.upper[axis])) {
                        result[axis] = rBounds.upper[axis];
                    }
                    result[axis] = std::clamp(result[axis], rBounds.lower[axis], rBounds.upper[axis]);
                }
                return result;
            };

            PointType first = SnapToFaceRectangle(r_first);
            PointType second = SnapToFaceRectangle(r_second);
            PointType face_normal{};
            face_normal[face.axis] = face.is_upper ? 1.0 : -1.0;
            const Vector3d desired_direction = Math::Cross(face_normal, source_normal);
            const double alignment = Math::Dot(second - first, desired_direction);
            if (rTolerance.IsZeroLength(alignment) || alignment < 0.0) { std::swap(first, second); }
            rContours.segments[face_index].push_back({ first, second, source_normal, SourceTriangle });
        }
    }
    Triangulate(rPolygon, rSurface);
}

}  // namespace queso::embedding::detail
