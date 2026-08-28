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
#include <array>
#include <cmath>
#include <tuple>
#include <vector>

//// Project includes
#include "queso/embedding/cell_face_closure.h"
#include "queso/includes/numerical_guards.hpp"

namespace queso::embedding::detail {
namespace {

    constexpr double NormalComponentTolerance = 1e-14;

    using Point2d = std::array<double, 2>;

    enum class EdgeSide { Upper, Lower, Vertical };

    struct FaceFrame
    {
        CellFace face;
        bool switches_axes;
        IndexType u_axis;
        IndexType v_axis;
        double plane;
        Point2d lower;
        Point2d upper;

        [[nodiscard]] Point2d Project(PointView rPoint) const noexcept
        { return { rPoint[u_axis], rPoint[v_axis] }; }

        [[nodiscard]] PointType Lift(const Point2d& rPoint) const noexcept
        {
            PointType result{};
            result[face.axis] = plane;
            result[u_axis] = rPoint[0];
            result[v_axis] = rPoint[1];
            return result;
        }
    };

    struct PlaneEdge
    {
        Point2d first;
        Point2d second;
        Point2d source_normal;
        EdgeSide side;
    };

    struct ActiveEdge
    {
        const PlaneEdge* p_edge;
        double middle;
        double left;
        double right;
    };

    [[nodiscard]] FaceFrame MakeFrame(const BoundingBoxType& rBounds, CellFace Face, bool SwitchAxes) noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(Face.axis < 3, "Cell-face axis is out-of-bounds.");
        const IndexType u_axis = (Face.axis + (SwitchAxes ? 2 : 1)) % 3;
        const IndexType v_axis = (Face.axis + (SwitchAxes ? 1 : 2)) % 3;
        return { Face,
                 SwitchAxes,
                 u_axis,
                 v_axis,
                 Face.is_upper ? rBounds.upper[Face.axis] : rBounds.lower[Face.axis],
                 { rBounds.lower[u_axis], rBounds.lower[v_axis] },
                 { rBounds.upper[u_axis], rBounds.upper[v_axis] } };
    }

    [[nodiscard]] Point2d
        SnapToRectangle(Point2d Point, const FaceFrame& rFrame, const GeometryTolerance& rTolerance) noexcept
    {
        for (IndexType axis = 0; axis < 2; ++axis) {
            if (rTolerance.CoordinatesAreSame(Point[axis], rFrame.lower[axis])) {
                Point[axis] = rFrame.lower[axis];
            } else if (rTolerance.CoordinatesAreSame(Point[axis], rFrame.upper[axis])) {
                Point[axis] = rFrame.upper[axis];
            }
            Point[axis] = std::clamp(Point[axis], rFrame.lower[axis], rFrame.upper[axis]);
        }
        return Point;
    }

    [[nodiscard]] std::vector<PlaneEdge> ProjectAndClassifyEdges(
        std::span<const CellFaceSegment> rSegments,
        const FaceFrame& rFrame,
        const GeometryTolerance& rTolerance
    )
    {
        // In the local (u, v) face frame:
        //
        //       v
        //       ^   upper boundary: projected source normal points +v
        //       |      --------
        //       |      filled strip
        //       |      --------
        //       |   lower boundary: projected source normal points -v
        //       +-----------------> u
        //
        // Endpoints are ordered by u so every non-vertical edge can be interpolated inside a slab.
        std::vector<PlaneEdge> result;
        result.reserve(rSegments.size());
        for (const auto& r_segment : rSegments) {
            Point2d first = SnapToRectangle(rFrame.Project(r_segment.first), rFrame, rTolerance);
            Point2d second = SnapToRectangle(rFrame.Project(r_segment.second), rFrame, rTolerance);
            if (std::tie(second[0], second[1]) < std::tie(first[0], first[1])) { std::swap(first, second); }
            if (std::hypot(second[0] - first[0], second[1] - first[1]) <= rTolerance.ZeroLength()) { continue; }
            const Point2d normal{ r_segment.source_normal[rFrame.u_axis], r_segment.source_normal[rFrame.v_axis] };
            EdgeSide side = EdgeSide::Vertical;
            if (normal[1] > NormalComponentTolerance) {
                side = EdgeSide::Upper;
            } else if (normal[1] < -NormalComponentTolerance) {
                side = EdgeSide::Lower;
            }
            result.push_back({ first, second, normal, side });
        }
        std::sort(result.begin(), result.end(), [](const PlaneEdge& rLeft, const PlaneEdge& rRight) {
            return std::tie(rLeft.first, rLeft.second, rLeft.side) < std::tie(rRight.first, rRight.second, rRight.side);
        });
        result.erase(
            std::unique(
                result.begin(),
                result.end(),
                [&rTolerance](const PlaneEdge& rLeft, const PlaneEdge& rRight) {
                    return rLeft.side == rRight.side && rTolerance.CoordinatesAreSame(rLeft.first[0], rRight.first[0])
                           && rTolerance.CoordinatesAreSame(rLeft.first[1], rRight.first[1])
                           && rTolerance.CoordinatesAreSame(rLeft.second[0], rRight.second[0])
                           && rTolerance.CoordinatesAreSame(rLeft.second[1], rRight.second[1]);
                }
            ),
            result.end()
        );
        return result;
    }

    [[nodiscard]] double Interpolate(const PlaneEdge& rEdge, double U) noexcept
    {
        const double delta = rEdge.second[0] - rEdge.first[0];
        if (!numerical_guards::IsSafeDivisor(delta)) { return 0.5 * (rEdge.first[1] + rEdge.second[1]); }
        const double parameter = (U - rEdge.first[0]) / delta;
        return rEdge.first[1] + parameter * (rEdge.second[1] - rEdge.first[1]);
    }

    [[nodiscard]] bool ClassifyStripPoint(
        const FaceFrame& rFrame,
        const Point2d& rPoint,
        const GeometryTolerance& rTolerance,
        const std::function<bool(PointView)>& rIsInside
    )
    {
        // Move the query into the cell so it never lies on the cap plane or source contour.
        PointType point = rFrame.Lift(rPoint);
        point[rFrame.face.axis] += (rFrame.face.is_upper ? -4.0 : 4.0) * rTolerance.SnapDistance();
        return rIsInside(point);
    }

    void AddTriangle(
        TriangleMesh& rMesh,
        const FaceFrame& rFrame,
        const Point2d& rFirst,
        const Point2d& rSecond,
        const Point2d& rThird,
        const GeometryTolerance& rTolerance
    )
    {
        // Swapping u and v reverses the local frame. Reverse the triangle order when required so its
        // stored normal always points out of the cell.
        const double signed_area =
            0.5
            * ((rSecond[0] - rFirst[0]) * (rThird[1] - rFirst[1]) - (rSecond[1] - rFirst[1]) * (rThird[0] - rFirst[0]));
        if (std::abs(signed_area) <= rTolerance.ZeroArea()) { return; }
        PointType normal{};
        normal[rFrame.face.axis] = rFrame.face.is_upper ? 1.0 : -1.0;
        const PointType first = rFrame.Lift(rFirst);
        const PointType second = rFrame.Lift(rSecond);
        const PointType third = rFrame.Lift(rThird);
        const IndexType first_index = rMesh.AddVertex(first);
        if (rFrame.face.is_upper != rFrame.switches_axes) {
            const IndexType second_index = rMesh.AddVertex(second);
            const IndexType third_index = rMesh.AddVertex(third);
            rMesh.AddTriangle({ first_index, second_index, third_index }, normal);
        } else {
            const IndexType second_index = rMesh.AddVertex(third);
            const IndexType third_index = rMesh.AddVertex(second);
            rMesh.AddTriangle({ first_index, second_index, third_index }, normal);
        }
    }

    void AddStrip(
        TriangleMesh& rMesh,
        const FaceFrame& rFrame,
        double LeftU,
        double RightU,
        double LowerLeft,
        double LowerRight,
        double UpperLeft,
        double UpperRight,
        const GeometryTolerance& rTolerance
    )
    {
        // A valid material strip has its upper boundary at or above its lower boundary at both endpoints.
        if (UpperLeft < LowerLeft - rTolerance.SnapDistance() || UpperRight < LowerRight - rTolerance.SnapDistance()) {
            return;
        }
        const Point2d lower_left{ LeftU, LowerLeft };
        const Point2d lower_right{ RightU, LowerRight };
        const Point2d upper_right{ RightU, UpperRight };
        const Point2d upper_left{ LeftU, UpperLeft };
        AddTriangle(rMesh, rFrame, lower_left, lower_right, upper_right, rTolerance);
        AddTriangle(rMesh, rFrame, lower_left, upper_right, upper_left, rTolerance);
    }

}  // namespace

TriangleMesh BuildCellFaceClosure(
    std::span<const CellFaceSegment> rSegments,
    const BoundingBoxType& rCellBounds,
    CellFace Face,
    const GeometryTolerance& rTolerance,
    const std::function<bool(PointView)>& rIsInside,
    bool SwitchAxes
)
{
    const FaceFrame frame = MakeFrame(rCellBounds, Face, SwitchAxes);
    const std::vector<PlaneEdge> edges = ProjectAndClassifyEdges(rSegments, frame, rTolerance);

    // Every endpoint u-coordinate splits the face into slabs:
    //
    //   u0       u1       u2
    //   |--------|--------|
    //
    // Within one open slab, non-vertical edges cannot begin, end, or cross, so their v-order is stable.
    std::vector<double> events{ frame.lower[0], frame.upper[0] };
    events.reserve(events.size() + 2 * edges.size());
    for (const PlaneEdge& r_edge : edges) {
        events.push_back(r_edge.first[0]);
        events.push_back(r_edge.second[0]);
    }
    std::sort(events.begin(), events.end());
    events.erase(
        std::unique(
            events.begin(),
            events.end(),
            [&rTolerance](double Left, double Right) { return rTolerance.CoordinatesAreSame(Left, Right); }
        ),
        events.end()
    );

    TriangleMesh result;
    for (IndexType event = 0; event + 1 < events.size(); ++event) {
        const double left_u = events[event];
        const double right_u = events[event + 1];
        if (right_u - left_u <= rTolerance.ZeroLength()) { continue; }
        const double middle_u = 0.5 * (left_u + right_u);
        std::vector<ActiveEdge> active;
        for (const PlaneEdge& r_edge : edges) {
            if (r_edge.side == EdgeSide::Vertical || r_edge.second[0] - r_edge.first[0] <= rTolerance.ZeroLength()
                || middle_u <= r_edge.first[0] || middle_u >= r_edge.second[0]) {
                continue;
            }
            active.push_back(
                { &r_edge, Interpolate(r_edge, middle_u), Interpolate(r_edge, left_u), Interpolate(r_edge, right_u) }
            );
        }
        std::sort(active.begin(), active.end(), [](const ActiveEdge& rLeft, const ActiveEdge& rRight) {
            return rLeft.middle < rRight.middle;
        });

        if (active.empty()) {
            const Point2d sample{ middle_u, 0.5 * (frame.lower[1] + frame.upper[1]) };
            if (ClassifyStripPoint(frame, sample, rTolerance, rIsInside)) {
                AddStrip(
                    result,
                    frame,
                    left_u,
                    right_u,
                    frame.lower[1],
                    frame.lower[1],
                    frame.upper[1],
                    frame.upper[1],
                    rTolerance
                );
            }
            continue;
        }

        // Sweep bottom to top. A lower edge enters material and an upper edge leaves it. Coincident
        // transitions are one event; ambiguous orientation is resolved with a point just above the event.
        bool is_inside = active.front().p_edge->side == EdgeSide::Upper;
        double lower_left = frame.lower[1];
        double lower_right = frame.lower[1];
        for (IndexType i = 0; i < active.size();) {
            IndexType end = i + 1;
            while (end < active.size() && rTolerance.CoordinatesAreSame(active[i].middle, active[end].middle)) {
                ++end;
            }
            if (is_inside) {
                AddStrip(
                    result, frame, left_u, right_u, lower_left, lower_right, active[i].left, active[i].right, rTolerance
                );
            }

            bool has_upper = false;
            bool has_lower = false;
            for (IndexType edge = i; edge < end; ++edge) {
                has_upper = has_upper || active[edge].p_edge->side == EdgeSide::Upper;
                has_lower = has_lower || active[edge].p_edge->side == EdgeSide::Lower;
            }
            if (has_upper != has_lower) {
                is_inside = has_lower;
            } else {
                const Point2d sample{ middle_u,
                                      std::min(frame.upper[1], active[i].middle + 4.0 * rTolerance.SnapDistance()) };
                is_inside = ClassifyStripPoint(frame, sample, rTolerance, rIsInside);
            }
            lower_left = active[i].left;
            lower_right = active[i].right;
            i = end;
        }
        if (is_inside) {
            AddStrip(
                result, frame, left_u, right_u, lower_left, lower_right, frame.upper[1], frame.upper[1], rTolerance
            );
        }
    }
    return result;
}

}  // namespace queso::embedding::detail
