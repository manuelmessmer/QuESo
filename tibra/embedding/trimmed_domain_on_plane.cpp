// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes

//// Project includes
#include "embedding/trimmed_domain_on_plane.h"

/// Project includes

namespace tibra
{

bool TrimmedDomainOnPlane::InsertEdge(const Point2DType& rV1, const Point2DType& rV2, OrientationType rOrientation ){
    // Get unique vertex indices.
    auto indices = GetVertexID(rV1, rV2, rOrientation);
    // Only insert if indices are not the same.
    if (indices.first != indices.second) {
        InsertVertex(rV1, indices.first, rOrientation);
        InsertVertex(rV2, indices.second, rOrientation);
        auto& r_edges = GetEdges(rOrientation);
        r_edges.push_back(Egde2D(indices.first, indices.second));
        return true;
    }
    return false;
}

void TrimmedDomainOnPlane::InsertEdge(const Point3DType& rV1, const Point3DType& rV2, const Point3DType &rNormal)
{
    // Insantiate 2D points.
    Point2DType v1{};
    Point2DType v2{};

    // Make sure edge is oriented along DIRINDEX1 such that x2 > x1.
    if (rV2[DIRINDEX1] > rV1[DIRINDEX1]) {
        v1 = {rV1[DIRINDEX1], rV1[DIRINDEX2]};
        v2 = {rV2[DIRINDEX1], rV2[DIRINDEX2]};
    }
    else {
        v1 = {rV2[DIRINDEX1], rV2[DIRINDEX2]};
        v2 = {rV1[DIRINDEX1], rV1[DIRINDEX2]};
    }

    if ((rNormal[DIRINDEX2] > EPS3) ) { // Positive oriented
        InsertEdge(v1, v2, Orientation::Positive);
    }
    else if ((rNormal[DIRINDEX2] < -EPS3) ) { // Negative oriented
        InsertEdge(v1, v2, Orientation::Negative);
    }
    else { // Vertical oriented
        InsertEdge(v1, v2, Orientation::Vertical);
    }
}

void TrimmedDomainOnPlane::CollectEdgesOnPlane(const TriangleMesh &rTriangleMesh)
{
    for (IndexType triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id)
    {
        const auto &P1 = rTriangleMesh.P1(triangle_id);
        const auto &P2 = rTriangleMesh.P2(triangle_id);
        const auto &P3 = rTriangleMesh.P3(triangle_id);
        const auto &normal = rTriangleMesh.Normal(triangle_id);

        const bool is_p1 = IsPointOnPlane(P1);
        const bool is_p2 = IsPointOnPlane(P2);
        const bool is_p3 = IsPointOnPlane(P3);
        if (is_p1 && is_p2) {
            InsertEdge(P1, P2, normal);
        }
        else if (is_p2 && is_p3) {
            InsertEdge(P2, P3, normal);
        }
        else if (is_p3 && is_p1) {
            InsertEdge(P3, P1, normal);
        }
    }
}

bool TrimmedDomainOnPlane::IsPointOnPlane(const Point3DType &rPoint, double Tolerance) const
{
    if (mUpperBoundary) {
        if ((rPoint[DIRINDEX3] < mUpperBound[DIRINDEX3] + Tolerance)
            && (rPoint[DIRINDEX3] > mUpperBound[DIRINDEX3] - Tolerance)) {
            return true;
        }
        return false;
    }
    else {
        if ((rPoint[DIRINDEX3] < mLowerBound[DIRINDEX3] + Tolerance)
            && (rPoint[DIRINDEX3] > mLowerBound[DIRINDEX3] - Tolerance)) {
            return true;
        }
        return false;
    }
}

} // End namespace tibra