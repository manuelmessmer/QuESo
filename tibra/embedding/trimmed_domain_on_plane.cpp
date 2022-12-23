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
    const auto& edges_on_planes = rTriangleMesh.GetEdgesOnPlane();
    const IndexType plane_index = DIRINDEX3*2UL + static_cast<IndexType>(mUpperBoundary);

    const auto& r_vertices = rTriangleMesh.GetVertices();
    const auto& edges_on_plane = edges_on_planes[plane_index];

    for( const auto& edge : edges_on_plane ){
        const IndexType vertex_index_1 = std::get<0>(edge);
        const IndexType vertex_index_2 = std::get<1>(edge);
        const IndexType normal_index = std::get<2>(edge);

        const auto& P1 = r_vertices[vertex_index_1];
        const auto& P2 = r_vertices[vertex_index_2];
        const auto& normal = rTriangleMesh.Normal(normal_index);
        InsertEdge(P1, P2, normal);
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