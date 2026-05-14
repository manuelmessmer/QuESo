/*
  ____        ______  _____
 / __ \      |  ____|/ ____|
| |  | |_   _| |__  | (___   ___
| |  | | | | |  __|  \___ \ / _ \
| |__| | |_| | |____ ____) | (_) |
 \___\_\\__,_|______|_____/ \___/
        Quadrature for Embedded Solids

 License:    BSD 4-Clause License
             See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE

 Authors:    Manuel Messmer
*/

#ifndef CLIPPED_TRIANGLE_MESH_INCLUDE_HPP
#define CLIPPED_TRIANGLE_MESH_INCLUDE_HPP

//// STL includes
#include <array>
#include <concepts>
#include <ranges>
#include <utility>
#include <vector>

//// Project includes
#include "queso/containers/triangle_mesh_view.hpp"
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  ClippedTriangleMesh
 * @author Manuel Messmer
 * @brief  Triangle mesh container with additional edges on clipped planes.
 * @details ClippedTriangleMesh uses composition: it owns an internal TriangleMesh for triangle
 *          and vertex data, and adds a separate per-plane edge container to store clipped
 *          boundary edges.
 * @see    TriangleMesh.
 **/
class ClippedTriangleMesh
{
    ///@name Type Definitions
    ///@{
public:
    using ClippedTag = void;

    /// @brief Lightweight view of an edge on a clipping plane.
    struct EdgeProxy
    {
        PointView mP1;
        PointView mP2;
        Vector3dView mNormal;
    };

private:
    struct Edge
    {
        IndexType V1;
        IndexType V2;
        IndexType Normal;
    };

    using EdgeContainer = std::vector<Edge>;
    using EdgesPerPlane = std::array<EdgeContainer, 6>;

public:
    ///@}
    ///@name Life Cycle
    ///@{

    ClippedTriangleMesh()
    {}

	/// Takes ownership of the triangle mesh.
    ClippedTriangleMesh(TriangleMesh &&rTriangleMesh) : mTriangleMesh(std::move(rTriangleMesh))
    {}

	/// Takes ownership of the triangle mesh.
    ClippedTriangleMesh(Unique<TriangleMesh> &&pTriangleMesh) : mTriangleMesh(std::move(*pTriangleMesh))
    {}

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns a read-only mesh view.
    /// @return Triangle mesh view.
    TriangleMeshView MeshView() const
    { return TriangleMeshView(mTriangleMesh); }

    /// @brief Returns mutable access to the underlying mesh.
    /// @return Reference to owned triangle mesh.
    TriangleMesh &Mesh() noexcept
    { return mTriangleMesh; }

    /// @brief Returns const access to the underlying mesh.
    /// @return Const reference to owned triangle mesh.
    const TriangleMesh &Mesh() const noexcept
    { return mTriangleMesh; }

    /// @brief Reserves triangles container.
    /// @param Size Expected number of triangles.
    void Reserve(IndexType Size)
    { mTriangleMesh.Reserve(Size); }

    /// @brief Reserves per-plane edge storage.
    /// @param Size Expected number of edges per plane.
    void ReserveEdgesOnPlane(IndexType Size)
    {
        for (auto &edges : mEdgesOnPlanes) { edges.reserve(Size); }
    }

    /// @brief Appends one vertex.
    /// @param rVertex Vertex coordinates.
    /// @return Index of the new vertex.
    IndexType AddVertex(const Vector3d &rVertex)
    { return mTriangleMesh.AddVertex(rVertex); }

    /// @brief Appends one triangle with normal.
    /// @param rTriangle Triangle vertex indices.
    /// @param rNormal Triangle normal.
    void AddTriangle(const Vector3i &rTriangle, const Vector3d &rNormal)
    { mTriangleMesh.AddTriangle(rTriangle, rNormal); }

    /// @brief Returns number of triangles.
    /// @return Number of triangles.
    IndexType NumOfTriangles() const noexcept
    { return mTriangleMesh.NumOfTriangles(); }

    /// @brief Returns number of vertices.
    /// @return Number of vertices.
    IndexType NumOfVertices() const noexcept
    { return mTriangleMesh.NumOfVertices(); }

    /// @brief Returns a triangle proxy for the given triangle id.
    /// @tparam Mode Triangle mode.
    /// @param tid Triangle index.
    /// @return Triangle proxy.
    template<class Mode = WithNormals>
    TriangleProxy<Mode> Triangle(IndexType tid) const
    { return mTriangleMesh.Triangle<Mode>(tid); }

    /// @brief Returns lazy range over all triangles.
    /// @tparam Mode Triangle mode.
    /// @return Range of triangle proxies.
    template<class Mode = WithNormals>
    auto Triangles() const
    { return mTriangleMesh.Triangles<Mode>(); }

    /// @brief Returns lazy range over selected triangles.
    /// @tparam Mode Triangle mode.
    /// @tparam TIndexRange Range of triangle ids.
    /// @param indices Triangle ids.
    /// @return Range of triangle proxies.
    template<class Mode = WithNormals, std::ranges::viewable_range TIndexRange>
        requires std::convertible_to<std::ranges::range_value_t<TIndexRange>, IndexType>
    auto Triangles(TIndexRange &&indices) const
    { return mTriangleMesh.Triangles<Mode>(std::forward<TIndexRange>(indices)); }

    /// @brief Returns one triangle normal.
    /// @param tid Triangle index.
    /// @return Const normal reference.
    const Vector3d &Normal(IndexType tid) const noexcept(NOTDEBUG)
    { return mTriangleMesh.Normal(tid); }

    /// @brief Returns all vertices as a span.
    /// @return Read-only span of vertices.
    std::span<const Vector3d> Vertices() const noexcept(NOTDEBUG)
    { return mTriangleMesh.Vertices(); }

    /// @brief Add edge on a clipping plane.
    /// @param PlaneIndex Signed axis of plane.
    /// @param V1 Vertex index of first endpoint.
    /// @param V2 Vertex index of second endpoint.
    /// @param Normal Index of associated triangle normal.
    void AddEdgeOnPlane(SignedAxis PlaneIndex, IndexType V1, IndexType V2, IndexType Normal)
    {
        auto &edges = mEdgesOnPlanes[static_cast<IndexType>(PlaneIndex)];
        edges.push_back({ V1, V2, Normal });
    }

    /// @brief Returns edges on a given plane.
    /// @param PlaneIndex Signed axis of plane.
    /// @return Lazy range of edge proxies.
    auto Edges(SignedAxis PlaneIndex) const
    {
        const auto &edges = mEdgesOnPlanes[static_cast<IndexType>(PlaneIndex)];

        return edges | std::views::transform([this](const Edge &edge) {
            const auto &p1 = mTriangleMesh.Vertex(edge.V1);
            const auto &p2 = mTriangleMesh.Vertex(edge.V2);
            const auto &normal = mTriangleMesh.Normal(edge.Normal);
            return EdgeProxy{ PointView(p1), PointView(p2), Vector3dView(normal) };
        });
    }

    /// @brief Returns number of edges on a given plane.
    /// @param PlaneIndex Signed axis of plane.
    /// @return Number of edges.
    IndexType NumberOfEdges(SignedAxis PlaneIndex) const
    { return mEdgesOnPlanes[static_cast<IndexType>(PlaneIndex)].size(); }

private:
    ///@}
    ///@name Private Members
    ///@{

    TriangleMesh mTriangleMesh{};
    EdgesPerPlane mEdgesOnPlanes{};
    ///@}
};

///@}

}// namespace queso

#endif

