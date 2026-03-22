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

#pragma once

//// STL includes
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

//// Project includes
#include "queso/containers/triangle_mesh_concepts.hpp"
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

// Forward declaration.
class TriangleMeshView;

/**
 * @class  TriangleMesh
 * @author Manuel Messmer
 * @brief  Lightweight owning container for triangle connectivity, vertex coordinates, and per-triangle normals.
 * @see    TriangleMeshView.
 **/
class TriangleMesh
{
    ///@name Operations
    ///@{

public:
    /// @brief Reserves memory for vertices and triangles.
    /// @param NumTriangles Expected number of triangles.
    void Reserve(SizeType NumTriangles)
    {
        mTriangles.reserve(NumTriangles);
        mNormals.reserve(NumTriangles);
        mVertices.reserve((NumTriangles + 1) / 2);
    }

    /// @brief Appends a new vertex.
    /// @param rVertex Vertex coordinates.
    /// @return Index of new vertex.
    IndexType AddVertex(const Vector3d &rVertex)
    {
        const IndexType vid = NumOfVertices();
        mVertices.push_back(rVertex);
        return vid;
    }

    /// @brief Appends a new triangle with a precomputed normal.
    /// @param rTriangle Vertex indices of the triangle.
    /// @param rNormal Triangle normal vector.
    void AddTriangle(const Vector3i &rTriangle, const Vector3d &rNormal)
    {
        mTriangles.push_back(rTriangle);
        mNormals.push_back(rNormal);
    }

    /// @brief Clears all mesh data.
    void Clear()
    {
        mVertices.clear();
        mTriangles.clear();
        mNormals.clear();
    }

    /// @brief Removes a triangle.
    /// @details Swaps the to be deleted triangle with the last one and then pops the last element.
    /// @param tid Triangle index to remove.
    void RemoveTriangle(IndexType tid)
    {
        QuESo_ASSERT(tid < mTriangles.size(), "Index is out-of-bounds.");
        mTriangles[tid] = mTriangles.back();
        mNormals[tid] = mNormals.back();
        mTriangles.pop_back();
        mNormals.pop_back();
    }

    /// @brief Returns the number of triangles.
    /// @return Number of stored triangles.
    IndexType NumOfTriangles() const noexcept
    { return mTriangles.size(); }

    /// @brief Returns the number of vertices.
    /// @return Number of stored vertices.
    IndexType NumOfVertices() const noexcept
    { return mVertices.size(); }

    /// @brief Returns a triangle proxy for the given triangle id.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @param tid Triangle index.
    /// @return Triangle proxy in the requested mode.
    /// @note May only throw in debug mode.
    template<class Mode = WithNormals>
    TriangleProxy<Mode> Triangle(IndexType tid) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(tid < mTriangles.size(), "Index is out-of-bounds.");
        const auto &tri = mTriangles[tid];

        if constexpr (std::is_same_v<Mode, WithNormals>) {
            return { mVertices[tri[0]], mVertices[tri[1]], mVertices[tri[2]], mNormals[tid] };
        } else {
            return { mVertices[tri[0]], mVertices[tri[1]], mVertices[tri[2]] };
        }
    }

    /// @brief Returns the stored normal for the given triangle id.
    /// @param tid Triangle index.
    /// @return Const reference to the triangle normal.
    /// @note May only throw in debug mode.
    const Vector3d &Normal(IndexType tid) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(tid < mNormals.size(), "Index is out-of-bounds.");
        return mNormals[tid];
    }

    /// @brief Returns the vertex coordinates for a given vertex id.
    /// @param vid Vertex index.
    /// @return Const reference to the vertex coordinates.
    const Vector3d &Vertex(IndexType vid) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(vid < mVertices.size(), "Index is out-of-bounds.");
        return mVertices[vid];
    }

    /// @brief Returns a lazy range over all triangles.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @return View of triangle proxies over all triangle indices.
    template<class Mode = WithNormals>
    auto Triangles() const
    { return Triangles<Mode>(std::views::iota(IndexType{ 0 }, NumOfTriangles())); }

    /// @brief Returns a lazy range over a subset of triangles.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @tparam TIndexRange Range type containing triangle ids. Range type must be convertible to `IndexType`.
    /// @param Indices Range of triangle indices to visit.
    /// @return View of triangle proxies for the provided indices.
    template<class Mode = WithNormals, concepts::IndexRange TIndexRange>
    auto Triangles(TIndexRange &&Indices) const
    {
        std::span<const Vector3i> triangles = mTriangles;
        std::span<const Vector3d> vertices = mVertices;
        auto triangle_indices = std::views::all(std::forward<TIndexRange>(Indices));

        if constexpr (std::is_same_v<Mode, WithNormals>) {
            std::span<const Vector3d> normals = mNormals;
            return triangle_indices | std::views::transform([triangles, vertices, normals](IndexType i) {
                const auto &tri = triangles[i];
                return TriangleProxy<WithNormals>{ vertices[tri[0]], vertices[tri[1]], vertices[tri[2]], normals[i] };
            });
        } else {
            return triangle_indices | std::views::transform([triangles, vertices](IndexType i) {
                const auto &tri = triangles[i];
                return TriangleProxy<WithoutNormals>{ vertices[tri[0]], vertices[tri[1]], vertices[tri[2]] };
            });
        }
    }

    /// @brief Returns all stored normals as a span.
    /// @return Read-only span of triangle normals.
    std::span<const Vector3d> Normals() const noexcept
    { return mNormals; }

    /// @brief Returns this mesh as a TriangleMeshView.
    /// @details The function is declared here but defined in `triangle_mesh_view.hpp` to avoid
    ///          an unnecessary include dependency and circular include between the two headers.
    /// @note To use this function, the call site must include `triangle_mesh_view.hpp`.
    /// @return Triangle mesh view over this mesh.
    TriangleMeshView View() const;

    /// @brief Returns all stored vertices as a span.
    /// @return Read-only span of vertices.
    std::span<const Vector3d> Vertices() const noexcept
    { return mVertices; }

    /// @brief Returns all triangle index triplets as a span.
    /// @return Read-only span of triangle index triplets.
    std::span<const Vector3i> TriangleIndices() const noexcept
    { return mTriangles; }

    /// @brief Basic consistency check of this triangle mesh instance.
    void Check() const
    {
        if (mTriangles.size() != mNormals.size()) {
            QuESo_ERROR << "Number of Triangles and Normals in mesh do not match.\n";
        }

        for (IndexType i = 0; i < mTriangles.size(); ++i) {
            for (IndexType j = 0; j < 3; ++j) {
                if (mTriangles[i][j] >= mVertices.size()) { QuESo_ERROR << "Triangle/Vertex mismatch.\n"; }
            }
        }
    }

private:
    ///@}
    ///@name Private member variables
    ///@{

    std::vector<Vector3d> mVertices;
    std::vector<Vector3i> mTriangles;
    std::vector<Vector3d> mNormals;
    ///@}
};// End of TriangleMesh.

///@}

}// End namespace queso.
