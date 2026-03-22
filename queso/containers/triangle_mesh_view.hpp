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
#include <memory>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

//// Project includes
#include "queso/containers/triangle_mesh.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso {

/**
 * @class  TriangleMeshView
 * @author Manuel Messmer
 * @brief  Non-owning view that provides uniform triangle access for native and foreign mesh types.
 * @details TriangleMeshView adapts any triangle mesh implementation into a single non-owning
 *          interface that supports both native and non-native meshes. Native `TriangleMesh`
 *          instances use a fast path without virtual dispatch, while non-native meshes use a
 *          virtual-dispatch path via the internal function table.
 * @see    triangle_mesh.hpp and triangle_mesh_concepts.hpp.
 **/
class TriangleMeshView
{
    ///@name Life Cycle
    ///@{
public:
    enum class VisitToken { stop_loop, continue_loop };

private:
    /// Function table used for foreign mesh access.
    struct VTable
    {
        IndexType (*NumOfTriangles)(const void *);
        TriangleProxy<WithoutNormals> (*Triangle)(const void *, IndexType);
    };

    /// @brief Typed implementation of the function vtable.
    /// @tparam TMesh Foreign mesh type.
    /// @tparam TNumTriangles Callable returning the number of triangles.
    /// @tparam TToTriangle Callable returning TriangleProxy.
    template<class TMesh, class TNumTriangles, class TToTriangle>
    struct VTableImpl
    {
        inline static constexpr TNumTriangles mNumTriangles{};
        inline static constexpr TToTriangle mToTriangle{};

        static IndexType NumOfTriangles(const void *m)
        { return mNumTriangles(*static_cast<const TMesh *>(m)); }

        static TriangleProxy<WithoutNormals> TriangleWithoutNormal(const void *m, IndexType i)
        { return mToTriangle(*static_cast<const TMesh *>(m), i); }

        inline static constexpr VTable table{ &NumOfTriangles, &TriangleWithoutNormal };
    };

public:
    /// @brief Creates a view over a native `TriangleMesh`.
    /// @param mesh Native mesh instance.
    /// @see TriangleMesh.
    explicit TriangleMeshView(const TriangleMesh &mesh)
        : mMesh(static_cast<const void *>(&mesh)),
          mVTable(&GetVTable<TriangleMesh,
              decltype([](const TriangleMesh &m) { return m.NumOfTriangles(); }),
              decltype([](const TriangleMesh &m, IndexType i) { return m.Triangle<WithoutNormals>(i); })>()),
          mIsNativeMesh(true), mNormals(mesh.Normals())
    {}

    /// @brief Creates a view over a foreign mesh type.
    /// @tparam TMesh Foreign mesh type.
    /// @tparam TNumTriangles Callable returning the number of triangles.
    /// @tparam TToTriangle Callable converting an index to `TriangleProxy<WithoutNormals>`.
    /// @param rMesh Foreign mesh instance.
    /// @param NumTriangles Callable returning number of triangles.
    /// @param ToTriangle Callable returning one triangle proxy.
    /// @note Callables are treated as type-level stateless policies and are default-constructed
    ///       by the vtable implementation. Stateful or capturing callables are not supported.
    template<class TMesh, class TNumTriangles, class TToTriangle>
        requires concepts::ForeignTriangleMeshBinding<TMesh, TNumTriangles, TToTriangle>
    TriangleMeshView(const TMesh &rMesh, TNumTriangles NumTriangles, TToTriangle ToTriangle)
        : mMesh(static_cast<const void *>(&rMesh)), mVTable(&GetVTable<TMesh, TNumTriangles, TToTriangle>()),
          mNormals(std::make_shared<std::vector<Vector3d>>())
    {
        static_cast<void>(NumTriangles);
        static_cast<void>(ToTriangle);
        PopulateNormals();
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns number of triangles in the referenced mesh.
    /// @return Number of triangles.
    IndexType NumOfTriangles() const
    { return mVTable->NumOfTriangles(mMesh); }

    /// @brief Returns a triangle proxy for the given triangle id.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @param tid Triangle index.
    /// @return Triangle proxy in the requested mode.
    template<class Mode = WithNormals>
    TriangleProxy<Mode> Triangle(IndexType tid) const
    {
        QuESo_ASSERT(tid < NumOfTriangles(), "Index is out-of-bounds.");
        const auto tri = mVTable->Triangle(mMesh, tid);
        if constexpr (std::is_same_v<Mode, WithoutNormals>) {
            return tri;
        } else {
            const std::span<const Vector3d> normals = NormalsSpan();
            return TriangleProxy<WithNormals>{ tri.P1, tri.P2, tri.P3, normals[tid] };
        }
    }

    /// @brief Returns a lazy range over all triangles.
    /// @note Convenience API. Prefer `Visit`/`VisitEachTriangle` in performance-critical code.
    ///       `Visit` provides the native fast path, while this range-based access uses erased dispatch.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @return View of triangle proxies over all triangle indices.
    template<class Mode = WithNormals>
    auto Triangles() const
    { return Triangles<Mode>(std::views::iota(IndexType{ 0 }, NumOfTriangles())); }

    /// @brief Returns a lazy range over a subset of triangles.
    /// @note Convenience API. Prefer `Visit`/`VisitEachTriangle` in performance-critical code.
    ///       `Visit` provides the native fast path, while this range-based access uses erased dispatch.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @tparam TIndexRange Range type containing triangle ids.
    /// @param indices Range of triangle indices to visit.
    /// @return View of triangle proxies for the provided indices.
    template<class Mode = WithNormals, concepts::IndexRange TIndexRange>
    auto Triangles(TIndexRange &&indices) const
    {
        const auto mesh = mMesh;
        const auto fn = mVTable->Triangle;
        auto triangle_indices = std::views::all(std::forward<TIndexRange>(indices));

        if constexpr (std::is_same_v<Mode, WithoutNormals>) {
            return triangle_indices | std::views::transform([fn, mesh](IndexType i) { return fn(mesh, i); });
        } else {
            const std::span<const Vector3d> normals = NormalsSpan();
            return triangle_indices | std::views::transform([fn, mesh, normals](IndexType i) {
                const auto tri = fn(mesh, i);
                return TriangleProxy<WithNormals>{ tri.P1, tri.P2, tri.P3, normals[i] };
            });
        }
    }

    /// @brief Dispatches and provides all triangles as a range.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @tparam TCallBack Callable accepting the produced triangle range.
    /// @param rCallback Visitor callable.
    /// @return Return value of the visitor.
    template<class Mode = WithNormals, class TCallBack>
        requires concepts::MeshTrianglesCallBack<TriangleMeshView, Mode, TCallBack>
    decltype(auto) Visit(TCallBack &&rCallback) const
    { return Visit<Mode>(std::views::iota(IndexType{ 0 }, NumOfTriangles()), std::forward<TCallBack>(rCallback)); }

    /// @brief Dispatches and provides a specific set of triangles as a range.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @tparam TIndexRange Range type containing triangle ids.
    /// @tparam TCallBack Callable accepting the produced triangle range.
    /// @param rIndices Range of triangle indices to visit.
    /// @param rCallback Visitor callable.
    /// @return Return value of the visitor.
    template<class Mode = WithNormals, concepts::IndexRange TIndexRange, class TCallBack>
        requires concepts::MeshIndexedTrianglesCallBack<TriangleMeshView, Mode, TIndexRange, TCallBack>
    decltype(auto) Visit(TIndexRange &&rIndices, TCallBack &&rCallback) const
    {
        auto triangle_indices = std::views::all(std::forward<TIndexRange>(rIndices));
        return VisitDispatch([&](const auto &rMesh) -> decltype(auto) {
            return rCallback(rMesh.template Triangles<Mode>(triangle_indices));
        });
    }

    /// @brief Visits all triangles and invokes a callback per triangle.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @tparam TCallBack Callback type.
    /// @param rCallback Callback receiving one triangle proxy.
    template<class Mode = WithNormals, class TCallBack>
        requires concepts::TriangleCallBack<Mode, TCallBack, VisitToken>
    void VisitEachTriangle(TCallBack &&rCallback) const
    { VisitEachTriangle<Mode>(std::views::iota(IndexType{ 0 }, NumOfTriangles()), std::forward<TCallBack>(rCallback)); }

    /// @brief Visits a subset of triangles and invokes a callback per triangle.
    /// @tparam Mode Triangle access mode (`WithNormals` or `WithoutNormals`).
    /// @tparam TIndexRange Range type containing triangle ids.
    /// @tparam TCallBack Callback type.
    /// @param rIndices Range of triangle indices to visit.
    /// @param rCallback Callback receiving one triangle proxy. Callback may return a `VisitToken` to prematurely exit
    /// the loop.
    template<class Mode = WithNormals, concepts::IndexRange TIndexRange, class TCallBack>
        requires concepts::TriangleCallBack<Mode, TCallBack, VisitToken>
    void VisitEachTriangle(TIndexRange &&rIndices, TCallBack &&rCallback) const
    {
        auto triangle_indices = std::views::all(std::forward<TIndexRange>(rIndices));
        VisitDispatch([&](const auto &mesh) {
            for (const auto &tri : mesh.template Triangles<Mode>(triangle_indices)) {
                if constexpr (std::same_as<std::invoke_result_t<TCallBack &, TriangleProxy<Mode>>, VisitToken>) {
                    if (rCallback(tri) == VisitToken::stop_loop) { break; }
                } else {
                    rCallback(tri);
                }
            }
        });
    }

    ///@}

private:
    ///@name Private Operations
    ///@{
    /// @brief Dispatches to native mesh (fast path) or erased-view path.
    /// @details Use this dispatch for performance-sensitive traversal APIs.
    /// @tparam TCallBack Callable receiving either `const TriangleMesh&` or `const TriangleMeshView&`.
    /// @param rCallback Visitor callable.
    /// @return Return value of the visitor.
    template<class TCallBack>
        requires concepts::TriangleMeshCallBack<TCallBack, TriangleMesh, TriangleMeshView>
    decltype(auto) VisitDispatch(TCallBack &&rCallback) const
    {
        if (mIsNativeMesh) {
            const TriangleMesh &triangle_mesh = *static_cast<const TriangleMesh *>(mMesh);
            return rCallback(triangle_mesh);
        } else {
            return rCallback(*this);
        }
    }

    /// @brief Populates `mNormals` for foreign meshes.
    void PopulateNormals()
    {
        auto &r_normals = *std::get<std::shared_ptr<std::vector<Vector3d>>>(mNormals);
        r_normals.clear();
        r_normals.reserve(NumOfTriangles());

        VisitEachTriangle<WithoutNormals>(
            [&r_normals](const auto &r_triangle) { r_normals.push_back(TriangleUtilities::Normal(r_triangle)); });
    }

    /// @brief Returns normals as a span regardless of storage backend.
    std::span<const Vector3d> NormalsSpan() const
    {
        return std::visit(
            [](const auto &ns) -> std::span<const Vector3d> {
                if constexpr (std::is_same_v<std::decay_t<decltype(ns)>, std::shared_ptr<std::vector<Vector3d>>>) {
                    return std::span<const Vector3d>{ *ns };
                } else {
                    return ns;
                }
            },
            mNormals);
    }

    /// @brief Returns the static vtable for a foreign mesh binding.
    /// @tparam TMesh Foreign mesh type.
    /// @tparam TNumTriangles Callable returning the number of triangles.
    /// @tparam TToTriangle Callable returning one triangle proxy.
    /// @return Reference to static vtable instance.
    template<class TMesh, class TNumTriangles, class TToTriangle>
    static const VTable &GetVTable()
    { return VTableImpl<TMesh, TNumTriangles, TToTriangle>::table; }

    ///@}
    ///@name Private Members
    ///@{

    const void *mMesh{};
    const VTable *mVTable{};
    bool mIsNativeMesh = false;
    std::variant<std::shared_ptr<std::vector<Vector3d>>, std::span<const Vector3d>> mNormals;
    ///@}
};// End TriangleMeshView.

/// Defined here because this implementation only makes sense once `TriangleMeshView` exists
/// as a complete type in this header. Keeping it in `triangle_mesh.hpp` would require a
/// full include, introduce an unnecessary include dependency, and create a circular dependency.
inline TriangleMeshView TriangleMesh::View() const
{ return TriangleMeshView(*this); }

}// End namespace queso
