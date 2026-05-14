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
#include <concepts>
#include <functional>
#include <ranges>
#include <type_traits>
#include <utility>

//// Project includes
#include "queso/containers/triangle_proxies.hpp"
#include "queso/includes/define.hpp"

namespace queso::concepts {

/// Callable taking `const Mesh&` and returning a count convertible to `IndexType`.
template<typename Mesh, typename NumTriangles>
concept TriangleCountAccessor = requires(const Mesh &m, NumTriangles num_triangles) {
    { std::invoke(num_triangles, m) } -> std::convertible_to<IndexType>;
};

/// Callable taking (`const Mesh&`, `IndexType`) and returning `TriangleProxy<WithoutNormals>`.
template<typename Mesh, typename ToTriangle>
concept TriangleWithoutNormalsAccessor = requires(const Mesh &m, ToTriangle to_triangle, IndexType i) {
    { std::invoke(to_triangle, m, i) } -> std::same_as<TriangleProxy<WithoutNormals>>;
};

/// Empty, default-constructible callable policy type.
template<class TCallable>
concept StatelessCallablePolicy = std::default_initializable<TCallable> && std::is_empty_v<TCallable>;

/// Foreign mesh binding with valid accessors and stateless callable policies.
template<class TMesh, class TNumTriangles, class TToTriangle>
concept ForeignTriangleMeshBinding =
    TriangleCountAccessor<TMesh, TNumTriangles> && TriangleWithoutNormalsAccessor<TMesh, TToTriangle>
    && StatelessCallablePolicy<TNumTriangles> && StatelessCallablePolicy<TToTriangle>;

/// Viewable range whose values are convertible to `IndexType`.
template<class TRange>
concept IndexRange =
    std::ranges::viewable_range<TRange> && std::convertible_to<std::ranges::range_value_t<TRange>, IndexType>;

/// Callback invocable with both native mesh and mesh-view references.
template<class TCallBack, class TNativeMesh, class TMeshView>
concept TriangleMeshCallBack =
    std::invocable<TCallBack, const TNativeMesh &> && std::invocable<TCallBack, const TMeshView &>;

/// Viewable range whose value type is `TriangleProxy<Mode>`.
template<class Mode, class TRange>
concept TriangleProxyRange =
    std::ranges::viewable_range<TRange>
    && std::same_as<std::remove_cvref_t<std::ranges::range_value_t<TRange>>, TriangleProxy<Mode>>;

/// Triangle callback that returns `void` or a visit token.
template<class Mode, class TCallBack, class TVisitToken>
concept TriangleCallBack = std::invocable<TCallBack &, TriangleProxy<Mode>>
                           && (std::same_as<std::invoke_result_t<TCallBack &, TriangleProxy<Mode>>, void>
                               || std::same_as<std::invoke_result_t<TCallBack &, TriangleProxy<Mode>>, TVisitToken>);

/// `Triangles<Mode>()` return type satisfies `TriangleProxyRange`.
template<class TMesh, class Mode>
concept MeshTrianglesRange =
    TriangleProxyRange<Mode, decltype(std::declval<const TMesh &>().template Triangles<Mode>())>;

/// `Triangles<Mode>(indices)` return type satisfies `TriangleProxyRange`.
template<class TMesh, class Mode, class TIndexRange>
concept MeshIndexedTrianglesRange = TriangleProxyRange<Mode,
    decltype(std::declval<const TMesh &>().template Triangles<Mode>(std::declval<TIndexRange>()))>;

/// Callback invocable with `Triangles<Mode>()` range.
template<class TMesh, class Mode, class TCallBack>
concept MeshTrianglesCallBack =
    MeshTrianglesRange<TMesh, Mode>
    && std::invocable<TCallBack, decltype(std::declval<const TMesh &>().template Triangles<Mode>())>;

/// Callback invocable with `Triangles<Mode>(indices)` range.
template<class TMesh, class Mode, class TIndexRange, class TCallBack>
concept MeshIndexedTrianglesCallBack =
    MeshIndexedTrianglesRange<TMesh, Mode, TIndexRange>
    && std::invocable<TCallBack,
        decltype(std::declval<const TMesh &>().template Triangles<Mode>(std::declval<TIndexRange>()))>;

}// namespace queso::concepts

