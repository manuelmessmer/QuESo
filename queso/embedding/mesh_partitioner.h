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
#include <optional>
#include <span>
#include <type_traits>
#include <vector>

//// Project includes
#include "queso/containers/grid_indexer.hpp"
#include "queso/containers/triangle_mesh_view.hpp"
#include "queso/embedding/cell_surface_section.h"

namespace queso::embedding {

/// @brief Selects the cell product built by `MeshPartitioner`.
enum class CellProduct { Surface, Domain };

/// @brief Partitions one triangle mesh into unique background-grid cell and physical-face sections.
/// @details Each source triangle is first clipped to the closed outer grid AABB, then split at every intersected
///          internal grid plane. Pieces lying entirely on a grid face are assigned to that face's canonical section;
///          all other pieces are assigned to the cell containing their centroid. Results are merged in a deterministic
///          order before triangulation. Domain consumers receive explicit contour provenance; surface consumers and
///          physical faces receive plain triangle meshes. Section index lists remain valid after sections are taken,
///          and each section may be taken exactly once.
template<CellProduct TProduct>
class MeshPartitioner
{
public:
    using CellProductType = std::conditional_t<TProduct == CellProduct::Surface, TriangleMesh, CellSurfaceSection>;

    /// @brief Partitions one source mesh against one regular grid.
    /// @param rMesh Source mesh view, needed only during construction.
    /// @param rGridIndexer Authoritative grid topology and geometry.
    MeshPartitioner(const TriangleMeshView& rMesh, const GridIndexer& rGridIndexer);

    /// @brief Returns sorted indices of cells with non-empty sections.
    /// @return Stable non-empty cell-index span.
    [[nodiscard]] std::span<const IndexType> GetCellIndices() const noexcept
    { return mCellIndices; }

    /// @brief Returns sorted IDs of physical faces with non-empty sections.
    /// @return Stable non-empty face-ID span.
    [[nodiscard]] std::span<const GridFaceId> GetFaceIds() const noexcept
    { return mFaceIds; }

    /// @brief Transfers one non-empty cell product.
    /// @param CellIndex Cell listed by `GetCellIndices()`.
    /// @return Owned plain surface or closure input selected by TProduct.
    [[nodiscard]] CellProductType TakeCellProduct(IndexType CellIndex) noexcept(NOTDEBUG);

    /// @brief Transfers one non-empty canonical grid-face surface.
    /// @param Face Face listed by `GetFaceIds()`.
    /// @return Owned source surface lying on Face.
    [[nodiscard]] TriangleMesh TakeGridFaceSurface(GridFaceId Face) noexcept(NOTDEBUG);

private:
    std::vector<std::optional<CellProductType>> mCellProducts;  ///< Non-empty products indexed by cell.
    std::vector<std::optional<TriangleMesh>> mGridFaceSurfaces;  ///< Canonical source surfaces indexed by grid face.
    std::vector<IndexType> mCellIndices;  ///< Sorted non-empty cell targets for TProduct.
    std::vector<GridFaceId> mFaceIds;  ///< Sorted non-empty canonical grid-face targets.
};

extern template class MeshPartitioner<CellProduct::Surface>;
extern template class MeshPartitioner<CellProduct::Domain>;

}  // namespace queso::embedding
