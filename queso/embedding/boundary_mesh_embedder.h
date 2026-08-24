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
//
#pragma once

//// STL includes
#include <functional>
#include <optional>
#include <span>
#include <vector>

//// Project includes
#include "queso/containers/background_grid.hpp"
#include "queso/containers/boundary_integration_point.hpp"

namespace queso::embedding {

/// @class BoundaryMeshEmbedder
/// @brief Eagerly prepares a possibly open boundary mesh against finalized background-grid elements.
/// @details The background grid must be locked before construction. Geometry overlapping or within snap tolerance of
///          the closed grid AABB is retained. Cell-interior sections keep their containing cell. A face section prefers
///          its sole active adjacent cell and otherwise uses its canonical adjacent cell. Internal faces canonically
///          select the positive-side cell; exterior faces canonically select their only adjacent cell. Inactive
///          canonical parents remain valid, and changing a parent never changes geometry orientation.
///
///          Sections assigned to the same parent are combined. Parent indices are exposed once in ascending order,
///          and each cached section can subsequently be moved out exactly once.
class BoundaryMeshEmbedder
{
public:
    ///@name Type Definitions
    ///@{

    using BackgroundGridType = BackgroundGrid<IntegrationPoint, BoundaryIntegrationPoint>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructs the open-mesh query and eagerly prepares all assignable boundary sections.
    /// @details The mesh may extend outside the background grid. The constructor performs all spatial queries and
    ///          clipping; the iteration and consumption operations do not perform further geometry work.
    /// @param rMesh Non-owning boundary mesh view that must remain valid throughout construction.
    /// @param rBackgroundGrid Locked grid containing the finalized active elements; it must outlive this embedder.
    BoundaryMeshEmbedder(const TriangleMeshView& rMesh, const BackgroundGridType& rBackgroundGrid);

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns parent cell indices assigned during preparation in ascending order.
    /// @details Every assigned parent occurs exactly once, independently of preparation order. The list remains stable
    ///          after individual sections are consumed.
    /// @return Read-only span of zero-based cell indices.
    [[nodiscard]] std::span<const IndexType> GetParentCellIndices() const noexcept
    { return mParentCellIndices; }

    /// @brief Moves the prepared section assigned to one listed parent cell.
    /// @details A section can be consumed exactly once. CellIndex must occur in GetParentCellIndices(), and its section
    ///          must not have been consumed. No clipping or mesh query is performed here.
    /// @param CellIndex Zero-based parent cell index listed by GetParentCellIndices().
    /// @return Owned prepared boundary section.
    [[nodiscard]] TriangleMesh TakeBoundarySection(IndexType CellIndex) noexcept(NOTDEBUG);

    ///@}
private:
    ///@name Private Operations
    ///@{

    /// @brief Partitions the source mesh, assigns every section once, and builds the sorted parent list.
    /// @param rMesh Boundary mesh view needed only during preparation.
    void Prepare(const TriangleMeshView& rMesh);

    /// @brief Adds one non-empty section to a parent's combined boundary mesh.
    /// @param ParentCellIndex Zero-based destination parent cell index.
    /// @param rSection Owned section to initialize or append to the parent's cache entry.
    void AssignSection(IndexType ParentCellIndex, TriangleMesh&& rSection);

    /// @brief Selects the sole active adjacent cell of a face, or its canonical adjacent cell otherwise.
    /// @param Face Physical grid face whose parent is requested.
    /// @return Zero-based parent cell index; exterior faces select their only adjacent cell.
    [[nodiscard]] IndexType GetFaceParent(GridFaceId Face) const;

    ///@}
    ///@name Private Members
    ///@{

    /// @brief Non-owning finalized grid reference used for topology and active-element ownership lookup.
    std::reference_wrapper<const BackgroundGridType> mBackgroundGrid;
    /// @brief Per-cell sections whose empty entries represent unassigned or already-consumed parents.
    std::vector<std::optional<TriangleMesh>> mSectionCache;
    /// @brief Stable assigned parent indices in deterministic ascending order.
    std::vector<IndexType> mParentCellIndices;

    ///@}
};

}  // namespace queso::embedding
