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
#include <functional>
#include <optional>

//// Project includes
#include "queso/embedding/domain_topology.hpp"
#include "queso/embedding/flood_fill.h"
#include "queso/embedding/mesh_operator.h"
#include "queso/embedding/trimmed_domain.h"

namespace queso::embedding {

/// @brief Prepares topology and cell surface sections for one closed domain mesh and background grid.
class DomainMeshEmbedder
{
public:
    /// @brief Constructs the mesh operator and eagerly prepares reusable topology and cell sections.
    /// @param rMesh Closed domain mesh view.
    /// @param rGridIndexer Authoritative background-grid indexer.
    DomainMeshEmbedder(const TriangleMeshView& rMesh, const GridIndexer& rGridIndexer);

    /// @brief Classifies all cells using prepared topology and flood fill only.
    /// @param Options Flood-fill execution options used on the first call.
    /// @return Cached states ordered by cell index.
    [[nodiscard]] ElementStates Classify(FloodFillOptions Options = {});

    /// @brief Consumes one prepared section and constructs its self-contained trimmed domain.
    /// @param CellIndex Zero-based trimmed-cell index.
    /// @param MinNumberOfBoundaryTriangles Minimum completed-boundary triangle count.
    /// @return Trimmed domain by value.
    [[nodiscard]] TrimmedDomain MakeTrimmedDomain(IndexType CellIndex, IndexType MinNumberOfBoundaryTriangles);

private:
    /// @brief Eagerly constructs cell and face topology required by flood fill and trimmed-domain construction.
    void PrepareTopology();

    /// @brief Partitions the source mesh against the grid and stores one consumable section for each trimmed cell.
    void PrepareCellTopology();

    /// @brief Records one blocked grid-face section and its orientation-derived votes for adjacent untrimmed cells.
    /// @param Face Grid face represented by rSection.
    /// @param rSection Oriented source surface section on Face.
    void PrepareFaceSection(GridFaceId Face, TriangleMesh&& rSection);

    /// @brief Marks every grid face adjacent to a trimmed cell as blocked.
    void PrepareFaceTopology();

    /// @brief Resolves remaining blocked-face traversal votes from adjacent trimmed-cell local surface sections.
    void PrepareFaceVotes();

    MeshOperator mMeshOperator;  ///< Owns the immutable source-mesh query and cell-scale clipping tolerance.
    /// Authoritative grid indexer. Must outlive this embedder.
    std::reference_wrapper<const GridIndexer> mGridIndexer;
    detail::DomainTopology mTopology;  ///< Prepared cell/face topology and consumable trimmed-cell sections.
    std::optional<ElementStates> mStates;  ///< Flood-fill result cached after the first Classify() call.
};

}  // namespace queso::embedding
