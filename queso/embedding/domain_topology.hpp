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
#include <cstdint>
#include <functional>
#include <optional>
#include <vector>

//// Project includes
#include "queso/containers/grid_indexer.hpp"
#include "queso/embedding/cell_surface_section.h"

namespace queso::embedding::detail {

/// @brief Traversal state of one canonical physical grid face.
enum class GridFaceState : std::uint8_t {
    traversable,  ///< Flood fill may traverse between the adjacent cells.
    blocked  ///< Flood fill cannot traverse and accumulates the directed vote for its current side.
};

/// @brief Stores one traversal state per physical grid face.
/// @details Opposite cell-local directions resolving to the same physical face share one canonical state. Interior
///          faces are initially traversable; exterior faces are always blocked. The referenced GridIndexer must
///          outlive this object.
class GridFaceStates
{
public:
    /// @brief Creates face states for one grid.
    /// @param rGridIndexer Grid defining the canonical physical faces and exterior boundaries.
    explicit GridFaceStates(const GridIndexer& rGridIndexer)
        : mGridIndexer(rGridIndexer), mStates(rGridIndexer.NumberOfFaces(), GridFaceState::traversable)
    {
        const IndexType number_of_cells = rGridIndexer.NumberOfElements();
        for (IndexType cell_index = 0; cell_index < number_of_cells; ++cell_index) {
            const Vector3i indices = rGridIndexer.GetMatrixIndicesFromVectorIndex(cell_index);
            const Vector3i number_of_elements = rGridIndexer.ElementCounts();
            if (indices[0] == 0) { Get(cell_index, GridIndexer::Direction::x_backward) = GridFaceState::blocked; }
            if (indices[0] + 1 == number_of_elements[0]) {
                Get(cell_index, GridIndexer::Direction::x_forward) = GridFaceState::blocked;
            }
            if (indices[1] == 0) { Get(cell_index, GridIndexer::Direction::y_backward) = GridFaceState::blocked; }
            if (indices[1] + 1 == number_of_elements[1]) {
                Get(cell_index, GridIndexer::Direction::y_forward) = GridFaceState::blocked;
            }
            if (indices[2] == 0) { Get(cell_index, GridIndexer::Direction::z_backward) = GridFaceState::blocked; }
            if (indices[2] + 1 == number_of_elements[2]) {
                Get(cell_index, GridIndexer::Direction::z_forward) = GridFaceState::blocked;
            }
        }
    }

    /// @brief Returns the state of a canonical physical face.
    /// @param Face Physical grid-face identity.
    /// @return Current traversal state.
    [[nodiscard]] GridFaceState Get(GridFaceId Face) const noexcept(NOTDEBUG)
    { return mStates[Face.value]; }

    /// @brief Returns mutable access to the state of a canonical physical face.
    /// @param Face Physical grid-face identity.
    /// @return Mutable traversal state.
    GridFaceState& Get(GridFaceId Face) noexcept(NOTDEBUG)
    { return mStates[Face.value]; }

    /// @brief Returns the state of a face addressed from one adjacent cell.
    /// @param CellIndex Adjacent cell index.
    /// @param FaceDirection Cell-local direction identifying the face.
    /// @return Current traversal state of the resolved physical face.
    [[nodiscard]] GridFaceState Get(IndexType CellIndex, GridIndexer::Direction FaceDirection) const noexcept(NOTDEBUG)
    { return Get(mGridIndexer.get().GetFace(CellIndex, FaceDirection)); }

    /// @brief Returns mutable access to a face addressed from one adjacent cell.
    /// @param CellIndex Adjacent cell index.
    /// @param FaceDirection Cell-local direction identifying the face.
    /// @return Mutable traversal state of the resolved physical face.
    GridFaceState& Get(IndexType CellIndex, GridIndexer::Direction FaceDirection) noexcept(NOTDEBUG)
    { return Get(mGridIndexer.get().GetFace(CellIndex, FaceDirection)); }

private:
    std::reference_wrapper<const GridIndexer> mGridIndexer;
    std::vector<GridFaceState> mStates;
};

/// @brief Stores one signed vote for each directed side of every physical grid face.
/// @details A zero vote provides no classification evidence. When flood fill encounters a blocked face, a positive
///          vote favors inside classification and a negative vote favors outside classification for that component.
///          Each physical face stores separate values for its negative and positive adjacent sides. Votes on exterior
///          sides are unused because exterior-connected components are classified as outside. The referenced
///          GridIndexer must outlive this object.
class GridFaceVotes
{
public:
    /// @brief Creates zero-initialized directed votes for one grid.
    /// @param rGridIndexer Grid defining the canonical physical faces and their adjacent sides.
    explicit GridFaceVotes(const GridIndexer& rGridIndexer)
        : mGridIndexer(rGridIndexer), mVotes(rGridIndexer.NumberOfFaces(), { 0, 0 })
    {}

    /// @brief Returns the vote for one side of a canonical physical face.
    /// @param Face Physical grid-face identity.
    /// @param Side Adjacent side whose vote is requested.
    /// @return Signed classification vote.
    [[nodiscard]] int Get(GridFaceId Face, GridFaceSide Side) const noexcept(NOTDEBUG)
    { return mVotes[Face.value][static_cast<IndexType>(Side)]; }

    /// @brief Returns mutable access to the vote for one side of a canonical physical face.
    /// @param Face Physical grid-face identity.
    /// @param Side Adjacent side whose vote is requested.
    /// @return Mutable signed classification vote.
    int& Get(GridFaceId Face, GridFaceSide Side) noexcept(NOTDEBUG)
    { return mVotes[Face.value][static_cast<IndexType>(Side)]; }

    /// @brief Returns the vote of the side adjacent to one cell-local face.
    /// @param CellIndex Adjacent cell index.
    /// @param FaceDirection Cell-local direction identifying the face.
    /// @return Signed classification vote for the resolved adjacent side.
    [[nodiscard]] int Get(IndexType CellIndex, GridIndexer::Direction FaceDirection) const noexcept(NOTDEBUG)
    {
        const auto& r_grid_indexer = mGridIndexer.get();
        const GridFaceId face = r_grid_indexer.GetFace(CellIndex, FaceDirection);
        return Get(face, r_grid_indexer.GetAdjacentSide(face, CellIndex));
    }

    /// @brief Returns mutable access to the vote of the side adjacent to one cell-local face.
    /// @param CellIndex Adjacent cell index.
    /// @param FaceDirection Cell-local direction identifying the face.
    /// @return Mutable signed classification vote for the resolved adjacent side.
    int& Get(IndexType CellIndex, GridIndexer::Direction FaceDirection) noexcept(NOTDEBUG)
    {
        const auto& r_grid_indexer = mGridIndexer.get();
        const GridFaceId face = r_grid_indexer.GetFace(CellIndex, FaceDirection);
        return Get(face, r_grid_indexer.GetAdjacentSide(face, CellIndex));
    }

private:
    std::reference_wrapper<const GridIndexer> mGridIndexer;
    std::vector<std::array<int, 2>> mVotes;
};

/// @brief Prepared trimmed-cell state and optional exact interior section for one domain cell.
struct DomainCellTopology
{
    /// @brief Whether this cell contains mesh boundary geometry and is excluded from flood fill.
    bool is_trimmed = false;
    /// @brief Explicit source-surface section used by the replacement closure path.
    std::optional<CellSurfaceSection> surface_section;
    /// @brief Whether the frozen input has been moved into a TrimmedDomain.
    bool input_consumed = false;
};

/// @brief Storage for cell, canonical face, and directed vote topology for one domain mesh and grid.
/// @details Construction initializes all cells as untrimmed, all interior faces as traversable, all exterior faces as
///          blocked, and all directed votes as zero. DomainMeshEmbedder then populates mesh-derived trimmed cells and
///          sections, blocks intersected and trimmed-cell-adjacent faces, and assigns oriented face votes before
///          production flood fill. Prepared trimmed cells have a boundary section; synthetic flood-fill tests may
///          omit it when no trimmed-domain construction is required.
class DomainTopology
{
public:
    explicit DomainTopology(const GridIndexer& rGridIndexer)
        : cell_topology(rGridIndexer.NumberOfElements()), face_topology(rGridIndexer), face_votes(rGridIndexer),
          mGridIndexer(rGridIndexer)
    {}

    /// @brief Returns the authoritative grid defining all stored topology.
    /// @return Grid indexer that must outlive this topology.
    [[nodiscard]] const GridIndexer& GetGridIndexer() const noexcept
    { return mGridIndexer.get(); }

    std::vector<DomainCellTopology> cell_topology;  ///< Per-cell trimmed state and optional boundary section.
    GridFaceStates face_topology;  ///< Canonical physical-face traversal barriers.
    GridFaceVotes face_votes;  ///< Directed classification votes on canonical physical faces.

private:
    std::reference_wrapper<const GridIndexer> mGridIndexer;  ///< Authoritative non-owning grid reference.
};

}  // namespace queso::embedding::detail
