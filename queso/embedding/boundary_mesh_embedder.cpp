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
//// STL includes
//// Project includes
#include "queso/embedding/boundary_mesh_embedder.h"
#include "queso/embedding/mesh_partitioner.h"
#include "queso/utilities/mesh_utilities.h"

namespace queso::embedding {

BoundaryMeshEmbedder::BoundaryMeshEmbedder(const TriangleMeshView& rMesh, const BackgroundGridType& rBackgroundGrid)
    : mBackgroundGrid(std::cref(rBackgroundGrid)), mSectionCache(rBackgroundGrid.GetGridIndexer().NumberOfElements())
{
    // Active-element addresses and the authoritative cell-index map must be stable throughout preparation.
    QuESo_ASSERT(rBackgroundGrid.ElementsAreLocked(), "BoundaryMeshEmbedder requires a locked BackgroundGrid.");
    Prepare(rMesh);
}

IndexType BoundaryMeshEmbedder::GetFaceParent(GridFaceId Face) const
{
    const auto& r_background_grid = mBackgroundGrid.get();
    const auto& r_grid_indexer = r_background_grid.GetGridIndexer();
    const auto negative = r_grid_indexer.GetAdjacentCell(Face, GridFaceSide::negative);
    const auto positive = r_grid_indexer.GetAdjacentCell(Face, GridFaceSide::positive);
    const bool negative_is_active = negative && r_background_grid.GetElementView(*negative + 1).has_value();
    const bool positive_is_active = positive && r_background_grid.GetElementView(*positive + 1).has_value();
    if (negative_is_active != positive_is_active) { return negative_is_active ? *negative : *positive; }
    return r_grid_indexer.GetCanonicalAdjacentCell(Face);
}

void BoundaryMeshEmbedder::AssignSection(IndexType ParentCellIndex, TriangleMesh&& rSection)
{
    QuESo_ASSERT(ParentCellIndex < mSectionCache.size(), "Boundary parent cell index is out-of-bounds.");
    if (rSection.NumOfTriangles() == 0) { return; }
    auto& r_cached_section = mSectionCache[ParentCellIndex];
    if (r_cached_section.has_value()) {
        MeshUtilities::Append(*r_cached_section, rSection);
    } else {
        r_cached_section.emplace(std::move(rSection));
    }
}

void BoundaryMeshEmbedder::Prepare(const TriangleMeshView& rMesh)
{
    const auto& r_background_grid = mBackgroundGrid.get();
    const auto& r_grid_indexer = r_background_grid.GetGridIndexer();
    MeshPartitioner<CellProduct::Surface> partitioner(rMesh, r_grid_indexer);
    for (const IndexType cell_index : partitioner.GetCellIndices()) {
        AssignSection(cell_index, partitioner.TakeCellProduct(cell_index));
    }
    for (const GridFaceId face : partitioner.GetFaceIds()) {
        AssignSection(GetFaceParent(face), partitioner.TakeGridFaceSurface(face));
    }

    for (IndexType cell_index = 0; cell_index < mSectionCache.size(); ++cell_index) {
        if (mSectionCache[cell_index].has_value()) { mParentCellIndices.push_back(cell_index); }
    }
}

TriangleMesh BoundaryMeshEmbedder::TakeBoundarySection(IndexType CellIndex) noexcept(NOTDEBUG)
{
    QuESo_ASSERT(CellIndex < mSectionCache.size(), "Boundary cell index is out-of-bounds.");
    QuESo_ASSERT(mSectionCache[CellIndex].has_value(), "Boundary section is missing, unlisted, or already consumed.");
    TriangleMesh section = std::move(*mSectionCache[CellIndex]);
    mSectionCache[CellIndex].reset();
    return section;
}

}  // namespace queso::embedding
