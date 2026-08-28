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

//// STL includes
#include <array>
#include <cstdint>
#include <optional>

//// Project includes
#include "queso/embedding/domain_mesh_embedder.h"
#include "queso/embedding/local_surface_classifier.h"
#include "queso/embedding/mesh_partitioner.h"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso::embedding {

DomainMeshEmbedder::DomainMeshEmbedder(const TriangleMeshView& rMesh, const GridIndexer& rGridIndexer)
    : mMeshOperator(rMesh, rGridIndexer.GetGeometryTolerance()), mGridIndexer(std::cref(rGridIndexer)),
      mTopology(rGridIndexer)
{ PrepareTopology(); }

void DomainMeshEmbedder::PrepareTopology()
{
    PrepareCellTopology();
    PrepareFaceTopology();
    PrepareFaceVotes();
}

void DomainMeshEmbedder::PrepareCellTopology()
{
    const GridIndexer& r_grid_indexer = mGridIndexer.get();
    MeshPartitioner<CellProduct::Domain> partitioner(mMeshOperator.Query().MeshView(), r_grid_indexer);
    for (const IndexType cell_index : partitioner.GetCellIndices()) {
        auto& r_cell = mTopology.cell_topology[cell_index];
        r_cell.is_trimmed = true;
        r_cell.surface_section.emplace(partitioner.TakeCellProduct(cell_index));
    }
    for (const GridFaceId face : partitioner.GetFaceIds()) {
        PrepareFaceSection(face, partitioner.TakeGridFaceSurface(face));
    }
}

void DomainMeshEmbedder::PrepareFaceSection(GridFaceId Face, TriangleMesh&& rSection)
{
    if (rSection.NumOfTriangles() == 0) { return; }

    const auto& r_grid_indexer = mGridIndexer.get();
    mTopology.face_topology.Get(Face) = detail::GridFaceState::blocked;
    Vector3d weighted_normal{};
    for (const auto& r_triangle : rSection.Triangles<WithNormals>()) {
        weighted_normal += TriangleUtilities::Area(r_triangle) * r_triangle.Normal;
    }
    const double normal_component = weighted_normal[r_grid_indexer.GetFaceAxis(Face)];
    QuESo_ASSERT(
        std::abs(normal_component) > r_grid_indexer.GetGeometryTolerance().ZeroArea(),
        "A face section must have a non-zero normal component."
    );
    const int negative_side_vote = normal_component > 0.0 ? 1 : -1;
    if (r_grid_indexer.GetAdjacentCell(Face, GridFaceSide::negative)) {
        mTopology.face_votes.Get(Face, GridFaceSide::negative) = negative_side_vote;
    }
    if (r_grid_indexer.GetAdjacentCell(Face, GridFaceSide::positive)) {
        mTopology.face_votes.Get(Face, GridFaceSide::positive) = -negative_side_vote;
    }
}

void DomainMeshEmbedder::PrepareFaceTopology()
{
    using Direction = GridIndexer::Direction;
    constexpr std::array forward_directions{ Direction::x_forward, Direction::y_forward, Direction::z_forward };
    const GridIndexer& r_grid_indexer = mGridIndexer.get();
    const auto number_of_cells = static_cast<std::int64_t>(r_grid_indexer.NumberOfElements());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t raw_index = 0; raw_index < number_of_cells; ++raw_index) {
        const IndexType cell_index = static_cast<IndexType>(raw_index);
        for (const Direction direction : forward_directions) {
            if (r_grid_indexer.IsEnd(cell_index, direction)) { continue; }
            const auto [next_index, index_info] = r_grid_indexer.GetNextIndex(cell_index, direction);
            QuESo_ASSERT(index_info == GridIndexer::IndexInfo::middle, "Expected a physical neighboring grid cell.");
            if (mTopology.cell_topology[cell_index].is_trimmed || mTopology.cell_topology[next_index].is_trimmed) {
                mTopology.face_topology.Get(cell_index, direction) = detail::GridFaceState::blocked;
            }
        }
    }
}

void DomainMeshEmbedder::PrepareFaceVotes()
{
    using Direction = GridIndexer::Direction;
    const GridIndexer& r_grid_indexer = mGridIndexer.get();
    const auto number_of_cells = static_cast<std::int64_t>(r_grid_indexer.NumberOfElements());
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (std::int64_t raw_index = 0; raw_index < number_of_cells; ++raw_index) {
        const IndexType cell_index = static_cast<IndexType>(raw_index);
        if (mTopology.cell_topology[cell_index].is_trimmed) { continue; }
        const BoundingBoxType bounds = r_grid_indexer.GetBoundingBoxXYZFromIndex(cell_index);
        const PointType center = 0.5 * (bounds.lower + bounds.upper);
        for (const Direction direction : EnumRange<Direction>()) {
            if (r_grid_indexer.IsEnd(cell_index, direction)
                || mTopology.face_topology.Get(cell_index, direction) == detail::GridFaceState::traversable
                || mTopology.face_votes.Get(cell_index, direction) != 0) {
                continue;
            }
            const auto [next_index, index_info] = r_grid_indexer.GetNextIndex(cell_index, direction);
            QuESo_ASSERT(index_info == GridIndexer::IndexInfo::middle, "Expected a physical neighboring grid cell.");
            const auto& r_neighbor = mTopology.cell_topology[next_index];
            if (!r_neighbor.is_trimmed || !r_neighbor.surface_section) { continue; }
            mTopology.face_votes.Get(cell_index, direction) =
                detail::ClassifyOnBoundedSide(
                    center,
                    r_neighbor.surface_section->SurfaceView(),
                    r_neighbor.surface_section->GetGeometryTolerance().ZeroLength()
                ) == detail::LocalSurfaceClassification::Inside
                    ? 1
                    : -1;
        }
    }
}

ElementStates DomainMeshEmbedder::Classify(FloodFillOptions Options)
{
    if (!mStates) { mStates = detail::FloodFill(mTopology, Options); }
    return *mStates;
}

TrimmedDomain DomainMeshEmbedder::MakeTrimmedDomain(IndexType CellIndex, IndexType MinNumberOfBoundaryTriangles)
{
    QuESo_ASSERT(CellIndex < mTopology.cell_topology.size(), "Cell index is out-of-bounds.");
    auto& r_cell = mTopology.cell_topology[CellIndex];
    QuESo_ASSERT(r_cell.is_trimmed, "Trimmed domains can be constructed only for trimmed cells.");
    QuESo_ASSERT(!r_cell.input_consumed, "The cell surface section has already been consumed.");
    QuESo_ASSERT(mStates.has_value(), "Trimmed-domain construction requires completed classification.");
    QuESo_ASSERT(r_cell.surface_section.has_value(), "Trimmed cell has no prepared surface section.");

    CellSurfaceSection section = std::move(*r_cell.surface_section);
    r_cell.surface_section.reset();
    r_cell.input_consumed = true;
    const auto GlobalIsInside = [this](PointView rPoint) { return mMeshOperator.IsInside(rPoint); };
    TrimmedDomain best_domain(std::move(section), GlobalIsInside, MinNumberOfBoundaryTriangles);
    double best_quality = MeshUtilities::EstimateQuality(best_domain.GetBoundaryMesh());
    constexpr double quality_retry_threshold = 1e-6;
    constexpr std::array<double, 3> expansion_factors{ 1.0, 3.0, 5.0 };
    const GridIndexer& r_grid_indexer = mGridIndexer.get();
    const BoundingBoxType exact_bounds = r_grid_indexer.GetBoundingBoxXYZFromIndex(CellIndex);
    const GeometryTolerance& r_tolerance = r_grid_indexer.GetGeometryTolerance();
    const auto ConsiderExpansion = [&](double factor) {
        BoundingBoxType expanded_bounds = exact_bounds;
        for (IndexType axis = 0; axis < 3; ++axis) {
            expanded_bounds.lower[axis] -= factor * r_tolerance.SnapDistance();
            expanded_bounds.upper[axis] += factor * r_tolerance.SnapDistance();
        }
        const auto candidates = mMeshOperator.Query().GetAabbCandidates(expanded_bounds.lower, expanded_bounds.upper);
        auto expanded_section = mMeshOperator.ClipCellSurfaceSection(candidates, expanded_bounds);
        if (expanded_section.SurfaceView().NumOfTriangles() == 0) { return; }
        TrimmedDomain candidate(std::move(expanded_section), GlobalIsInside, MinNumberOfBoundaryTriangles);
        const double quality = MeshUtilities::EstimateQuality(candidate.GetBoundaryMesh());
        if (quality < best_quality) {
            best_quality = quality;
            best_domain = std::move(candidate);
        }
    };

    if (best_quality < quality_retry_threshold) { return best_domain; }
    for (IndexType i = 0; i < expansion_factors.size(); ++i) {
        ConsiderExpansion(expansion_factors[i]);
        if (best_quality < quality_retry_threshold) { break; }
    }
    return best_domain;
}

}  // namespace queso::embedding
