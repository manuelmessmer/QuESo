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
#include <span>

//// Project includes
#include "queso/embedding/cell_surface_section.h"
#include "queso/embedding/mesh_query.h"

namespace queso::embedding {

/// @class MeshOperator
/// @brief Provides closed-mesh point classification and clipping over one immutable mesh query.
/// @details The operator owns its query and acceleration structure but not the referenced triangle mesh. The mesh must
///          outlive the operator. Completed trimmed domains are self-contained.
class MeshOperator
{
public:
    ///@name Life Cycle
    ///@{

    /// @brief Constructs a closed-mesh operator and acceleration structure for the supplied mesh.
    /// @param rMesh Triangle mesh view whose referenced data must outlive this operator.
    /// @param rTolerance Immutable geometry tolerance for the query coordinate scale.
    MeshOperator(const TriangleMeshView& rMesh, const GeometryTolerance& rTolerance)
        : mGeometryTolerance(rTolerance), mQuery(rMesh, MeshQueryMode::Closed, mGeometryTolerance)
    {}

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns the immutable mesh query.
    /// @return Mesh query reference.
    [[nodiscard]] const MeshQuery& Query() const noexcept
    { return mQuery; }

    /// @brief Classifies a point against the closed mesh using deterministic ray voting.
    /// @details Points outside the mesh AABB and queries without five conclusive rays classify as outside.
    /// @param rPoint Query point.
    /// @return True when a majority of conclusive rays classify the point as inside.
    [[nodiscard]] bool IsInside(PointView rPoint) const;

    /// @brief Clips supplied candidates to an explicit domain-cell surface section.
    /// @param rCandidateIds Conservative source-triangle candidates.
    /// @param rBounds Exact original cell bounds.
    /// @return Cell surface with oriented perimeter provenance; aligned source triangles are excluded.
    [[nodiscard]] CellSurfaceSection
        ClipCellSurfaceSection(std::span<const IndexType> rCandidateIds, const BoundingBoxType& rBounds) const;

    ///@}

private:
    ///@name Private Members
    ///@{

    GeometryTolerance mGeometryTolerance;  ///< Immutable cell-scale geometry tolerance.
    /// @brief Immutable query and acceleration structure over the non-owning source mesh view.
    MeshQuery mQuery;

    ///@}
};

}  // namespace queso::embedding
