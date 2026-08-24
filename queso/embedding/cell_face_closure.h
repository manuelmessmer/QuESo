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
#include <span>

//// Project includes
#include "queso/embedding/cell_surface_section.h"

namespace queso::embedding::detail {

/// @brief Builds one cell-face closure with the original upper/lower strip algorithm.
/// @details The function owns all temporary state. The classifier is used only for strips without a source crossing and
///          for tolerance-ambiguous coincident upper/lower transitions.
/// @param rSegments Source-surface traces on one cell face.
/// @param rCellBounds Exact original cell bounds.
/// @param Face Face to close.
/// @param rTolerance Authoritative cell geometry tolerance.
/// @param rIsInside Stateless classifier over the owned clipped source surface.
/// @param SwitchAxes Whether to swap the two in-plane projection axes for diagnostic retry.
/// @return Outward-oriented face closure mesh.
[[nodiscard]] TriangleMesh BuildCellFaceClosure(
    std::span<const CellFaceSegment> rSegments,
    const BoundingBoxType& rCellBounds,
    CellFace Face,
    const GeometryTolerance& rTolerance,
    const std::function<bool(PointView)>& rIsInside,
    bool SwitchAxes = false
);

}  // namespace queso::embedding::detail
