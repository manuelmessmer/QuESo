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

//// Project includes
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/containers/trimmed_element.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/includes/checks.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing::TrimmedElementTestHelpers {

using IntegrationPointType = IntegrationPoint;
using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
using TrimmedElementType = TrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;

/// Embedded box used to create a deterministic trimmed element without external STL data.
constexpr BoundingBoxType MakeEmbeddedBoxBounds()
{ return MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }); }

/// Cell that intersects the generated unit box and contains both active and inactive regions.
constexpr BoundingBoxType MakeCellBoundsXYZ()
{ return MakeBox({ 0.5, 0.0, 0.0 }, { 1.5, 1.0, 1.0 }); }

/// Parametric bounds chosen to give a clean DetJ = (1/2)^3 = 0.125.
constexpr BoundingBoxType MakeCellBoundsUVW()
{ return MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 }); }

constexpr double ReferenceDetJ()
{
    // Delta xyz = {1.0, 1.0, 1.0}, Delta uvw = {2.0, 2.0, 2.0}
    return (1.0 / 2.0) * (1.0 / 2.0) * (1.0 / 2.0);
}

/// Builds a generated box mesh and constructs a TrimmedElement for the pre-selected trimmed cell.
/// Checks that the chosen cell is actually trimmed before constructing the element.
inline TrimmedElementType MakeTrimmedElement()
{
    constexpr auto embedded_box = MakeEmbeddedBoxBounds();
    auto mesh = MeshUtilities::MakeMeshBox(embedded_box.lower, embedded_box.upper);
    BRepOperator brep_op(mesh);

    constexpr auto xyz = MakeCellBoundsXYZ();
    QuESo_CHECK_EQUAL(brep_op.GetIntersectionState(xyz.lower, xyz.upper), IntersectionState::trimmed);

    auto p_domain = brep_op.pGetTrimmedDomain(xyz.lower, xyz.upper, 0.0, 100);
    QuESo_CHECK(p_domain != nullptr);

    return TrimmedElementType(7, ElementBounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() }, std::move(*p_domain));
}

}  // namespace queso::Testing::TrimmedElementTestHelpers
