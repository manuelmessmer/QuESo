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
#include "queso/embedding/mesh_operator.h"
#include "queso/embedding/trimmed_domain.h"
#include "queso/includes/checks.hpp"
#include "queso/tests/cpp_tests/trimmed_domain_test_helpers.hpp"
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
    static const TriangleMesh mesh = MeshUtilities::MakeMeshBox(embedded_box.lower, embedded_box.upper);
    constexpr auto xyz = MakeCellBoundsXYZ();
    const embedding::MeshOperator mesh_operator(mesh.View(), TrimmedDomainTestHelpers::MakeTolerance(xyz));
    auto p_domain = TrimmedDomainTestHelpers::MakeTrimmedDomain(mesh_operator, xyz, 100);
    QuESo_CHECK(p_domain != nullptr);

    return TrimmedElementType(7, ElementBounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() }, std::move(*p_domain));
}

}  // namespace queso::Testing::TrimmedElementTestHelpers
