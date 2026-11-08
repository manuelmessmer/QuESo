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

//// External includes
#include <doctest/doctest.h>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/tests/cpp/helpers/trimmed_element_test_helpers.hpp"
#include "queso/utilities/mesh_utilities.h"

// This suite tests the TrimmedElement container API on a representative trimmed cell. Related coverage:
// test_moment_fitting_compute.cpp checks public moment-fitting Compute behavior on trimmed elements, and
// test_point_elimination.cpp checks reduced rules on larger trimmed geometries.

namespace queso::testing {

namespace {

    using namespace TrimmedElementTestHelpers;

    double SumWeights(const std::vector<BoundaryIntegrationPointType>& rIps)
    {
        double sum = 0.0;
        for (const auto& rIp : rIps) { sum += rIp.Weight(); }
        return sum;
    }

}  // namespace

TEST_SUITE("TrimmedElementTestSuite")
{

    TEST_CASE("TrimmedElementTestSuite - Construction And Defaults")
    {
        // Verifies basic construction state for a representative trimmed element.


        const auto element = MakeTrimmedElement();

        CHECK((element.GetId()) == (7UL));
        CHECK_FALSE(!element.IsTrimmed());
    }

    TEST_CASE("TrimmedElementTestSuite - Cell Bounds Vs Active Domain Bounds")
    {
        // Verifies cell bounds remain exact and active-domain bounds are mapped consistently.


        const auto element = MakeTrimmedElement();
        constexpr auto xyz = MakeCellBoundsXYZ();
        constexpr auto uvw = MakeCellBoundsUVW();

        // Cell bounds must equal the construction bounds exactly.
        CheckPointNear(element.GetCellBounds<CoordinateSpace::global>().lower, xyz.lower, 1e-12);
        CheckPointNear(element.GetCellBounds<CoordinateSpace::global>().upper, xyz.upper, 1e-12);
        CheckPointNear(element.GetCellBounds<CoordinateSpace::parametric>().lower, uvw.lower, 1e-12);
        CheckPointNear(element.GetCellBounds<CoordinateSpace::parametric>().upper, uvw.upper, 1e-12);

        // Active-domain bounds (global) must be contained within the cell bounds.
        const auto active_global = element.GetActiveDomainBounds<CoordinateSpace::global>();
        for (IndexType i = 0; i < 3; ++i) {
            CHECK_FALSE(active_global.lower[i] < xyz.lower[i] - 1e-10);
            CHECK_FALSE(active_global.upper[i] > xyz.upper[i] + 1e-10);
        }
        const auto active_parametric = element.GetActiveDomainBounds<CoordinateSpace::parametric>();
        constexpr ElementBounds bounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() };
        CheckPointNear(active_global.lower, mapping::ToGlobal(active_parametric.lower, bounds), 1e-12);
        CheckPointNear(active_global.upper, mapping::ToGlobal(active_parametric.upper, bounds), 1e-12);
    }

    TEST_CASE("TrimmedElementTestSuite - Det J")
    {
        // Verifies the element Jacobian determinant for the chosen affine cell mapping.


        const auto element = MakeTrimmedElement();
        CheckNear(element.DetJ(), ReferenceDetJ(), 1e-12);
    }

    TEST_CASE("TrimmedElementTestSuite - Integration Points")
    {
        // Verifies storage and global mapping/scaling of element integration points.


        auto element = MakeTrimmedElement();
        constexpr IntegrationPointType ip_a({ -0.5, 0.0, 0.5 }, 1.5);
        constexpr IntegrationPointType ip_b({ 0.3, -0.7, 0.1 }, 0.75);

        element.GetIntegrationPoints().push_back(ip_a);
        element.GetIntegrationPoints().push_back(ip_b);

        // Parametric IPs are returned as-is.
        const auto& r_param = static_cast<const TrimmedElementType&>(element).GetIntegrationPoints();
        REQUIRE((r_param.size()) == (2UL));
        CheckPointNear(r_param[0].Point(), ip_a.Point(), 1e-12);
        CheckPointNear(r_param[1].Point(), ip_b.Point(), 1e-12);

        // Global IPs must be mapped by the element-local mapper and weights scaled by DetJ.
        const auto global_ips =
            static_cast<const TrimmedElementType&>(element).GetIntegrationPoints<CoordinateSpace::global>();
        auto it = global_ips.begin();
        const auto global_a = *it++;
        const auto global_b = *it;

        const double det_j = element.DetJ();
        constexpr ElementBounds bounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() };

        CheckPointNear(global_a.Point(), mapping::ToGlobal(ip_a.Point(), bounds), 1e-12);
        CheckNear(global_a.Weight(), ip_a.Weight() * det_j, 1e-12);
        CheckPointNear(global_b.Point(), mapping::ToGlobal(ip_b.Point(), bounds), 1e-12);
        CheckNear(global_b.Weight(), ip_b.Weight() * det_j, 1e-12);
    }

    TEST_CASE("TrimmedElementTestSuite - Boundary Mesh And Boundary Ips")
    {
        // Verifies boundary mesh and boundary integration points are available in both coordinate spaces.


        const auto element = MakeTrimmedElement();

        // Trimmed boundary mesh must be non-empty and have positive area.
        const auto boundary_mesh = element.GetActiveDomainBoundaryMesh();
        CHECK_FALSE(!(boundary_mesh.NumOfTriangles() > 0UL));
        CHECK_FALSE(!(MeshUtilities::Area(boundary_mesh) > 0.0));

        // Global boundary IPs must carry positive total weight.
        const auto bips_global =
            element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::global>();
        CHECK_FALSE(bips_global.empty());
        CHECK_FALSE(!(SumWeights(bips_global) > 0.0));

        // Parametric boundary IPs must also have positive total weight.
        const auto bips_parametric =
            element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::parametric>();
        CHECK_FALSE(bips_parametric.empty());
        CHECK_FALSE(!(SumWeights(bips_parametric) > 0.0));
    }

    TEST_CASE("TrimmedElementTestSuite - Is Inside Active Domain")
    {
        // Verifies active-domain point classification for representative inside/outside points.


        const auto element = MakeTrimmedElement();

        // Inside the embedded unit box and inside the trimmed cell.
        CHECK_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>({ 0.75, 0.5, 0.5 }));

        // Outside the embedded unit box but still inside the trimmed cell.
        CHECK_FALSE(element.IsInsideActiveDomain<CoordinateSpace::global>({ 1.25, 0.5, 0.5 }));
    }

}  // TEST_SUITE

}  // namespace queso::testing
