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
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/tests/cpp_tests/trimmed_element_test_helpers.hpp"
#include "queso/utilities/mesh_utilities.h"

// This suite tests the TrimmedElement container API on a representative trimmed cell. Related coverage:
// test_moment_fitting_compute.cpp checks public moment-fitting Compute behavior on trimmed elements, and
// test_point_elimination.cpp checks reduced rules on larger trimmed geometries.

namespace queso::Testing {

namespace {

    using namespace TrimmedElementTestHelpers;

    double SumWeights(const std::vector<BoundaryIntegrationPointType>& rIps)
    {
        double sum = 0.0;
        for (const auto& rIp : rIps) { sum += rIp.Weight(); }
        return sum;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(TrimmedElementTestSuite)

BOOST_AUTO_TEST_CASE(ConstructionAndDefaults)
{
    // Verifies basic construction state for a representative trimmed element.
    QuESo_INFO << "Testing :: TrimmedElement :: ConstructionAndDefaults" << std::endl;

    const auto element = MakeTrimmedElement();

    QuESo_CHECK_EQUAL(element.GetId(), 7UL);
    QuESo_CHECK_IS_FALSE(!element.IsTrimmed());
}

BOOST_AUTO_TEST_CASE(CellBoundsVsActiveDomainBounds)
{
    // Verifies cell bounds remain exact and active-domain bounds are mapped consistently.
    QuESo_INFO << "Testing :: TrimmedElement :: CellBoundsVsActiveDomainBounds" << std::endl;

    const auto element = MakeTrimmedElement();
    constexpr auto xyz = MakeCellBoundsXYZ();
    constexpr auto uvw = MakeCellBoundsUVW();

    // Cell bounds must equal the construction bounds exactly.
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::global>().lower, xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::global>().upper, xyz.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::parametric>().lower, uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::parametric>().upper, uvw.upper, 1e-12);

    // Active-domain bounds (global) must be contained within the cell bounds.
    const auto active_global = element.GetActiveDomainBounds<CoordinateSpace::global>();
    for (IndexType i = 0; i < 3; ++i) {
        QuESo_CHECK_IS_FALSE(active_global.lower[i] < xyz.lower[i] - 1e-10);
        QuESo_CHECK_IS_FALSE(active_global.upper[i] > xyz.upper[i] + 1e-10);
    }
    const auto active_parametric = element.GetActiveDomainBounds<CoordinateSpace::parametric>();
    constexpr ElementBounds bounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() };
    QuESo_CHECK_POINT_NEAR(active_global.lower, mapping::ToGlobal(active_parametric.lower, bounds), 1e-12);
    QuESo_CHECK_POINT_NEAR(active_global.upper, mapping::ToGlobal(active_parametric.upper, bounds), 1e-12);
}

BOOST_AUTO_TEST_CASE(DetJ)
{
    // Verifies the element Jacobian determinant for the chosen affine cell mapping.
    QuESo_INFO << "Testing :: TrimmedElement :: DetJ" << std::endl;

    const auto element = MakeTrimmedElement();
    QuESo_CHECK_NEAR(element.DetJ(), ReferenceDetJ(), 1e-12);
}

BOOST_AUTO_TEST_CASE(IntegrationPoints)
{
    // Verifies storage and global mapping/scaling of element integration points.
    QuESo_INFO << "Testing :: TrimmedElement :: IntegrationPoints" << std::endl;

    auto element = MakeTrimmedElement();
    constexpr IntegrationPointType ip_a({ -0.5, 0.0, 0.5 }, 1.5);
    constexpr IntegrationPointType ip_b({ 0.3, -0.7, 0.1 }, 0.75);

    element.GetIntegrationPoints().push_back(ip_a);
    element.GetIntegrationPoints().push_back(ip_b);

    // Parametric IPs are returned as-is.
    const auto& r_param = static_cast<const TrimmedElementType&>(element).GetIntegrationPoints();
    QuESo_CHECK_EQUAL(r_param.size(), 2UL);
    QuESo_CHECK_POINT_NEAR(r_param[0].Point(), ip_a.Point(), 1e-12);
    QuESo_CHECK_POINT_NEAR(r_param[1].Point(), ip_b.Point(), 1e-12);

    // Global IPs must be mapped by the element-local mapper and weights scaled by DetJ.
    const auto global_ips =
        static_cast<const TrimmedElementType&>(element).GetIntegrationPoints<CoordinateSpace::global>();
    auto it = global_ips.begin();
    const auto global_a = *it++;
    const auto global_b = *it;

    const double det_j = element.DetJ();
    constexpr ElementBounds bounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() };

    QuESo_CHECK_POINT_NEAR(global_a.Point(), mapping::ToGlobal(ip_a.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(global_a.Weight(), ip_a.Weight() * det_j, 1e-12);
    QuESo_CHECK_POINT_NEAR(global_b.Point(), mapping::ToGlobal(ip_b.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(global_b.Weight(), ip_b.Weight() * det_j, 1e-12);
}

BOOST_AUTO_TEST_CASE(BoundaryMeshAndBoundaryIps)
{
    // Verifies boundary mesh and boundary integration points are available in both coordinate spaces.
    QuESo_INFO << "Testing :: TrimmedElement :: BoundaryMeshAndBoundaryIps" << std::endl;

    const auto element = MakeTrimmedElement();

    // Trimmed boundary mesh must be non-empty and have positive area.
    const auto boundary_mesh = element.GetActiveDomainBoundaryMesh();
    QuESo_CHECK_IS_FALSE(!(boundary_mesh.NumOfTriangles() > 0UL));
    QuESo_CHECK_IS_FALSE(!(MeshUtilities::Area(boundary_mesh) > 0.0));

    // Global boundary IPs must carry positive total weight.
    const auto bips_global =
        element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::global>();
    QuESo_CHECK_IS_FALSE(bips_global.empty());
    QuESo_CHECK_IS_FALSE(!(SumWeights(bips_global) > 0.0));

    // Parametric boundary IPs must also have positive total weight.
    const auto bips_parametric =
        element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::parametric>();
    QuESo_CHECK_IS_FALSE(bips_parametric.empty());
    QuESo_CHECK_IS_FALSE(!(SumWeights(bips_parametric) > 0.0));
}

BOOST_AUTO_TEST_CASE(IsInsideActiveDomain)
{
    // Verifies active-domain point classification for representative inside/outside points.
    QuESo_INFO << "Testing :: TrimmedElement :: IsInsideActiveDomain" << std::endl;

    const auto element = MakeTrimmedElement();

    // Inside the embedded unit box and inside the trimmed cell.
    QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>({ 0.75, 0.5, 0.5 }));

    // Outside the embedded unit box but still inside the trimmed cell.
    QuESo_CHECK_IS_FALSE(element.IsInsideActiveDomain<CoordinateSpace::global>({ 1.25, 0.5, 0.5 }));
}

BOOST_AUTO_TEST_SUITE_END()

}  // End namespace queso::Testing
