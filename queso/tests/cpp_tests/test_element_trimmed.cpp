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
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/containers/trimmed_element.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/includes/checks.hpp"
#include "queso/io/io_utilities.h"
#include "queso/tests/cpp_tests/global_config.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {

namespace {

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using TrimmedElementType = TrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;

    /// Cell that intersects the cylinder surface (cylinder radius ~1, centered at origin).
    /// Corners at r = sqrt(0.5) ≈ 0.71 (inside) and r = sqrt(2.5) ≈ 1.58 (outside).
    constexpr BoundingBoxType MakeCellBoundsXYZ()
    { return MakeBox({ 0.5, -0.5, 0.0 }, { 1.5, 0.5, 1.0 }); }

    /// Parametric bounds chosen to give a clean DetJ = (1/2)^3 = 0.125.
    constexpr BoundingBoxType MakeCellBoundsUVW()
    { return MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 }); }

    constexpr double ReferenceDetJ()
    {
        // Δxyz = {1.0, 1.0, 1.0}, Δuvw = {2.0, 2.0, 2.0}
        return (1.0 / 2.0) * (1.0 / 2.0) * (1.0 / 2.0);
    }

    /// Reads cylinder.stl and constructs a TrimmedElement for the pre-selected trimmed cell.
    /// Requires that the chosen cell is actually trimmed — asserted via BOOST_REQUIRE.
    TrimmedElementType MakeTrimmedElement()
    {
        const std::string stl_path = GlobalConfig::GetInstance().BaseDir + "/data/cylinder.stl";
        TriangleMesh mesh{};
        IO::ReadMeshFromSTL(mesh, stl_path.c_str());
        BRepOperator brep_op(mesh);

        constexpr auto xyz = MakeCellBoundsXYZ();
        BOOST_REQUIRE_EQUAL(brep_op.GetIntersectionState(xyz.lower, xyz.upper), IntersectionState::trimmed);

        auto p_domain = brep_op.pGetTrimmedDomain(xyz.lower, xyz.upper, 0.0, 100);
        BOOST_REQUIRE(p_domain != nullptr);

        return TrimmedElementType(7, ElementBounds{ MakeCellBoundsXYZ(), MakeCellBoundsUVW() }, std::move(*p_domain));
    }

    double SumWeights(const std::vector<BoundaryIntegrationPointType>& rIps)
    {
        double sum = 0.0;
        for (const auto& rIp : rIps) { sum += rIp.Weight(); }
        return sum;
    }

}// namespace

BOOST_AUTO_TEST_SUITE(TrimmedElementTestSuite)

BOOST_AUTO_TEST_CASE(ConstructionAndDefaults)
{
    QuESo_INFO << "Testing :: TrimmedElement :: ConstructionAndDefaults" << std::endl;

    const auto element = MakeTrimmedElement();

    QuESo_CHECK_EQUAL(element.GetId(), 7UL);
    QuESo_CHECK_IS_FALSE(!element.IsTrimmed());
}

BOOST_AUTO_TEST_CASE(CellBoundsVsActiveDomainBounds)
{
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
    QuESo_INFO << "Testing :: TrimmedElement :: DetJ" << std::endl;

    const auto element = MakeTrimmedElement();
    QuESo_CHECK_NEAR(element.DetJ(), ReferenceDetJ(), 1e-12);
}

BOOST_AUTO_TEST_CASE(IntegrationPoints)
{
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
    QuESo_INFO << "Testing :: TrimmedElement :: IsInsideActiveDomain" << std::endl;

    const auto element = MakeTrimmedElement();

    // Well inside the cylinder (r = 0.5 < 1) and inside the cell.
    QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>({ 0.5, 0.0, 0.5 }));

    // Outside the cylinder (r = 1.3 > 1) but still inside the cell.
    QuESo_CHECK_IS_FALSE(element.IsInsideActiveDomain<CoordinateSpace::global>({ 1.3, 0.0, 0.5 }));
}

BOOST_AUTO_TEST_SUITE_END()

}// End namespace queso::Testing
