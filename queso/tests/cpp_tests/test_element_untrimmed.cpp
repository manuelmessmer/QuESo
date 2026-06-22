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
#include "queso/containers/untrimmed_element.hpp"
#include "queso/includes/checks.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {

namespace {

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using ElementType = UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;

    constexpr BoundingBoxType MakeBoundsXYZ()
    { return MakeBox({ -25.0, -111.44, 7.89 }, { 78.67, -35.68, 18.99 }); }

    constexpr BoundingBoxType MakeBoundsUVW()
    { return MakeBox({ -10.0, -2.2, 2.0 }, { 2.0, 10.0, 17.0 }); }

    constexpr double ReferenceDetJ(const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW)
    {
        return ((rBoundsXYZ.upper[0] - rBoundsXYZ.lower[0]) / (rBoundsUVW.upper[0] - rBoundsUVW.lower[0]))
               * ((rBoundsXYZ.upper[1] - rBoundsXYZ.lower[1]) / (rBoundsUVW.upper[1] - rBoundsUVW.lower[1]))
               * ((rBoundsXYZ.upper[2] - rBoundsXYZ.lower[2]) / (rBoundsUVW.upper[2] - rBoundsUVW.lower[2]));
    }

    constexpr double BoxSurfaceArea(const BoundingBoxType& rBounds)
    {
        const double dx = rBounds.upper[0] - rBounds.lower[0];
        const double dy = rBounds.upper[1] - rBounds.lower[1];
        const double dz = rBounds.upper[2] - rBounds.lower[2];
        return 2.0 * (dx * dy + dx * dz + dy * dz);
    }

    ElementType MakeElement()
    { return ElementType(17, ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() }); }

    double SumWeights(const std::vector<BoundaryIntegrationPointType>& rIps)
    {
        double sum = 0.0;
        for (const auto& rIp : rIps) { sum += rIp.Weight(); }
        return sum;
    }

}// namespace

BOOST_AUTO_TEST_SUITE(UntrimmedElementTestSuite)

BOOST_AUTO_TEST_CASE(ConstructionAndDefaults)
{
    QuESo_INFO << "Testing :: UntrimmedElement :: ConstructionAndDefaults" << std::endl;

    const auto element = MakeElement();

    QuESo_CHECK_EQUAL(element.GetId(), 17UL);
    QuESo_CHECK_IS_FALSE(element.IsTrimmed());
    QuESo_CHECK_EQUAL(element.GetValue<bool>(ElementValues::is_visited), false);
    QuESo_CHECK_NEAR(element.GetValue<double>(ElementValues::neighbor_coefficient), 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(CellBoundsAndActiveDomainBounds)
{
    QuESo_INFO << "Testing :: UntrimmedElement :: CellBoundsAndActiveDomainBounds" << std::endl;

    const auto element = MakeElement();
    constexpr auto bounds_xyz = MakeBoundsXYZ();
    constexpr auto bounds_uvw = MakeBoundsUVW();

    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::global>().lower, bounds_xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::global>().upper, bounds_xyz.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::parametric>().lower, bounds_uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetCellBounds<CoordinateSpace::parametric>().upper, bounds_uvw.upper, 1e-12);

    QuESo_CHECK_POINT_NEAR(element.GetActiveDomainBounds<CoordinateSpace::global>().lower, bounds_xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetActiveDomainBounds<CoordinateSpace::global>().upper, bounds_xyz.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetActiveDomainBounds<CoordinateSpace::parametric>().lower, bounds_uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(element.GetActiveDomainBounds<CoordinateSpace::parametric>().upper, bounds_uvw.upper, 1e-12);
}

BOOST_AUTO_TEST_CASE(IntegrationPoints)
{
    QuESo_INFO << "Testing :: UntrimmedElement :: IntegrationPoints" << std::endl;

    auto element = MakeElement();
    constexpr IntegrationPointType ip_a({ -2.5, 4.25, 9.1 }, 1.5);
    constexpr IntegrationPointType ip_b({ 1.0, -1.0, 2.0 }, 0.75);

    auto& r_integration_points = element.GetIntegrationPoints();
    r_integration_points.push_back(ip_a);
    r_integration_points.push_back(ip_b);

    const auto& r_parametric_points = static_cast<const ElementType&>(element).GetIntegrationPoints();
    QuESo_CHECK_EQUAL(r_parametric_points.size(), 2UL);
    QuESo_CHECK_POINT_NEAR(r_parametric_points[0].Point(), ip_a.Point(), 1e-12);
    QuESo_CHECK_POINT_NEAR(r_parametric_points[1].Point(), ip_b.Point(), 1e-12);

    const auto global_points = static_cast<const ElementType&>(element).GetIntegrationPoints<CoordinateSpace::global>();
    auto it = global_points.begin();

    const auto global_ip_a = *it++;
    const auto global_ip_b = *it;
    const double det_j = element.DetJ();
    constexpr ElementBounds bounds{ MakeBoundsXYZ(), MakeBoundsUVW() };

    QuESo_CHECK_POINT_NEAR(global_ip_a.Point(), mapping::ToGlobal(ip_a.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(global_ip_a.Weight(), ip_a.Weight() * det_j, 1e-12);
    QuESo_CHECK_POINT_NEAR(global_ip_b.Point(), mapping::ToGlobal(ip_b.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(global_ip_b.Weight(), ip_b.Weight() * det_j, 1e-12);
}

BOOST_AUTO_TEST_CASE(BoundaryMeshAndBoundaryIps)
{
    QuESo_INFO << "Testing :: UntrimmedElement :: BoundaryMeshAndBoundaryIps" << std::endl;

    const auto element = MakeElement();
    constexpr auto bounds_xyz = MakeBoundsXYZ();
    constexpr auto bounds_uvw = MakeBoundsUVW();

    const auto boundary_mesh = element.GetActiveDomainBoundaryMesh();
    QuESo_CHECK_IS_FALSE(!(boundary_mesh.NumOfTriangles() > 0UL));
    QuESo_CHECK_NEAR(MeshUtilities::Area(boundary_mesh), BoxSurfaceArea(bounds_xyz), 1e-12);

    const auto boundary_ips_global =
        element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::global>();
    QuESo_CHECK_NEAR(SumWeights(boundary_ips_global), BoxSurfaceArea(bounds_xyz), 1e-10);

    const auto boundary_ips_parametric =
        element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::parametric>();
    QuESo_CHECK_NEAR(SumWeights(boundary_ips_parametric), BoxSurfaceArea(bounds_uvw), 1e-10);
}

BOOST_AUTO_TEST_CASE(InsideAndIntersectionState)
{
    QuESo_INFO << "Testing :: UntrimmedElement :: InsideAndIntersectionState" << std::endl;

    const auto element = MakeElement();

    QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>(PointType{ 0.0, -80.0, 10.0 }));
    QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>(PointType{ -25.0, -111.44, 7.89 }));
    QuESo_CHECK_IS_FALSE(element.IsInsideActiveDomain<CoordinateSpace::global>(PointType{ -26.0, -80.0, 10.0 }));

    QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::parametric>(PointType{ -2.0, 4.0, 9.0 }));
    QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::parametric>(PointType{ -10.0, -2.2, 2.0 }));
    QuESo_CHECK_IS_FALSE(element.IsInsideActiveDomain<CoordinateSpace::parametric>(PointType{ -10.1, 4.0, 9.0 }));

    QuESo_CHECK_EQUAL(
        element.GetIntersectionState<CoordinateSpace::global>(
            PointType{ -20.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }, 0.0
        ),
        queso::IntersectionState::inside
    );
    QuESo_CHECK_EQUAL(
        element.GetIntersectionState<CoordinateSpace::global>(
            PointType{ 100.0, -100.0, 8.5 }, PointType{ 110.0, -90.0, 12.0 }, 0.0
        ),
        queso::IntersectionState::outside
    );
    QuESo_CHECK_EQUAL(
        element.GetIntersectionState<CoordinateSpace::global>(
            PointType{ -30.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }, 0.0
        ),
        queso::IntersectionState::trimmed
    );

    QuESo_CHECK_EQUAL(
        element.GetIntersectionState<CoordinateSpace::parametric>(
            PointType{ -8.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }, 0.0
        ),
        queso::IntersectionState::inside
    );
    QuESo_CHECK_EQUAL(
        element.GetIntersectionState<CoordinateSpace::parametric>(
            PointType{ 3.0, -1.0, 3.0 }, PointType{ 4.0, 8.0, 14.0 }, 0.0
        ),
        queso::IntersectionState::outside
    );
    QuESo_CHECK_EQUAL(
        element.GetIntersectionState<CoordinateSpace::parametric>(
            PointType{ -11.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }, 0.0
        ),
        queso::IntersectionState::trimmed
    );
}

BOOST_AUTO_TEST_CASE(DetJAndDataSet)
{
    QuESo_INFO << "Testing :: UntrimmedElement :: DetJAndDataSet" << std::endl;

    auto element = MakeElement();
    constexpr auto bounds_xyz = MakeBoundsXYZ();
    constexpr auto bounds_uvw = MakeBoundsUVW();

    QuESo_CHECK_NEAR(element.DetJ(), ReferenceDetJ(bounds_xyz, bounds_uvw), 1e-12);

    element.SetValue(ElementValues::is_visited, true);
    element.SetValue(ElementValues::neighbor_coefficient, 2.5);

    QuESo_CHECK_EQUAL(element.GetValue<bool>(ElementValues::is_visited), true);
    QuESo_CHECK_NEAR(element.GetValue<double>(ElementValues::neighbor_coefficient), 2.5, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

}// End namespace queso::Testing
