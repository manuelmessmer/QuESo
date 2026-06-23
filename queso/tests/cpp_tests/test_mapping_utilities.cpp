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
#include "queso/containers/triangle_proxies.hpp"
#include "queso/includes/checks.hpp"
#include "queso/utilities/mapping_utilities.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso::Testing {

namespace {
    Vector3d Normalize(PointView rValue)
    {
        const double norm = Math::Norm(rValue);
        return { rValue[0] / norm, rValue[1] / norm, rValue[2] / norm };
    }

    constexpr BoundingBoxType MakeBoundsXYZ()
    { return MakeBox({ -25.0, -111.44, 7.89 }, { 78.67, -35.68, 18.99 }); }

    constexpr BoundingBoxType MakeBoundsUVW()
    { return MakeBox({ -10.0, -2.2, 2.0 }, { 2.0, 10.0, 17.0 }); }

    constexpr ElementBounds MakeElementBounds()
    { return ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() }; }

    constexpr double ReferenceDetJ(const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW)
    {
        return ((rBoundsXYZ.upper[0] - rBoundsXYZ.lower[0]) / (rBoundsUVW.upper[0] - rBoundsUVW.lower[0]))
               * ((rBoundsXYZ.upper[1] - rBoundsXYZ.lower[1]) / (rBoundsUVW.upper[1] - rBoundsUVW.lower[1]))
               * ((rBoundsXYZ.upper[2] - rBoundsXYZ.lower[2]) / (rBoundsUVW.upper[2] - rBoundsUVW.lower[2]));
    }

}// namespace

BOOST_AUTO_TEST_SUITE(MappingUtilitiesTestSuite)

BOOST_AUTO_TEST_CASE(BoundsAndDetJ)
{
    QuESo_INFO << "Testing :: MappingUtilities :: BoundsAndDetJ" << std::endl;

    constexpr BoundingBoxType bounds_xyz = MakeBox({ -25.0, -111.44, 7.89 }, { 78.67, -35.68, 18.99 });
    constexpr BoundingBoxType bounds_uvw = MakeBox({ -10.0, -2.2, 2.0 }, { 2.0, 10.0, 17.0 });
    constexpr ElementBounds bounds{ bounds_xyz, bounds_uvw };
    constexpr auto r_global_bounds = bounds_xyz;
    constexpr auto r_parametric_bounds = bounds_uvw;

    QuESo_CHECK_POINT_NEAR(r_global_bounds.lower, bounds_xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(r_global_bounds.upper, bounds_xyz.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(r_parametric_bounds.lower, bounds_uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(r_parametric_bounds.upper, bounds_uvw.upper, 1e-12);
    QuESo_CHECK_NEAR(mapping::DetJ(bounds), ReferenceDetJ(bounds_xyz, bounds_uvw), 1e-12);
}

BOOST_AUTO_TEST_CASE(PointMapping)
{
    QuESo_INFO << "Testing :: MappingUtilities :: PointMapping" << std::endl;

    constexpr auto r_bounds_xyz = MakeBoundsXYZ();
    constexpr auto r_bounds_uvw = MakeBoundsUVW();
    constexpr auto bounds = MakeElementBounds();

    QuESo_CHECK_POINT_NEAR(mapping::ToParametric(r_bounds_xyz.lower, bounds), r_bounds_uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(mapping::ToParametric(r_bounds_xyz.upper, bounds), r_bounds_uvw.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(mapping::ToGlobal(r_bounds_uvw.lower, bounds), r_bounds_xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(mapping::ToGlobal(r_bounds_uvw.upper, bounds), r_bounds_xyz.upper, 1e-12);

    constexpr PointType point_global{ 17.25, -79.12, 11.34 };
    const PointType point_parametric = mapping::ToParametric(point_global, bounds);
    QuESo_CHECK_POINT_NEAR(mapping::ToGlobal(point_parametric, bounds), point_global, 1e-12);

    constexpr PointType point_parametric_2{ -2.5, 4.25, 9.1 };
    const PointType point_global_2 = mapping::ToGlobal(point_parametric_2, bounds);
    QuESo_CHECK_POINT_NEAR(mapping::ToParametric(point_global_2, bounds), point_parametric_2, 1e-12);
}

BOOST_AUTO_TEST_CASE(IntegrationPointMapping)
{
    QuESo_INFO << "Testing :: MappingUtilities :: IntegrationPointMapping" << std::endl;

    constexpr IntegrationPoint global_ip({ 17.25, -79.12, 11.34 }, 3.7);
    constexpr auto bounds = MakeElementBounds();

    const auto parametric_ip = mapping::ToParametric(global_ip, bounds);
    QuESo_CHECK_POINT_NEAR(parametric_ip.Point(), mapping::ToParametric(global_ip.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(parametric_ip.Weight(), global_ip.Weight() / mapping::DetJ(bounds), 1e-12);

    const auto global_ip_roundtrip = mapping::ToGlobal(parametric_ip, bounds);
    QuESo_CHECK_POINT_NEAR(global_ip_roundtrip.Point(), global_ip.Point(), 1e-12);
    QuESo_CHECK_NEAR(global_ip_roundtrip.Weight(), global_ip.Weight(), 1e-12);
}

BOOST_AUTO_TEST_CASE(BoundaryIntegrationPointRoundTrip)
{
    QuESo_INFO << "Testing :: MappingUtilities :: BoundaryIntegrationPointRoundTrip" << std::endl;

    constexpr Vector3d normal{ 2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0 };
    constexpr BoundaryIntegrationPoint global_ip({ 17.25, -79.12, 11.34 }, 3.7, normal);
    constexpr auto bounds = MakeElementBounds();

    const auto parametric_ip = mapping::ToParametric(global_ip, bounds);
    const auto global_ip_roundtrip = mapping::ToGlobal(parametric_ip, bounds);

    QuESo_CHECK_POINT_NEAR(global_ip_roundtrip.Point(), global_ip.Point(), 1e-12);
    QuESo_CHECK_POINT_NEAR(global_ip_roundtrip.Normal(), global_ip.Normal(), 1e-12);
    QuESo_CHECK_NEAR(global_ip_roundtrip.Weight(), global_ip.Weight(), 1e-12);
}

BOOST_AUTO_TEST_CASE(BoundaryIntegrationPointMatchesMappedTriangle)
{
    QuESo_INFO << "Testing :: MappingUtilities :: BoundaryIntegrationPointMatchesMappedTriangle" << std::endl;

    constexpr auto bounds = MakeElementBounds();
    constexpr PointType p1{ -25.0, -111.44, 7.89 };
    constexpr PointType p2{ 40.0, -91.0, 9.2 };
    constexpr PointType p3{ 5.0, -50.0, 17.4 };

    const Vector3d global_normal = Normalize(Math::Cross(p2 - p1, p3 - p1));
    const TriangleProxy<WithNormals> global_triangle{ p1, p2, p3, global_normal };

    const PointType p1_parametric = mapping::ToParametric(p1, bounds);
    const PointType p2_parametric = mapping::ToParametric(p2, bounds);
    const PointType p3_parametric = mapping::ToParametric(p3, bounds);
    const Vector3d parametric_normal =
        Normalize(Math::Cross(p2_parametric - p1_parametric, p3_parametric - p1_parametric));
    const TriangleProxy<WithNormals> parametric_triangle{
        p1_parametric, p2_parametric, p3_parametric, parametric_normal
    };

    constexpr IndexType method = 3;
    const auto global_ips = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPoint>(global_triangle, method);
    const auto parametric_reference_ips =
        TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPoint>(parametric_triangle, method);

    QuESo_CHECK_EQUAL(global_ips.size(), parametric_reference_ips.size());

    double parametric_area = 0.0;
    for (IndexType i = 0; i < global_ips.size(); ++i) {
        const auto mapped_ip = mapping::ToParametric(global_ips[i], bounds);
        QuESo_CHECK_POINT_NEAR(mapped_ip.Point(), parametric_reference_ips[i].Point(), 1e-12);
        QuESo_CHECK_POINT_NEAR(mapped_ip.Normal(), parametric_reference_ips[i].Normal(), 1e-12);
        QuESo_CHECK_NEAR(mapped_ip.Weight(), parametric_reference_ips[i].Weight(), 1e-12);
        parametric_area += mapped_ip.Weight();
    }
    QuESo_CHECK_NEAR(parametric_area, TriangleUtilities::Area(parametric_triangle), 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

}// End namespace queso::Testing
