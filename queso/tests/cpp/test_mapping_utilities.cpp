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

//// STL includes
#include <utility>
#include <vector>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/containers/triangle_proxies.hpp"
#include "queso/utilities/mapping_utilities.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso::testing {

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

}  // namespace

TEST_SUITE("MappingUtilitiesTestSuite")
{

    TEST_CASE("MappingUtilitiesTestSuite - Degenerate Source Bounds Are Rejected When Diagnostics Are Enabled")
    {
        if (!QuESo_DIAGNOSTICS_ENABLED()) { return; }
        for (IndexType axis = 0; axis < 3; ++axis) {
            auto bounds = MakeElementBounds();
            bounds.global.upper[axis] = bounds.global.lower[axis];
            CHECK_THROWS_AS((void)mapping::ToParametric(PointType{}, bounds), Exception);
            bounds = MakeElementBounds();
            bounds.parametric.upper[axis] = bounds.parametric.lower[axis];
            CHECK_THROWS_AS((void)mapping::ToGlobal(PointType{}, bounds), Exception);
            CHECK_THROWS_AS((void)mapping::DetJ(bounds), Exception);
        }
    }

    TEST_CASE("MappingUtilitiesTestSuite - Bounds And Det J")
    {


        constexpr BoundingBoxType bounds_xyz = MakeBox({ -25.0, -111.44, 7.89 }, { 78.67, -35.68, 18.99 });
        constexpr BoundingBoxType bounds_uvw = MakeBox({ -10.0, -2.2, 2.0 }, { 2.0, 10.0, 17.0 });
        constexpr ElementBounds bounds{ bounds_xyz, bounds_uvw };
        constexpr auto r_global_bounds = bounds_xyz;
        constexpr auto r_parametric_bounds = bounds_uvw;

        CheckPointNear(r_global_bounds.lower, bounds_xyz.lower, 1e-12);
        CheckPointNear(r_global_bounds.upper, bounds_xyz.upper, 1e-12);
        CheckPointNear(r_parametric_bounds.lower, bounds_uvw.lower, 1e-12);
        CheckPointNear(r_parametric_bounds.upper, bounds_uvw.upper, 1e-12);
        CheckNear(mapping::DetJ(bounds), ReferenceDetJ(bounds_xyz, bounds_uvw), 1e-12);
    }

    TEST_CASE("MappingUtilitiesTestSuite - Point Mapping")
    {


        constexpr auto r_bounds_xyz = MakeBoundsXYZ();
        constexpr auto r_bounds_uvw = MakeBoundsUVW();
        constexpr auto bounds = MakeElementBounds();

        CheckPointNear(mapping::ToParametric(r_bounds_xyz.lower, bounds), r_bounds_uvw.lower, 1e-12);
        CheckPointNear(mapping::ToParametric(r_bounds_xyz.upper, bounds), r_bounds_uvw.upper, 1e-12);
        CheckPointNear(mapping::ToGlobal(r_bounds_uvw.lower, bounds), r_bounds_xyz.lower, 1e-12);
        CheckPointNear(mapping::ToGlobal(r_bounds_uvw.upper, bounds), r_bounds_xyz.upper, 1e-12);

        constexpr PointType point_global{ 17.25, -79.12, 11.34 };
        const PointType point_parametric = mapping::ToParametric(point_global, bounds);
        CheckPointNear(mapping::ToGlobal(point_parametric, bounds), point_global, 1e-12);

        constexpr PointType point_parametric_2{ -2.5, 4.25, 9.1 };
        const PointType point_global_2 = mapping::ToGlobal(point_parametric_2, bounds);
        CheckPointNear(mapping::ToParametric(point_global_2, bounds), point_parametric_2, 1e-12);
    }

    TEST_CASE("MappingUtilitiesTestSuite - Integration Point Mapping")
    {


        constexpr IntegrationPoint global_ip({ 17.25, -79.12, 11.34 }, 3.7);
        constexpr auto bounds = MakeElementBounds();

        const auto parametric_ip = mapping::ToParametric(global_ip, bounds);
        CheckPointNear(parametric_ip.Point(), mapping::ToParametric(global_ip.Point(), bounds), 1e-12);
        CheckNear(parametric_ip.Weight(), global_ip.Weight() / mapping::DetJ(bounds), 1e-12);

        const auto global_ip_roundtrip = mapping::ToGlobal(parametric_ip, bounds);
        CheckPointNear(global_ip_roundtrip.Point(), global_ip.Point(), 1e-12);
        CheckNear(global_ip_roundtrip.Weight(), global_ip.Weight(), 1e-12);
    }

    TEST_CASE("MappingUtilitiesTestSuite - Boundary Integration Point Round Trip")
    {


        constexpr Vector3d normal{ 2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0 };
        constexpr BoundaryIntegrationPoint global_ip({ 17.25, -79.12, 11.34 }, 3.7, normal);
        constexpr auto bounds = MakeElementBounds();

        const auto parametric_ip = mapping::ToParametric(global_ip, bounds);
        const auto global_ip_roundtrip = mapping::ToGlobal(parametric_ip, bounds);

        CheckPointNear(global_ip_roundtrip.Point(), global_ip.Point(), 1e-12);
        CheckPointNear(global_ip_roundtrip.Normal(), global_ip.Normal(), 1e-12);
        CheckNear(global_ip_roundtrip.Weight(), global_ip.Weight(), 1e-12);
    }

    TEST_CASE("MappingUtilitiesTestSuite - Boundary Integration Point Matches Mapped Triangle")
    {


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

        constexpr auto order = triangle_utilities::TriangleQuadratureOrder::Degree6;
        std::vector<BoundaryIntegrationPoint> global_ips;
        triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(global_triangle, order, [&](auto&& rPoint) {
            global_ips.push_back(std::move(rPoint));
        });
        std::vector<BoundaryIntegrationPoint> parametric_reference_ips;
        triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(parametric_triangle, order, [&](auto&& rPoint) {
            parametric_reference_ips.push_back(std::move(rPoint));
        });

        REQUIRE((global_ips.size()) == (parametric_reference_ips.size()));

        double parametric_area = 0.0;
        for (IndexType i = 0; i < global_ips.size(); ++i) {
            const auto mapped_ip = mapping::ToParametric(global_ips[i], bounds);
            CheckPointNear(mapped_ip.Point(), parametric_reference_ips[i].Point(), 1e-12);
            CheckPointNear(mapped_ip.Normal(), parametric_reference_ips[i].Normal(), 1e-12);
            CheckNear(mapped_ip.Weight(), parametric_reference_ips[i].Weight(), 1e-12);
            parametric_area += mapped_ip.Weight();
        }
        CheckNear(parametric_area, triangle_utilities::Area(parametric_triangle), 1e-12);
    }

}  // TEST_SUITE

}  // namespace queso::testing
