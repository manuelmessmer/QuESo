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
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/containers/untrimmed_element.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::testing {

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

}  // namespace

TEST_SUITE("UntrimmedElementTestSuite")
{

    TEST_CASE("UntrimmedElementTestSuite - Construction And Defaults")
    {


        const auto element = MakeElement();

        CHECK((element.GetId()) == (17UL));
        CHECK_FALSE(element.IsTrimmed());
        CHECK((element.GetValue<bool>(ElementValues::is_visited)) == (false));
        CheckNear(element.GetValue<double>(ElementValues::neighbor_coefficient), 1.0, 1e-12);
    }

    TEST_CASE("UntrimmedElementTestSuite - Cell Bounds And Active Domain Bounds")
    {


        const auto element = MakeElement();
        constexpr auto bounds_xyz = MakeBoundsXYZ();
        constexpr auto bounds_uvw = MakeBoundsUVW();

        CheckPointNear(element.GetCellBounds<CoordinateSpace::global>().lower, bounds_xyz.lower, 1e-12);
        CheckPointNear(element.GetCellBounds<CoordinateSpace::global>().upper, bounds_xyz.upper, 1e-12);
        CheckPointNear(element.GetCellBounds<CoordinateSpace::parametric>().lower, bounds_uvw.lower, 1e-12);
        CheckPointNear(element.GetCellBounds<CoordinateSpace::parametric>().upper, bounds_uvw.upper, 1e-12);

        CheckPointNear(element.GetActiveDomainBounds<CoordinateSpace::global>().lower, bounds_xyz.lower, 1e-12);
        CheckPointNear(element.GetActiveDomainBounds<CoordinateSpace::global>().upper, bounds_xyz.upper, 1e-12);
        CheckPointNear(element.GetActiveDomainBounds<CoordinateSpace::parametric>().lower, bounds_uvw.lower, 1e-12);
        CheckPointNear(element.GetActiveDomainBounds<CoordinateSpace::parametric>().upper, bounds_uvw.upper, 1e-12);
    }

    TEST_CASE("UntrimmedElementTestSuite - Integration Points")
    {


        auto element = MakeElement();
        constexpr IntegrationPointType ip_a({ -2.5, 4.25, 9.1 }, 1.5);
        constexpr IntegrationPointType ip_b({ 1.0, -1.0, 2.0 }, 0.75);

        auto& r_integration_points = element.GetIntegrationPoints();
        r_integration_points.push_back(ip_a);
        r_integration_points.push_back(ip_b);

        const auto& r_parametric_points = static_cast<const ElementType&>(element).GetIntegrationPoints();
        REQUIRE((r_parametric_points.size()) == (2UL));
        CheckPointNear(r_parametric_points[0].Point(), ip_a.Point(), 1e-12);
        CheckPointNear(r_parametric_points[1].Point(), ip_b.Point(), 1e-12);

        const auto global_points =
            static_cast<const ElementType&>(element).GetIntegrationPoints<CoordinateSpace::global>();
        auto it = global_points.begin();

        const auto global_ip_a = *it++;
        const auto global_ip_b = *it;
        const double det_j = element.DetJ();
        constexpr ElementBounds bounds{ MakeBoundsXYZ(), MakeBoundsUVW() };

        CheckPointNear(global_ip_a.Point(), mapping::ToGlobal(ip_a.Point(), bounds), 1e-12);
        CheckNear(global_ip_a.Weight(), ip_a.Weight() * det_j, 1e-12);
        CheckPointNear(global_ip_b.Point(), mapping::ToGlobal(ip_b.Point(), bounds), 1e-12);
        CheckNear(global_ip_b.Weight(), ip_b.Weight() * det_j, 1e-12);
    }

    TEST_CASE("UntrimmedElementTestSuite - Boundary Mesh And Boundary Ips")
    {


        const auto element = MakeElement();
        constexpr auto bounds_xyz = MakeBoundsXYZ();
        constexpr auto bounds_uvw = MakeBoundsUVW();

        const auto boundary_mesh = element.GetActiveDomainBoundaryMesh();
        CHECK_FALSE(!(boundary_mesh.NumOfTriangles() > 0UL));
        CheckNear(MeshUtilities::Area(boundary_mesh), BoxSurfaceArea(bounds_xyz), 1e-12);

        const auto boundary_ips_global =
            element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::global>();
        CheckNear(SumWeights(boundary_ips_global), BoxSurfaceArea(bounds_xyz), 1e-10);

        const auto boundary_ips_parametric =
            element.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::parametric>();
        CheckNear(SumWeights(boundary_ips_parametric), BoxSurfaceArea(bounds_uvw), 1e-10);
    }

    TEST_CASE("UntrimmedElementTestSuite - Inside And Intersection State")
    {


        const auto element = MakeElement();

        CHECK_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>(PointType{ 0.0, -80.0, 10.0 }));
        CHECK_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::global>(PointType{ -25.0, -111.44, 7.89 }));
        CHECK_FALSE(element.IsInsideActiveDomain<CoordinateSpace::global>(PointType{ -26.0, -80.0, 10.0 }));

        CHECK_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::parametric>(PointType{ -2.0, 4.0, 9.0 }));
        CHECK_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::parametric>(PointType{ -10.0, -2.2, 2.0 }));
        CHECK_FALSE(element.IsInsideActiveDomain<CoordinateSpace::parametric>(PointType{ -10.1, 4.0, 9.0 }));

        CHECK(
            (element.GetIntersectionState<CoordinateSpace::global>(
                PointType{ -20.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }
            ))
            == (queso::IntersectionState::inside)
        );
        CHECK(
            (element.GetIntersectionState<CoordinateSpace::global>(
                PointType{ 100.0, -100.0, 8.5 }, PointType{ 110.0, -90.0, 12.0 }
            ))
            == (queso::IntersectionState::outside)
        );
        CHECK(
            (element.GetIntersectionState<CoordinateSpace::global>(
                PointType{ -30.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }
            ))
            == (queso::IntersectionState::trimmed)
        );

        CHECK(
            (element.GetIntersectionState<CoordinateSpace::parametric>(
                PointType{ -8.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }
            ))
            == (queso::IntersectionState::inside)
        );
        CHECK(
            (element.GetIntersectionState<CoordinateSpace::parametric>(
                PointType{ 3.0, -1.0, 3.0 }, PointType{ 4.0, 8.0, 14.0 }
            ))
            == (queso::IntersectionState::outside)
        );
        CHECK(
            (element.GetIntersectionState<CoordinateSpace::parametric>(
                PointType{ -11.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }
            ))
            == (queso::IntersectionState::trimmed)
        );
    }

    TEST_CASE("UntrimmedElementTestSuite - Det JAnd Data Set")
    {


        auto element = MakeElement();
        constexpr auto bounds_xyz = MakeBoundsXYZ();
        constexpr auto bounds_uvw = MakeBoundsUVW();

        CheckNear(element.DetJ(), ReferenceDetJ(bounds_xyz, bounds_uvw), 1e-12);

        element.SetValue(ElementValues::is_visited, true);
        element.SetValue(ElementValues::neighbor_coefficient, 2.5);

        CHECK((element.GetValue<bool>(ElementValues::is_visited)) == (true));
        CheckNear(element.GetValue<double>(ElementValues::neighbor_coefficient), 2.5, 1e-12);
    }

}  // TEST_SUITE

}  // namespace queso::testing
