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
#include "queso/containers/cell_domain.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::testing {

namespace {

    constexpr BoundingBoxType MakeBoundsXYZ()
    { return MakeBox({ -25.0, -111.44, 7.89 }, { 78.67, -35.68, 18.99 }); }

    constexpr BoundingBoxType MakeBoundsUVW()
    { return MakeBox({ -10.0, -2.2, 2.0 }, { 2.0, 10.0, 17.0 }); }

    constexpr double BoxSurfaceArea(const BoundingBoxType& rBounds)
    {
        const double dx = rBounds.upper[0] - rBounds.lower[0];
        const double dy = rBounds.upper[1] - rBounds.lower[1];
        const double dz = rBounds.upper[2] - rBounds.lower[2];
        return 2.0 * (dx * dy + dx * dz + dy * dz);
    }

    double SumWeights(const std::vector<BoundaryIntegrationPoint>& rIps)
    {
        double sum = 0.0;
        for (const auto& rIp : rIps) { sum += rIp.Weight(); }
        return sum;
    }

}  // namespace

TEST_SUITE("CellDomainTestSuite")
{

    TEST_CASE("CellDomainTestSuite - Bounds")
    {


        const CellDomain cell_domain{};
        constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };
        constexpr auto bounds_xyz = MakeBoundsXYZ();
        constexpr auto bounds_uvw = MakeBoundsUVW();

        CheckPointNear(cell_domain.GetBounds<CoordinateSpace::global>(bounds).lower, bounds_xyz.lower, 1e-12);
        CheckPointNear(cell_domain.GetBounds<CoordinateSpace::global>(bounds).upper, bounds_xyz.upper, 1e-12);
        CheckPointNear(cell_domain.GetBounds<CoordinateSpace::parametric>(bounds).lower, bounds_uvw.lower, 1e-12);
        CheckPointNear(cell_domain.GetBounds<CoordinateSpace::parametric>(bounds).upper, bounds_uvw.upper, 1e-12);
    }

    TEST_CASE("CellDomainTestSuite - Boundary Mesh And Boundary Ips Global")
    {


        const CellDomain cell_domain{};
        constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };
        constexpr auto bounds_xyz = MakeBoundsXYZ();

        const auto boundary_mesh = cell_domain.GetBoundaryMesh(bounds);
        CHECK_FALSE(!(boundary_mesh.NumOfTriangles() > 0UL));
        CheckNear(MeshUtilities::Area(boundary_mesh), BoxSurfaceArea(bounds_xyz), 1e-12);

        const auto boundary_ips = cell_domain.GetBoundaryIps<BoundaryIntegrationPoint, CoordinateSpace::global>(bounds);
        CheckNear(SumWeights(boundary_ips), BoxSurfaceArea(bounds_xyz), 1e-10);
    }

    TEST_CASE("CellDomainTestSuite - Boundary Ips Parametric")
    {


        const CellDomain cell_domain{};
        constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };
        constexpr auto bounds_uvw = MakeBoundsUVW();

        const auto boundary_ips =
            cell_domain.GetBoundaryIps<BoundaryIntegrationPoint, CoordinateSpace::parametric>(bounds);
        CheckNear(SumWeights(boundary_ips), BoxSurfaceArea(bounds_uvw), 1e-10);
    }

    TEST_CASE("CellDomainTestSuite - Is Inside")
    {


        const CellDomain cell_domain{};
        constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };

        CHECK_FALSE(!cell_domain.IsInside<CoordinateSpace::global>(PointType{ 0.0, -80.0, 10.0 }, bounds));
        CHECK_FALSE(!cell_domain.IsInside<CoordinateSpace::global>(PointType{ -25.0, -111.44, 7.89 }, bounds));
        CHECK_FALSE(cell_domain.IsInside<CoordinateSpace::global>(PointType{ -26.0, -80.0, 10.0 }, bounds));

        CHECK_FALSE(!cell_domain.IsInside<CoordinateSpace::parametric>(PointType{ -2.0, 4.0, 9.0 }, bounds));
        CHECK_FALSE(!cell_domain.IsInside<CoordinateSpace::parametric>(PointType{ -10.0, -2.2, 2.0 }, bounds));
        CHECK_FALSE(cell_domain.IsInside<CoordinateSpace::parametric>(PointType{ -10.1, 4.0, 9.0 }, bounds));
    }

    TEST_CASE("CellDomainTestSuite - Intersection State Checks")
    {


        const CellDomain cell_domain{};
        constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };

        CHECK(
            (cell_domain.GetIntersectionState<CoordinateSpace::global>(
                PointType{ -20.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }, bounds
            ))
            == (queso::IntersectionState::inside)
        );
        CHECK(
            (cell_domain.GetIntersectionState<CoordinateSpace::global>(
                PointType{ 100.0, -100.0, 8.5 }, PointType{ 110.0, -90.0, 12.0 }, bounds
            ))
            == (queso::IntersectionState::outside)
        );
        CHECK(
            (cell_domain.GetIntersectionState<CoordinateSpace::global>(
                PointType{ -30.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }, bounds
            ))
            == (queso::IntersectionState::trimmed)
        );

        CHECK(
            (cell_domain.GetIntersectionState<CoordinateSpace::parametric>(
                PointType{ -8.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }, bounds
            ))
            == (queso::IntersectionState::inside)
        );
        CHECK(
            (cell_domain.GetIntersectionState<CoordinateSpace::parametric>(
                PointType{ 3.0, -1.0, 3.0 }, PointType{ 4.0, 8.0, 14.0 }, bounds
            ))
            == (queso::IntersectionState::outside)
        );
        CHECK(
            (cell_domain.GetIntersectionState<CoordinateSpace::parametric>(
                PointType{ -11.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }, bounds
            ))
            == (queso::IntersectionState::trimmed)
        );
    }

}  // TEST_SUITE

}  // namespace queso::testing
