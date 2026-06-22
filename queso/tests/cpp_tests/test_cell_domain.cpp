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
#include "queso/containers/cell_domain.hpp"
#include "queso/includes/checks.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {

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

}// namespace

BOOST_AUTO_TEST_SUITE(CellDomainTestSuite)

BOOST_AUTO_TEST_CASE(Bounds)
{
    QuESo_INFO << "Testing :: CellDomain :: Bounds" << std::endl;

    const CellDomain cell_domain{};
    constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };
    constexpr auto bounds_xyz = MakeBoundsXYZ();
    constexpr auto bounds_uvw = MakeBoundsUVW();

    QuESo_CHECK_POINT_NEAR(cell_domain.GetBounds<CoordinateSpace::global>(bounds).lower, bounds_xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(cell_domain.GetBounds<CoordinateSpace::global>(bounds).upper, bounds_xyz.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(cell_domain.GetBounds<CoordinateSpace::parametric>(bounds).lower, bounds_uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(cell_domain.GetBounds<CoordinateSpace::parametric>(bounds).upper, bounds_uvw.upper, 1e-12);
}

BOOST_AUTO_TEST_CASE(BoundaryMeshAndBoundaryIpsGlobal)
{
    QuESo_INFO << "Testing :: CellDomain :: BoundaryMeshAndBoundaryIpsGlobal" << std::endl;

    const CellDomain cell_domain{};
    constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };
    constexpr auto bounds_xyz = MakeBoundsXYZ();

    const auto boundary_mesh = cell_domain.GetBoundaryMesh(bounds);
    QuESo_CHECK_IS_FALSE(!(boundary_mesh.NumOfTriangles() > 0UL));
    QuESo_CHECK_NEAR(MeshUtilities::Area(boundary_mesh), BoxSurfaceArea(bounds_xyz), 1e-12);

    const auto boundary_ips = cell_domain.GetBoundaryIps<BoundaryIntegrationPoint, CoordinateSpace::global>(bounds);
    QuESo_CHECK_NEAR(SumWeights(boundary_ips), BoxSurfaceArea(bounds_xyz), 1e-10);
}

BOOST_AUTO_TEST_CASE(BoundaryIpsParametric)
{
    QuESo_INFO << "Testing :: CellDomain :: BoundaryIpsParametric" << std::endl;

    const CellDomain cell_domain{};
    constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };
    constexpr auto bounds_uvw = MakeBoundsUVW();

    const auto boundary_ips = cell_domain.GetBoundaryIps<BoundaryIntegrationPoint, CoordinateSpace::parametric>(bounds);
    QuESo_CHECK_NEAR(SumWeights(boundary_ips), BoxSurfaceArea(bounds_uvw), 1e-10);
}

BOOST_AUTO_TEST_CASE(IsInside)
{
    QuESo_INFO << "Testing :: CellDomain :: IsInside" << std::endl;

    const CellDomain cell_domain{};
    constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };

    QuESo_CHECK_IS_FALSE(!cell_domain.IsInside<CoordinateSpace::global>(PointType{ 0.0, -80.0, 10.0 }, bounds));
    QuESo_CHECK_IS_FALSE(!cell_domain.IsInside<CoordinateSpace::global>(PointType{ -25.0, -111.44, 7.89 }, bounds));
    QuESo_CHECK_IS_FALSE(cell_domain.IsInside<CoordinateSpace::global>(PointType{ -26.0, -80.0, 10.0 }, bounds));

    QuESo_CHECK_IS_FALSE(!cell_domain.IsInside<CoordinateSpace::parametric>(PointType{ -2.0, 4.0, 9.0 }, bounds));
    QuESo_CHECK_IS_FALSE(!cell_domain.IsInside<CoordinateSpace::parametric>(PointType{ -10.0, -2.2, 2.0 }, bounds));
    QuESo_CHECK_IS_FALSE(cell_domain.IsInside<CoordinateSpace::parametric>(PointType{ -10.1, 4.0, 9.0 }, bounds));
}

BOOST_AUTO_TEST_CASE(IntersectionStateChecks)
{
    QuESo_INFO << "Testing :: CellDomain :: IntersectionState" << std::endl;

    const CellDomain cell_domain{};
    constexpr auto bounds = ElementBounds{ MakeBoundsXYZ(), MakeBoundsUVW() };

    QuESo_CHECK_EQUAL(
        cell_domain.GetIntersectionState<CoordinateSpace::global>(
            PointType{ -20.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }, bounds, 0.0
        ),
        queso::IntersectionState::inside
    );
    QuESo_CHECK_EQUAL(
        cell_domain.GetIntersectionState<CoordinateSpace::global>(
            PointType{ 100.0, -100.0, 8.5 }, PointType{ 110.0, -90.0, 12.0 }, bounds, 0.0
        ),
        queso::IntersectionState::outside
    );
    QuESo_CHECK_EQUAL(
        cell_domain.GetIntersectionState<CoordinateSpace::global>(
            PointType{ -30.0, -100.0, 8.5 }, PointType{ 10.0, -70.0, 12.0 }, bounds, 0.0
        ),
        queso::IntersectionState::trimmed
    );

    QuESo_CHECK_EQUAL(
        cell_domain.GetIntersectionState<CoordinateSpace::parametric>(
            PointType{ -8.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }, bounds, 0.0
        ),
        queso::IntersectionState::inside
    );
    QuESo_CHECK_EQUAL(
        cell_domain.GetIntersectionState<CoordinateSpace::parametric>(
            PointType{ 3.0, -1.0, 3.0 }, PointType{ 4.0, 8.0, 14.0 }, bounds, 0.0
        ),
        queso::IntersectionState::outside
    );
    QuESo_CHECK_EQUAL(
        cell_domain.GetIntersectionState<CoordinateSpace::parametric>(
            PointType{ -11.0, -1.0, 3.0 }, PointType{ 0.0, 8.0, 14.0 }, bounds, 0.0
        ),
        queso::IntersectionState::trimmed
    );
}

BOOST_AUTO_TEST_SUITE_END()

}// End namespace queso::Testing
