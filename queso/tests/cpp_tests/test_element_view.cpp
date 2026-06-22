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
#include "queso/containers/untrimmed_element.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/includes/checks.hpp"
#include "queso/io/io_utilities.h"
#include "queso/tests/cpp_tests/global_config.hpp"

namespace queso::Testing {

namespace {

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using UntrimmedElementType = UntrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;
    using TrimmedElementType = TrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;

    // ---------------------------------------------------------------------------
    // Untrimmed fixture
    // ---------------------------------------------------------------------------

    constexpr BoundingBoxType MakeUntrimmedBoundsXYZ()
    { return MakeBox({ -25.0, -111.44, 7.89 }, { 78.67, -35.68, 18.99 }); }

    constexpr BoundingBoxType MakeUntrimmedBoundsUVW()
    { return MakeBox({ -10.0, -2.2, 2.0 }, { 2.0, 10.0, 17.0 }); }

    UntrimmedElementType MakeUntrimmedElement()
    { return UntrimmedElementType(17, ElementBounds{ MakeUntrimmedBoundsXYZ(), MakeUntrimmedBoundsUVW() }); }

    // ---------------------------------------------------------------------------
    // Trimmed fixture  (cylinder cell — same as test_element_trimmed.cpp)
    // ---------------------------------------------------------------------------

    constexpr BoundingBoxType MakeTrimmedCellBoundsXYZ()
    { return MakeBox({ 0.5, -0.5, 0.0 }, { 1.5, 0.5, 1.0 }); }

    constexpr BoundingBoxType MakeTrimmedCellBoundsUVW()
    { return MakeBox({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 }); }

    TrimmedElementType MakeTrimmedElement()
    {
        const std::string stl_path = GlobalConfig::GetInstance().BaseDir + "/data/cylinder.stl";
        TriangleMesh mesh{};
        IO::ReadMeshFromSTL(mesh, stl_path.c_str());
        BRepOperator brep_op(mesh);

        constexpr auto xyz = MakeTrimmedCellBoundsXYZ();
        BOOST_REQUIRE_EQUAL(brep_op.GetIntersectionState(xyz.lower, xyz.upper), IntersectionState::trimmed);

        auto p_domain = brep_op.pGetTrimmedDomain(xyz.lower, xyz.upper, 0.0, 100);
        BOOST_REQUIRE(p_domain != nullptr);

        return TrimmedElementType(
            23, ElementBounds{ MakeTrimmedCellBoundsXYZ(), MakeTrimmedCellBoundsUVW() }, std::move(*p_domain)
        );
    }

}// namespace

BOOST_AUTO_TEST_SUITE(ElementViewTestSuite)

BOOST_AUTO_TEST_CASE(IdentityAndTrimmedFlag)
{
    QuESo_INFO << "Testing :: ElementView :: IdentityAndTrimmedFlag" << std::endl;

    const auto untrimmed = MakeUntrimmedElement();
    const auto view_u = untrimmed.View();
    QuESo_CHECK_EQUAL(view_u.GetId(), 17UL);
    QuESo_CHECK_IS_FALSE(view_u.IsTrimmed());

    const auto trimmed = MakeTrimmedElement();
    const auto view_t = trimmed.View();
    QuESo_CHECK_EQUAL(view_t.GetId(), 23UL);
    QuESo_CHECK_IS_FALSE(!view_t.IsTrimmed());
}

BOOST_AUTO_TEST_CASE(CellBounds)
{
    QuESo_INFO << "Testing :: ElementView :: CellBounds" << std::endl;

    const auto element = MakeUntrimmedElement();
    const auto view = element.View();
    constexpr auto xyz = MakeUntrimmedBoundsXYZ();
    constexpr auto uvw = MakeUntrimmedBoundsUVW();

    QuESo_CHECK_POINT_NEAR(view.GetCellBounds<CoordinateSpace::global>().lower, xyz.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(view.GetCellBounds<CoordinateSpace::global>().upper, xyz.upper, 1e-12);
    QuESo_CHECK_POINT_NEAR(view.GetCellBounds<CoordinateSpace::parametric>().lower, uvw.lower, 1e-12);
    QuESo_CHECK_POINT_NEAR(view.GetCellBounds<CoordinateSpace::parametric>().upper, uvw.upper, 1e-12);
}

BOOST_AUTO_TEST_CASE(DetJ)
{
    QuESo_INFO << "Testing :: ElementView :: DetJ" << std::endl;

    const auto untrimmed = MakeUntrimmedElement();
    const auto view_u = untrimmed.View();
    QuESo_CHECK_NEAR(view_u.DetJ(), untrimmed.DetJ(), 1e-15);

    const auto trimmed = MakeTrimmedElement();
    const auto view_t = trimmed.View();
    QuESo_CHECK_NEAR(view_t.DetJ(), trimmed.DetJ(), 1e-15);
}

BOOST_AUTO_TEST_CASE(GetIntegrationPointsParametric)
{
    QuESo_INFO << "Testing :: ElementView :: GetIntegrationPointsParametric" << std::endl;

    auto element = MakeUntrimmedElement();
    constexpr IntegrationPointType ip_a({ -2.5, 4.25, 9.1 }, 1.5);
    constexpr IntegrationPointType ip_b({ 1.0, -1.0, 2.0 }, 0.75);
    constexpr IntegrationPointType ip_c({ -9.0, 8.0, 15.0 }, 2.0);

    element.GetIntegrationPoints().push_back(ip_a);
    element.GetIntegrationPoints().push_back(ip_b);
    element.GetIntegrationPoints().push_back(ip_c);

    const auto view = element.View();
    const auto& r_param = view.GetIntegrationPoints<CoordinateSpace::parametric>();

    QuESo_CHECK_EQUAL(r_param.size(), 3UL);
    QuESo_CHECK_POINT_NEAR(r_param[0].Point(), ip_a.Point(), 1e-12);
    QuESo_CHECK_NEAR(r_param[0].Weight(), ip_a.Weight(), 1e-12);
    QuESo_CHECK_POINT_NEAR(r_param[1].Point(), ip_b.Point(), 1e-12);
    QuESo_CHECK_NEAR(r_param[1].Weight(), ip_b.Weight(), 1e-12);
    QuESo_CHECK_POINT_NEAR(r_param[2].Point(), ip_c.Point(), 1e-12);
    QuESo_CHECK_NEAR(r_param[2].Weight(), ip_c.Weight(), 1e-12);
}

BOOST_AUTO_TEST_CASE(GetIntegrationPointsGlobal)
{
    QuESo_INFO << "Testing :: ElementView :: GetIntegrationPointsGlobal" << std::endl;

    auto element = MakeUntrimmedElement();
    constexpr IntegrationPointType ip_a({ -2.5, 4.25, 9.1 }, 1.5);
    constexpr IntegrationPointType ip_b({ 1.0, -1.0, 2.0 }, 0.75);

    element.GetIntegrationPoints().push_back(ip_a);
    element.GetIntegrationPoints().push_back(ip_b);

    const auto view = element.View();
    const auto global_ips = view.GetIntegrationPoints<CoordinateSpace::global>();
    auto it = global_ips.begin();
    const auto global_a = *it++;
    const auto global_b = *it;

    const double det_j = view.DetJ();
    constexpr ElementBounds bounds{ MakeUntrimmedBoundsXYZ(), MakeUntrimmedBoundsUVW() };

    QuESo_CHECK_POINT_NEAR(global_a.Point(), mapping::ToGlobal(ip_a.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(global_a.Weight(), ip_a.Weight() * det_j, 1e-12);
    QuESo_CHECK_POINT_NEAR(global_b.Point(), mapping::ToGlobal(ip_b.Point(), bounds), 1e-12);
    QuESo_CHECK_NEAR(global_b.Weight(), ip_b.Weight() * det_j, 1e-12);
}

BOOST_AUTO_TEST_CASE(ViewIsCopyable)
{
    QuESo_INFO << "Testing :: ElementView :: ViewIsCopyable" << std::endl;

    auto element = MakeUntrimmedElement();
    constexpr IntegrationPointType ip({ -2.5, 4.25, 9.1 }, 1.5);
    element.GetIntegrationPoints().push_back(ip);

    const auto view_original = element.View();
    const auto view_copy = view_original;// NOLINT(performance-unnecessary-copy-initialization)

    QuESo_CHECK_EQUAL(view_copy.GetId(), view_original.GetId());
    QuESo_CHECK_IS_FALSE(view_copy.IsTrimmed() != view_original.IsTrimmed());
    QuESo_CHECK_NEAR(view_copy.DetJ(), view_original.DetJ(), 1e-15);
    QuESo_CHECK_POINT_NEAR(
        view_copy.GetCellBounds<CoordinateSpace::global>().lower,
        view_original.GetCellBounds<CoordinateSpace::global>().lower,
        1e-15
    );

    const auto& r_param_orig = view_original.GetIntegrationPoints<CoordinateSpace::parametric>();
    const auto& r_param_copy = view_copy.GetIntegrationPoints<CoordinateSpace::parametric>();
    QuESo_CHECK_EQUAL(r_param_orig.size(), r_param_copy.size());
    QuESo_CHECK_POINT_NEAR(r_param_orig[0].Point(), r_param_copy[0].Point(), 1e-15);
}

BOOST_AUTO_TEST_SUITE_END()

}// End namespace queso::Testing
