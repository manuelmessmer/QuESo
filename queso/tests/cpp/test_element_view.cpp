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
#include "queso/containers/trimmed_element.hpp"
#include "queso/containers/untrimmed_element.hpp"
#include "queso/embedding/mesh_operator.h"
#include "queso/io/io_utilities.h"
#include "queso/tests/cpp/helpers/test_data.hpp"
#include "queso/tests/cpp/helpers/trimmed_domain_test_helpers.hpp"

namespace queso::testing {

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
        const std::string stl_path = TestDataDirectory() + "/primitives/cylinder.stl";
        static const TriangleMesh mesh = [&]() {
            TriangleMesh result;
            IO::ReadMeshFromSTL(result, stl_path.c_str());
            return result;
        }();
        constexpr auto xyz = MakeTrimmedCellBoundsXYZ();
        const embedding::MeshOperator mesh_operator(mesh.View(), TrimmedDomainTestHelpers::MakeTolerance(xyz));
        auto p_domain = TrimmedDomainTestHelpers::MakeTrimmedDomain(mesh_operator, xyz, 100);
        REQUIRE(p_domain != nullptr);

        return TrimmedElementType(
            23, ElementBounds{ MakeTrimmedCellBoundsXYZ(), MakeTrimmedCellBoundsUVW() }, std::move(*p_domain)
        );
    }

}  // namespace

TEST_SUITE("ElementViewTestSuite")
{

    TEST_CASE("ElementViewTestSuite - Identity And Trimmed Flag")
    {


        const auto untrimmed = MakeUntrimmedElement();
        const auto view_u = untrimmed.View();
        CHECK((view_u.GetId()) == (17UL));
        CHECK_FALSE(view_u.IsTrimmed());

        const auto trimmed = MakeTrimmedElement();
        const auto view_t = trimmed.View();
        CHECK((view_t.GetId()) == (23UL));
        CHECK_FALSE(!view_t.IsTrimmed());
    }

    TEST_CASE("ElementViewTestSuite - Cell Bounds")
    {


        const auto element = MakeUntrimmedElement();
        const auto view = element.View();
        constexpr auto xyz = MakeUntrimmedBoundsXYZ();
        constexpr auto uvw = MakeUntrimmedBoundsUVW();

        CheckPointNear(view.GetCellBounds<CoordinateSpace::global>().lower, xyz.lower, 1e-12);
        CheckPointNear(view.GetCellBounds<CoordinateSpace::global>().upper, xyz.upper, 1e-12);
        CheckPointNear(view.GetCellBounds<CoordinateSpace::parametric>().lower, uvw.lower, 1e-12);
        CheckPointNear(view.GetCellBounds<CoordinateSpace::parametric>().upper, uvw.upper, 1e-12);
    }

    TEST_CASE("ElementViewTestSuite - Det J")
    {


        const auto untrimmed = MakeUntrimmedElement();
        const auto view_u = untrimmed.View();
        CheckNear(view_u.DetJ(), untrimmed.DetJ(), 1e-14);

        const auto trimmed = MakeTrimmedElement();
        const auto view_t = trimmed.View();
        CheckNear(view_t.DetJ(), trimmed.DetJ(), 1e-14);
    }

    TEST_CASE("ElementViewTestSuite - Get Integration Points Parametric")
    {


        auto element = MakeUntrimmedElement();
        constexpr IntegrationPointType ip_a({ -2.5, 4.25, 9.1 }, 1.5);
        constexpr IntegrationPointType ip_b({ 1.0, -1.0, 2.0 }, 0.75);
        constexpr IntegrationPointType ip_c({ -9.0, 8.0, 15.0 }, 2.0);

        element.GetIntegrationPoints().push_back(ip_a);
        element.GetIntegrationPoints().push_back(ip_b);
        element.GetIntegrationPoints().push_back(ip_c);

        const auto view = element.View();
        const auto& r_param = view.GetIntegrationPoints<CoordinateSpace::parametric>();

        REQUIRE((r_param.size()) == (3UL));
        CheckPointNear(r_param[0].Point(), ip_a.Point(), 1e-12);
        CheckNear(r_param[0].Weight(), ip_a.Weight(), 1e-12);
        CheckPointNear(r_param[1].Point(), ip_b.Point(), 1e-12);
        CheckNear(r_param[1].Weight(), ip_b.Weight(), 1e-12);
        CheckPointNear(r_param[2].Point(), ip_c.Point(), 1e-12);
        CheckNear(r_param[2].Weight(), ip_c.Weight(), 1e-12);
    }

    TEST_CASE("ElementViewTestSuite - Get Integration Points Global")
    {


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

        CheckPointNear(global_a.Point(), mapping::ToGlobal(ip_a.Point(), bounds), 1e-12);
        CheckNear(global_a.Weight(), ip_a.Weight() * det_j, 1e-12);
        CheckPointNear(global_b.Point(), mapping::ToGlobal(ip_b.Point(), bounds), 1e-12);
        CheckNear(global_b.Weight(), ip_b.Weight() * det_j, 1e-12);
    }

    TEST_CASE("ElementViewTestSuite - View Is Copyable")
    {


        auto element = MakeUntrimmedElement();
        constexpr IntegrationPointType ip({ -2.5, 4.25, 9.1 }, 1.5);
        element.GetIntegrationPoints().push_back(ip);

        const auto view_original = element.View();
        const auto view_copy = view_original;  // NOLINT(performance-unnecessary-copy-initialization)

        CHECK((view_copy.GetId()) == (view_original.GetId()));
        CHECK_FALSE(view_copy.IsTrimmed() != view_original.IsTrimmed());
        CheckNear(view_copy.DetJ(), view_original.DetJ(), 1e-15);
        CheckPointNear(
            view_copy.GetCellBounds<CoordinateSpace::global>().lower,
            view_original.GetCellBounds<CoordinateSpace::global>().lower,
            1e-15
        );

        const auto& r_param_orig = view_original.GetIntegrationPoints<CoordinateSpace::parametric>();
        const auto& r_param_copy = view_copy.GetIntegrationPoints<CoordinateSpace::parametric>();
        REQUIRE((r_param_orig.size()) == (r_param_copy.size()));
        CheckPointNear(r_param_orig[0].Point(), r_param_copy[0].Point(), 1e-15);
    }

}  // TEST_SUITE

}  // namespace queso::testing
