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

//// STL includes
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

//// External includes
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/includes/checks.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso::Testing {

using triangle_utilities::TriangleQuadratureOrder;

namespace {

    constexpr std::array orders{ TriangleQuadratureOrder::Degree1,
                                 TriangleQuadratureOrder::Degree2,
                                 TriangleQuadratureOrder::Degree4,
                                 TriangleQuadratureOrder::Degree6 };
    constexpr std::array<SizeType, 4> rule_sizes{ 1, 3, 6, 12 };
    constexpr PointType p1{ 1.0, 2.0, 3.0 };
    constexpr PointType p2{ 3.0, 2.0, 4.0 };
    constexpr PointType p3{ 1.0, 5.0, 4.0 };
    constexpr TriangleProxy<WithoutNormals> triangle{ p1, p2, p3 };

    struct NonConstructibleBoundaryPoint : BoundaryIntegrationPoint
    {
    };

    struct LvalueOnlyCallback
    {
        SizeType calls = 0;
        void operator()(BoundaryIntegrationPoint&&) &
        { ++calls; }
    };

    struct RvalueOnlyCallback
    {
        void operator()(BoundaryIntegrationPoint&&) &&
        {}
    };

    template<class TPoint, class TCallback>
    constexpr bool can_visit_boundary_ips_v = requires(TCallback&& rCallback) {
        triangle_utilities::VisitBoundaryIps<TPoint>(triangle, orders[0], std::forward<TCallback>(rCallback));
    };

    static_assert(can_visit_boundary_ips_v<BoundaryIntegrationPoint, LvalueOnlyCallback>);
    static_assert(can_visit_boundary_ips_v<BoundaryIntegrationPoint, LvalueOnlyCallback&>);
    static_assert(!can_visit_boundary_ips_v<BoundaryIntegrationPoint, RvalueOnlyCallback>);
    static_assert(!can_visit_boundary_ips_v<BoundaryIntegrationPoint, decltype([](int) {})>);
    static_assert(concepts::BoundaryIntegrationPoint<NonConstructibleBoundaryPoint>);
    static_assert(!can_visit_boundary_ips_v<NonConstructibleBoundaryPoint, decltype([](auto&&) {})>);
    static_assert(can_visit_boundary_ips_v<BoundaryIntegrationPoint, decltype(&BoundaryIntegrationPoint::Normal)>);

    constexpr bool CheckConstantRules()
    {
        for (SizeType i = 0; i < orders.size(); ++i) {
            const auto points = triangle_utilities::GetIntegrationPoints(orders[i]);
            const auto rule = triangle_utilities::detail::GetRule(orders[i]);
            if (triangle_utilities::NumberOfIntegrationPoints(orders[i]) != rule_sizes[i]) { return false; }
            if (points.size() != rule_sizes[i] || rule.shape_functions.size() != points.size()) { return false; }
            for (SizeType j = 0; j < points.size(); ++j) {
                if (rule.shape_functions[j] != triangle_utilities::detail::ShapeFunctionValues(points[j].Point())) {
                    return false;
                }
            }
        }
        return true;
    }

    static_assert(CheckConstantRules());
    static_assert(triangle_utilities::GetIntegrationPoints(TriangleQuadratureOrder::Degree1)[0].Weight() == 0.5);
    static_assert(!std::is_invocable_v<decltype(triangle_utilities::GetIntegrationPoints), IndexType>);
    static_assert(std::is_nothrow_invocable_r_v<
                  Vector3d,
                  decltype(triangle_utilities::Normal),
                  const TriangleProxy<WithoutNormals>&>);
    static_assert(!std::is_invocable_v<decltype(triangle_utilities::Normal), const TriangleProxy<WithNormals>&>);
    static_assert(
        triangle_utilities::detail::ShapeFunctionValues(PointType{ 0.25, 0.5, 0.0 }) == std::array{ 0.25, 0.25, 0.5 }
    );
    constexpr PointType origin{ 0.0, 0.0, 0.0 };
    constexpr PointType x_vertex{ 3.0, 0.0, 0.0 };
    constexpr PointType y_vertex{ 0.0, 3.0, 0.0 };
    static_assert(
        triangle_utilities::Center(TriangleProxy<WithoutNormals>{ origin, x_vertex, y_vertex })
        == PointType{ 1.0, 1.0, 0.0 }
    );

    double Factorial(IndexType N)
    {
        double result = 1.0;
        for (IndexType i = 2; i <= N; ++i) { result *= i; }
        return result;
    }

}  // namespace

// This suite owns direct triangle geometry and quadrature contracts. Mesh/domain suites cover aggregation.
// `test_mapping_utilities.cpp` owns element-parametric transformations, and Python tests cover the bindings.
BOOST_AUTO_TEST_SUITE(TriangleUtilitiesTestSuite)

BOOST_AUTO_TEST_CASE(GeometryIsInvariantUnderVertexPermutations)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: GeometryIsInvariantUnderVertexPermutations" << std::endl;

    constexpr std::array vertices{ p1, p2, p3 };
    constexpr std::array<std::array<IndexType, 3>, 6> permutations{
        { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 0, 2, 1 }, { 2, 1, 0 }, { 1, 0, 2 } }
    };
    constexpr PointType expected_center{ 5.0 / 3.0, 3.0, 11.0 / 3.0 };
    constexpr Vector3d expected_normal{ -3.0 / 7.0, -2.0 / 7.0, 6.0 / 7.0 };
    const double expected_aspect_ratio =
        std::sqrt(13.0) * (std::sqrt(5.0) + std::sqrt(13.0) + std::sqrt(10.0)) / (14.0 * std::sqrt(3.0));

    // Distinct edge lengths and all vertex permutations exercise each normal edge-selection branch.
    for (SizeType i = 0; i < permutations.size(); ++i) {
        const auto& vertex_order = permutations[i];
        const TriangleProxy<WithoutNormals> permuted{ vertices[vertex_order[0]],
                                                      vertices[vertex_order[1]],
                                                      vertices[vertex_order[2]] };
        QuESo_CHECK_NEAR(triangle_utilities::Area(permuted), 3.5, 1e-14);
        QuESo_CHECK_POINT_NEAR(triangle_utilities::Center(permuted), expected_center, 1e-14);
        QuESo_CHECK_NEAR(triangle_utilities::AspectRatio(permuted), expected_aspect_ratio, 1e-14);
        const auto normal = triangle_utilities::Normal(permuted);
        const double sign = i < 3 ? 1.0 : -1.0;
        for (IndexType d = 0; d < 3; ++d) { QuESo_CHECK_NEAR(normal[d], sign * expected_normal[d], 1e-14); }
    }
}

BOOST_AUTO_TEST_CASE(AspectRatioIsOneForEquilateralTriangles)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: AspectRatioIsOneForEquilateralTriangles" << std::endl;

    const PointType top{ 1.5, 1.5 * std::sqrt(3.0), 0.0 };
    const TriangleProxy<WithoutNormals> equilateral{ origin, x_vertex, top };
    QuESo_CHECK_NEAR(triangle_utilities::AspectRatio(equilateral), 1.0, 1e-14);
}

BOOST_AUTO_TEST_CASE(ReferenceRulesIntegratePolynomialMoments)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: ReferenceRulesIntegratePolynomialMoments" << std::endl;

    constexpr std::array<IndexType, 4> degrees{ 1, 2, 4, 6 };
    for (SizeType i = 0; i < orders.size(); ++i) {
        const auto points = triangle_utilities::GetIntegrationPoints(orders[i]);
        QuESo_CHECK_EQUAL(points.size(), rule_sizes[i]);
        for (const auto& point : points) {
            BOOST_CHECK(point[0] >= 0.0 && point[1] >= 0.0 && point[0] + point[1] <= 1.0);
            QuESo_CHECK_EQUAL(point[2], 0.0);
            BOOST_CHECK_GT(point.Weight(), 0.0);
        }
        for (IndexType a = 0; a <= degrees[i]; ++a) {
            for (IndexType b = 0; a + b <= degrees[i]; ++b) {
                double moment = 0.0;
                for (const auto& point : points) {
                    moment += point.Weight() * std::pow(point[0], a) * std::pow(point[1], b);
                }
                const double exact_moment = Factorial(a) * Factorial(b) / Factorial(a + b + 2);
                QuESo_CHECK_NEAR(moment, exact_moment, 1e-14);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(VisitorMapsPointsInRuleOrder)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: VisitorMapsPointsInRuleOrder" << std::endl;

    constexpr Vector3d expected_normal{ -3.0 / 7.0, -2.0 / 7.0, 6.0 / 7.0 };
    // Deliberately disagree with winding: WithNormals must preserve the stored orientation.
    constexpr Vector3d stored_normal{ 1.0, 0.0, 0.0 };
    const TriangleProxy<WithNormals> with_stored_normal{ p1, p2, p3, stored_normal };

    for (SizeType i = 0; i < orders.size(); ++i) {
        const auto local = triangle_utilities::GetIntegrationPoints(orders[i]);
        const auto check_points = [&](const auto& rTriangle, const Vector3d& rExpectedNormal) {
            SizeType count = 0;
            double weight_sum = 0.0;
            triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(
                rTriangle, orders[i], [&](BoundaryIntegrationPoint&& rPoint) {
                    BOOST_REQUIRE_LT(count, local.size());
                    const double xi = local[count][0];
                    const double eta = local[count][1];
                    const PointType expected{ 1.0 + 2.0 * xi, 2.0 + 3.0 * eta, 3.0 + xi + eta };
                    QuESo_CHECK_POINT_NEAR(rPoint.Point(), expected, 1e-14);
                    QuESo_CHECK_POINT_NEAR(rPoint.Normal(), rExpectedNormal, 1e-14);
                    QuESo_CHECK_NEAR(rPoint.Weight(), 7.0 * local[count].Weight(), 1e-14);
                    weight_sum += rPoint.Weight();
                    ++count;
                }
            );
            QuESo_CHECK_EQUAL(count, rule_sizes[i]);
            QuESo_CHECK_NEAR(weight_sum, 3.5, 1e-14);
        };
        check_points(triangle, expected_normal);
        check_points(with_stored_normal, stored_normal);
        // Stored normals are passed through as-is, even if the caller supplies a non-unit vector.
        constexpr Vector3d scaled_normal{ 2.0, 0.0, 0.0 };
        const TriangleProxy<WithNormals> with_scaled_normal{ p1, p2, p3, scaled_normal };
        check_points(with_scaled_normal, scaled_normal);
    }
}

BOOST_AUTO_TEST_CASE(WindingAndDegenerateTriangles)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: WindingAndDegenerateTriangles" << std::endl;

    const TriangleProxy<WithoutNormals> reversed{ p1, p3, p2 };
    const auto center = triangle_utilities::Center(triangle);
    const auto forward_normal = triangle_utilities::Normal(triangle);
    SizeType reversed_count = 0;
    triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(
        reversed, TriangleQuadratureOrder::Degree1, [&](const auto& rPoint) {
            QuESo_CHECK_POINT_NEAR(rPoint.Point(), center, 1e-14);
            QuESo_CHECK_NEAR(rPoint.Weight(), 3.5, 1e-14);
            for (IndexType d = 0; d < 3; ++d) { QuESo_CHECK_NEAR(rPoint.Normal()[d], -forward_normal[d], 1e-14); }
            ++reversed_count;
        }
    );
    QuESo_CHECK_EQUAL(reversed_count, 1);

    constexpr PointType collinear{ 5.0, 2.0, 5.0 };
    const TriangleProxy<WithoutNormals> degenerate{ p1, p2, collinear };
    const TriangleProxy<WithoutNormals> collapsed{ p1, p1, p1 };
    for (const auto& degenerate_triangle : { degenerate, collapsed }) {
        QuESo_CHECK_EQUAL(triangle_utilities::Area(degenerate_triangle), 0.0);
        QuESo_CHECK_POINT_NEAR(triangle_utilities::Normal(degenerate_triangle), origin, 0.0);
        BOOST_CHECK_EQUAL(
            triangle_utilities::AspectRatio(degenerate_triangle), std::numeric_limits<double>::infinity()
        );
        for (SizeType i = 0; i < orders.size(); ++i) {
            SizeType count = 0;
            triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(
                degenerate_triangle, orders[i], [&](const auto& rPoint) {
                    QuESo_CHECK_EQUAL(rPoint.Weight(), 0.0);
                    QuESo_CHECK_POINT_NEAR(rPoint.Normal(), origin, 0.0);
                    ++count;
                }
            );
            QuESo_CHECK_EQUAL(count, rule_sizes[i]);
            constexpr Vector3d stored_normal{ 0.0, 0.0, 1.0 };
            const TriangleProxy<WithNormals> with_normal{
                degenerate_triangle.P1, degenerate_triangle.P2, degenerate_triangle.P3, stored_normal
            };
            count = 0;
            triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(
                with_normal, orders[i], [&](const auto& rPoint) {
                    QuESo_CHECK_EQUAL(rPoint.Weight(), 0.0);
                    QuESo_CHECK_POINT_NEAR(rPoint.Normal(), stored_normal, 0.0);
                    ++count;
                }
            );
            QuESo_CHECK_EQUAL(count, rule_sizes[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(InvalidOrdersAndCallbackExceptions)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: InvalidOrdersAndCallbackExceptions" << std::endl;

    for (const auto order : { static_cast<TriangleQuadratureOrder>(0),
                              static_cast<TriangleQuadratureOrder>(3),
                              static_cast<TriangleQuadratureOrder>(5),
                              static_cast<TriangleQuadratureOrder>(-1) }) {
        BOOST_CHECK_THROW((void)triangle_utilities::GetIntegrationPoints(order), Exception);
        BOOST_CHECK_THROW((void)triangle_utilities::NumberOfIntegrationPoints(order), Exception);
        SizeType count = 0;
        BOOST_CHECK_THROW(
            triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(triangle, order, [&](auto&&) { ++count; }),
            Exception
        );
        QuESo_CHECK_EQUAL(count, 0);
    }
    BOOST_CHECK_THROW(
        triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(
            triangle, TriangleQuadratureOrder::Degree1, [](auto&&) { throw std::runtime_error("callback failed"); }
        ),
        std::runtime_error
    );
}

BOOST_AUTO_TEST_CASE(CallbackInvocationUsesLvaluesAndSupportsMemberPointers)
{
    QuESo_INFO << "Testing :: TriangleUtilities :: CallbackInvocationUsesLvaluesAndSupportsMemberPointers" << std::endl;

    LvalueOnlyCallback callback;
    triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(triangle, orders[1], callback);
    QuESo_CHECK_EQUAL(callback.calls, 3);
    triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(triangle, orders[0], LvalueOnlyCallback{});
    triangle_utilities::VisitBoundaryIps<BoundaryIntegrationPoint>(
        triangle, orders[0], &BoundaryIntegrationPoint::Normal
    );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
