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
#include "queso/includes/checks.hpp"
#include "queso/quadrature/moment_fitting.hpp"
#include "queso/tests/cpp_tests/trimmed_element_test_helpers.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/mesh_utilities.h"

// This suite tests public moment-fitting Compute behavior on trimmed elements. Related coverage:
// test_element_trimmed.cpp checks the TrimmedElement container API, test_moment_fitting.cpp checks NNLS solve behavior,
// and test_point_elimination.cpp checks reduced rules on larger trimmed geometries.

namespace queso::Testing {

namespace {

    using namespace TrimmedElementTestHelpers;

    SizeType CountOutsideActiveDomain(const TrimmedElementType& rElement)
    {
        SizeType result = 0;
        for (const auto& rPoint : rElement.GetIntegrationPoints()) {
            if (!rElement.template IsInsideActiveDomain<CoordinateSpace::parametric>(rPoint.Point())) { ++result; }
        }
        return result;
    }

    double IntegratedVolume(const TrimmedElementType& rElement)
    {
        double result = 0.0;
        for (const auto& rPoint : rElement.GetIntegrationPoints()) { result += rPoint.Weight() * rElement.DetJ(); }
        return result;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(MomentFittingComputeTestSuite)

BOOST_AUTO_TEST_CASE(ReplacesExistingIntegrationPoints)
{
    // Verifies public Compute clears previous integration points and returns an active-domain quadrature rule.
    QuESo_INFO << "Testing :: Moment Fitting Compute :: ReplacesExistingIntegrationPoints" << std::endl;

    auto element = MakeTrimmedElement();
    constexpr IntegrationPointType sentinel({ -0.75, -0.75, -0.75 }, 123.0);
    element.GetIntegrationPoints().push_back(sentinel);

    const auto residual =
        quadrature::moment_fitting::Compute(element, { .integration_order = { 2, 2, 2 }, .residual = 1e-8 });

    QuESo_CHECK(residual.has_value());
    QuESo_CHECK_NEAR(residual.value(), 0.0, 1e-6);

    const auto& r_points = element.GetIntegrationPoints();
    QuESo_CHECK_IS_FALSE(r_points.empty());

    for (const auto& rPoint : r_points) {
        QuESo_CHECK_GT(rPoint.Weight(), EPS4);
        QuESo_CHECK_IS_FALSE(!element.IsInsideActiveDomain<CoordinateSpace::parametric>(rPoint.Point()));
        QuESo_CHECK_IS_FALSE(Math::Norm(rPoint.Point() - sentinel.Point()) < 1e-12);
    }

    const double reference_volume = MeshUtilities::Volume(element.GetActiveDomainBoundaryMesh());
    QuESo_CHECK_RELATIVE_NEAR(IntegratedVolume(element), reference_volume, 1e-4);
}

BOOST_AUTO_TEST_CASE(FictitiousDomainAppendsOutsidePoints)
{
    // Verifies fictitious-domain quadrature appends only outside-active-domain points.
    QuESo_INFO << "Testing :: Moment Fitting Compute :: FictitiousDomainAppendsOutsidePoints" << std::endl;

    auto element_without_alpha = MakeTrimmedElement();
    auto element_with_alpha = MakeTrimmedElement();

    const auto residual_without_alpha = quadrature::moment_fitting::Compute(
        element_without_alpha, { .integration_order = { 2, 2, 2 }, .residual = 1e-8 }
    );
    const auto residual_with_alpha = quadrature::moment_fitting::Compute(
        element_with_alpha, { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 0.25 }
    );

    QuESo_CHECK(residual_without_alpha.has_value());
    QuESo_CHECK(residual_with_alpha.has_value());
    QuESo_CHECK_NEAR(residual_without_alpha.value(), 0.0, 1e-6);
    QuESo_CHECK_NEAR(residual_with_alpha.value(), 0.0, 1e-6);

    const auto& r_points_without_alpha = element_without_alpha.GetIntegrationPoints();
    const auto& r_points_with_alpha = element_with_alpha.GetIntegrationPoints();

    QuESo_CHECK_EQUAL(CountOutsideActiveDomain(element_without_alpha), 0UL);
    QuESo_CHECK_GT(r_points_with_alpha.size(), r_points_without_alpha.size());
    QuESo_CHECK_EQUAL(
        CountOutsideActiveDomain(element_with_alpha), r_points_with_alpha.size() - r_points_without_alpha.size()
    );
}

BOOST_AUTO_TEST_CASE(RejectsInvalidFictitiousDomainAlpha)
{
    // Verifies public Compute rejects fictitious-domain weight scales outside 0 < alpha <= 1.
    QuESo_INFO << "Testing :: Moment Fitting Compute :: RejectsInvalidFictitiousDomainAlpha" << std::endl;

    auto negative_alpha_element = MakeTrimmedElement();
    BOOST_REQUIRE_THROW(
        quadrature::moment_fitting::Compute(
            negative_alpha_element,
            { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = -0.25 }
        ),
        std::exception
    );

    auto zero_alpha_element = MakeTrimmedElement();
    BOOST_REQUIRE_THROW(
        quadrature::moment_fitting::Compute(
            zero_alpha_element, { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 0.0 }
        ),
        std::exception
    );

    auto too_large_alpha_element = MakeTrimmedElement();
    BOOST_REQUIRE_THROW(
        quadrature::moment_fitting::Compute(
            too_large_alpha_element,
            { .integration_order = { 2, 2, 2 }, .residual = 1e-8, .fictitious_domain_alpha = 1.25 }
        ),
        std::exception
    );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
