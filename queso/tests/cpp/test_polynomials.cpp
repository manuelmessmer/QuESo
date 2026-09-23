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
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "queso/utilities/polynomial_utilities.hpp"

// This suite tests the templated Legendre polynomial utilities used by moment-fitting assembly. Related coverage:
// test_moment_fitting_assembly.cpp checks how these basis functions are assembled into matrix/RHS entries.

namespace queso::Testing {

/// @brief Checks numerical orthogonality of two different Legendre basis functions on an interval.
/// @details Uses tensor-product Gauss quadrature with enough points for the selected orders.
template<IndexType Order1, IndexType Order2>
void CheckLegendreOrthogonality()
{
    constexpr double lower_bound = 0.1;
    constexpr double upper_bound = 0.3;
    constexpr double length = upper_bound - lower_bound;

    const auto& ips_1 = IntegrationPointFactory1D::GetGauss(Order1 + 1, IntegrationMethod::gauss);
    const auto& ips_2 = IntegrationPointFactory1D::GetGauss(Order2 + 1, IntegrationMethod::gauss);
    double numerical_integral = 0.0;
    for (const auto& point1 : ips_1) {
        for (const auto& point2 : ips_2) {
            double position1 = point1[0] * length + lower_bound;
            double position2 = point2[0] * length + lower_bound;
            numerical_integral += polynomial::f_x<Order1>(position1, lower_bound, upper_bound) * point1[1]
                                  * polynomial::f_x<Order2>(position2, lower_bound, upper_bound) * point2[1];
        }
    }
    QuESo_CHECK_LT(std::abs(numerical_integral), 1e-12);
}

/// @brief Checks that the analytical antiderivative matches numerical Gauss integration.
template<IndexType Order>
void CheckLegendreIntegral()
{
    constexpr double lower_bound = 0.1;
    constexpr double upper_bound = 0.3;
    constexpr double length = upper_bound - lower_bound;

    const auto& ips = IntegrationPointFactory1D::GetGauss(Order + 1, IntegrationMethod::gauss);
    double numerical_integral = 0.0;
    for (const auto& point : ips) {
        double position = point[0] * length + lower_bound;
        numerical_integral += polynomial::f_x<Order>(position, lower_bound, upper_bound) * point[1] * length;
    }
    double analytical_int = polynomial::f_x_int<Order>(upper_bound, lower_bound, upper_bound)
                            - polynomial::f_x_int<Order>(lower_bound, lower_bound, upper_bound);
    QuESo_CHECK_LT(std::abs(analytical_int - numerical_integral), 1e-12);
}

BOOST_AUTO_TEST_SUITE(PolynomialTestSuite)

BOOST_AUTO_TEST_CASE(LegendreBasisFunctionsAreOrthogonal)
{
    QuESo_INFO << "Testing :: Test Polynomials :: Legendre Polynomials 1" << std::endl;
    CheckLegendreOrthogonality<0, 1>();
    CheckLegendreOrthogonality<0, 2>();
    CheckLegendreOrthogonality<0, 3>();
    CheckLegendreOrthogonality<0, 4>();
    CheckLegendreOrthogonality<0, 5>();
    CheckLegendreOrthogonality<0, 6>();
    CheckLegendreOrthogonality<0, 7>();
    CheckLegendreOrthogonality<0, 8>();
    CheckLegendreOrthogonality<1, 2>();
    CheckLegendreOrthogonality<1, 3>();
    CheckLegendreOrthogonality<1, 4>();
    CheckLegendreOrthogonality<1, 5>();
    CheckLegendreOrthogonality<1, 6>();
    CheckLegendreOrthogonality<1, 7>();
    CheckLegendreOrthogonality<1, 8>();
    CheckLegendreOrthogonality<2, 3>();
    CheckLegendreOrthogonality<2, 4>();
    CheckLegendreOrthogonality<2, 5>();
    CheckLegendreOrthogonality<2, 6>();
    CheckLegendreOrthogonality<2, 7>();
    CheckLegendreOrthogonality<2, 8>();
    CheckLegendreOrthogonality<3, 4>();
    CheckLegendreOrthogonality<3, 5>();
    CheckLegendreOrthogonality<3, 6>();
    CheckLegendreOrthogonality<3, 7>();
    CheckLegendreOrthogonality<3, 8>();
    CheckLegendreOrthogonality<4, 5>();
    CheckLegendreOrthogonality<4, 6>();
    CheckLegendreOrthogonality<4, 7>();
    CheckLegendreOrthogonality<4, 8>();
    CheckLegendreOrthogonality<5, 6>();
    CheckLegendreOrthogonality<5, 7>();
    CheckLegendreOrthogonality<5, 8>();
    CheckLegendreOrthogonality<6, 7>();
    CheckLegendreOrthogonality<6, 8>();
    CheckLegendreOrthogonality<7, 8>();
}

BOOST_AUTO_TEST_CASE(LegendreAntiderivativesMatchGaussIntegration)
{
    QuESo_INFO << "Testing :: Test Polynomials :: Legendre Polynomials 2" << std::endl;
    CheckLegendreIntegral<0>();
    CheckLegendreIntegral<1>();
    CheckLegendreIntegral<2>();
    CheckLegendreIntegral<3>();
    CheckLegendreIntegral<4>();
    CheckLegendreIntegral<5>();
    CheckLegendreIntegral<6>();
    CheckLegendreIntegral<7>();
    CheckLegendreIntegral<8>();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
