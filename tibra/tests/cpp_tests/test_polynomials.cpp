// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "utilities/polynomial_utilities.h"

namespace Testing{

BOOST_AUTO_TEST_SUITE( PolynomialTestSuite )

BOOST_AUTO_TEST_CASE(PolynomialsTestLegendrePolynomials1) {
    std::cout << "Testing :: Test Polynomials :: Legendre Polynomials 1" << std::endl;
    for(int order = 1; order <= 9; ++order){
        for( int order2 = 1; order2 <= 9; ++order2){
            if( order != order2){
                const auto ips_1 = IntegrationPointFactory1D::GetGauss(order, IntegrationPointFactory1D::IntegrationMethod::Gauss);
                const auto ips_2 = IntegrationPointFactory1D::GetGauss(order2, IntegrationPointFactory1D::IntegrationMethod::Gauss);
                double numerical_integral = 0.0;
                for( auto& point1 : *ips_1){
                    for( auto& point2 : *ips_2){
                        double position1 = point1[0]* 0.2 + 0.1;
                        double position2 = point2[0]* 0.2 + 0.1;
                        numerical_integral += Polynomial::f_x(position1, order-1, 0.1, 0.3)*point1[1]
                            *Polynomial::f_x(position2, order2-1, 0.1, 0.3)*point2[1];
                    }
                }
                BOOST_CHECK_LT(std::abs(numerical_integral), 1e-12);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(PolynomialsTestLegendrePolynomials2) {
    std::cout << "Testing :: Test Polynomials :: Legendre Polynomials 2" << std::endl;
    for(int order = 1; order <= 9; ++order){

        const auto ips = IntegrationPointFactory1D::GetGauss(order, IntegrationPointFactory1D::IntegrationMethod::Gauss);
        double numerical_integral = 0.0;
        for( auto& point : *ips){
            double position = point[0]* 0.2 + 0.1;
            numerical_integral += Polynomial::f_x(position, order-1, 0.1, 0.3)*point[1]*0.2;
        }
        double analytical_int = Polynomial::f_x_int(0.3,order-1, 0.1, 0.3) - Polynomial::f_x_int(0.1,order-1, 0.1, 0.3);
        BOOST_CHECK_LT(std::abs(analytical_int-numerical_integral), 1e-12);
    }
}

BOOST_AUTO_TEST_SUITE_END()

}// Namespace Testing