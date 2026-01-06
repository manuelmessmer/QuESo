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
#include "queso/solvers/nnls.h"
#include "queso/includes/logger.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( NNLSTestSuite )

BOOST_AUTO_TEST_CASE(NNLS_Rectangle_Matrix) {
    QuESo_INFO << "Testing :: Test nnls :: Rectangle Matrix" << std::endl;

    NNLS::MatrixType A = {
        0.1, 0.3, 0.2, 0.5,
        0.2, 0.7, 0.5, 0.3,
        3.0, 0.6, 1.1, 0.9};

    NNLS::VectorType b = {0.1, 0.1, 0.7, 0.3};

    NNLS::VectorType x(3);
    double Rnorm = NNLS::solve(A,b,x);

    NNLS::VectorType x_reference(3);
    x_reference[0] = 0.26527761625544183;
    x_reference[1] = 0.39411741902728925;
    x_reference[2] = 0.032491624806329486;

    QuESo_CHECK_NEAR(Rnorm, 0.5080235032184826, 1e-10);
    for( int i = 0; i < 3; ++i)
        QuESo_CHECK_NEAR(x[i], x_reference[i], 1e-10);
}

BOOST_AUTO_TEST_CASE(NNLS_Square_Matrix) {
    QuESo_INFO << "Testing :: Test nnls :: Square Matrix" << std::endl;

    NNLS::MatrixType A = {
        0.2, 0.3, 0.87, 0.22,
        0.5, 0.7, 0.5,  0.45,
        6.0, 0.6, 1.1,  0.9,
        0.3, 0.3, 0.2,  1.1};

    NNLS::VectorType b = {0.33, 0.12, 0.12, 0.3};

    NNLS::VectorType x(3);
    double Rnorm = NNLS::solve(A,b,x);

    NNLS::VectorType x_reference(4);
    x_reference[0] = 0.015692959503716887;
    x_reference[1] = 0.0347636417196625;
    x_reference[2] = 0.0404670415359823;
    x_reference[3] = 0.22225779341177274;

    QuESo_CHECK_NEAR(Rnorm, 0.0, 1e-10);
    for( int i = 0; i < 4; ++i)
        QuESo_CHECK_NEAR(x[i], x_reference[i], 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
