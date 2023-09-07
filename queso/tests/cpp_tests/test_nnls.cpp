// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "includes/checks.hpp"
#include "solvers/nnls.h"
#include "includes/logger.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( NNLSTestSuite )

BOOST_AUTO_TEST_CASE(NNLS_Rectangle_Matrix) {
    QuESo_INFO << "Testing :: Test nnls :: Rectangle Matrix" << std::endl;

    nnls::MatrixType A(4, nnls::VectorType(3, 0.0));
    A[0][0] = 0.1; A[0][1] = 0.2; A[0][2] = 3.0;
    A[1][0] = 0.3; A[1][1] = 0.7; A[1][2] = 0.6;
    A[2][0] = 0.2; A[2][1] = 0.5; A[2][2] = 1.1;
    A[3][0] = 0.5; A[3][1] = 0.3; A[3][2] = 0.9;

    nnls::VectorType b(4);
    b[0] = 0.1;
    b[1] = 0.1;
    b[2] = 0.7;
    b[3] = 0.3;

    nnls::VectorType x(3);
    double Rnorm = nnls::nnls(A,b,x);

    nnls::VectorType x_reference(3);
    x_reference[0] = 0.26527761625544183;
    x_reference[1] = 0.39411741902728925;
    x_reference[2] = 0.032491624806329486;

    QuESo_CHECK_NEAR(Rnorm, 0.5080235032184826, 1e-10);
    for( int i = 0; i < 3; ++i)
        QuESo_CHECK_NEAR(x[i], x_reference[i], 1e-10);
}

BOOST_AUTO_TEST_CASE(NNLS_Square_Matrix) {
    QuESo_INFO << "Testing :: Test nnls :: Square Matrix" << std::endl;

    nnls::MatrixType A(4, nnls::VectorType(4, 0.0));
    A[0][0] = 0.2;  A[0][1] = 0.5;  A[0][2] = 6.0; A[0][3] = 0.3;
    A[1][0] = 0.3;  A[1][1] = 0.7;  A[1][2] = 0.6; A[1][3] = 0.3;
    A[2][0] = 0.87; A[2][1] = 0.5;  A[2][2] = 1.1; A[2][3] = 0.2;
    A[3][0] = 0.22; A[3][1] = 0.45; A[3][2] = 0.9; A[3][3] = 1.1;

    nnls::VectorType b(4);
    b[0] = 0.33;
    b[1] = 0.12;
    b[2] = 0.12;
    b[3] = 0.3;

    nnls::VectorType x(3);
    double Rnorm = nnls::nnls(A,b,x);

    nnls::VectorType x_reference(4);
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