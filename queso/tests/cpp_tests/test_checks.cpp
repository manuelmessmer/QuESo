// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "queso/includes/checks.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( CheckTestSuite )

BOOST_AUTO_TEST_CASE(CheckTestCheckTrue) {
    QuESo_INFO << "Testing :: Test Checks :: Check true" << std::endl;

    QuESo_CHECK(true);
    QuESo_CHECK_IS_FALSE(false);

    BOOST_CHECK_THROW( QuESo_CHECK(false), std::exception  );
    BOOST_CHECK_THROW( QuESo_CHECK_IS_FALSE(true), std::exception  );
}

BOOST_AUTO_TEST_CASE(CheckTestCheckEqual) {
    QuESo_INFO << "Testing :: Test Checks :: Check equal" << std::endl;

    QuESo_CHECK_EQUAL(true, true);
    QuESo_CHECK_NOT_EQUAL(true, false);

    BOOST_CHECK_THROW( QuESo_CHECK_EQUAL(true, false), std::exception  );
    BOOST_CHECK_THROW( QuESo_CHECK_NOT_EQUAL(true, true), std::exception  );

    Vector3i a{1, 2, 3};
    Vector3i b{1, 2, 3};
    QuESo_CHECK_Vector3i_EQUAL(a, b);
    Vector3i b2{1, 3, 3};
    BOOST_CHECK_THROW( QuESo_CHECK_Vector3i_EQUAL(a, b2), std::exception );
}

BOOST_AUTO_TEST_CASE(CheckTestCheckLT) {
    QuESo_INFO << "Testing :: Test Checks :: Check less than" << std::endl;

    QuESo_CHECK_LT(0.1, 0.2);
    BOOST_CHECK_THROW( QuESo_CHECK_LT(0.2, 0.1), std::exception  );
}

BOOST_AUTO_TEST_CASE(CheckTestCheckGT) {
    QuESo_INFO << "Testing :: Test Checks :: Check greater than" << std::endl;

    QuESo_CHECK_GT(0.2, 0.1);
    BOOST_CHECK_THROW( QuESo_CHECK_GT(0.1, 0.2), std::exception  );
}

BOOST_AUTO_TEST_CASE(CheckTestNear) {
    QuESo_INFO << "Testing :: Test Checks :: Check near" << std::endl;
    double a = 1.0;
    double b = 1.0 + 0.9e-6;

    QuESo_CHECK_NEAR(a, b, 1e-6);
    BOOST_CHECK_THROW( QuESo_CHECK_NEAR(a, b, 1e-7), std::exception );

    a = -1.0;
    b = -1.0 + 0.9e-6;
    QuESo_CHECK_NEAR(a, b, 1e-6);
    BOOST_CHECK_THROW( QuESo_CHECK_NEAR(a, b, 1e-7), std::exception );

    PointType a_p{2.2, 3.3, -2.0};
    PointType b_p{2.2+9e-6, 3.3-9e-6, -2.0+9e-6};
    QuESo_CHECK_POINT_NEAR(a_p, b_p, 1e-5);
    PointType b2_p{2.2+0.9e-5, 3.3-0.9e-5, -2.0+0.9e-5};
    BOOST_CHECK_THROW( QuESo_CHECK_POINT_NEAR(a_p, b2_p, 8e-6), std::exception );
}

BOOST_AUTO_TEST_CASE(CheckTestRelativeNear) {
    QuESo_INFO << "Testing :: Test Checks :: Check relative near" << std::endl;

    double a = 100.0 + 0.009;
    double b = 100.0;
    QuESo_CHECK_RELATIVE_NEAR(a, b, 1e-4);
    BOOST_CHECK_THROW( QuESo_CHECK_RELATIVE_NEAR(a, b, 0.8e-4), std::exception );

    a = -100.0 + 0.009;
    b = -100.0;
    QuESo_CHECK_RELATIVE_NEAR(a, b, 1e-4);
    BOOST_CHECK_THROW( QuESo_CHECK_RELATIVE_NEAR(a, b, 0.8e-4), std::exception );

    a = 0.0 + 0.0009;
    b = 0.0;
    QuESo_CHECK_RELATIVE_NEAR(a, b, 1e-3);
    BOOST_CHECK_THROW( QuESo_CHECK_RELATIVE_NEAR(a, b, 1e-4), std::exception );

    PointType a_p{2.2, 3.3, -2.0};
    PointType b_p{2.2+9e-6, 3.3-9e-6, -2.0+9e-6};
    QuESo_CHECK_POINT_RELATIVE_NEAR(a_p, b_p, 1e-5);
    PointType b2_p{2.2+0.9e-5, 3.3-0.9e-5, -2.0+0.9e-5};
    BOOST_CHECK_THROW( QuESo_CHECK_POINT_RELATIVE_NEAR(a_p, b2_p, 1e-6), std::exception );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
