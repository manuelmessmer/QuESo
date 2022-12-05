// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "utilities/values.h"

namespace tibra {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ParamtersTestSuite )

BOOST_AUTO_TEST_CASE(ParameterDefaultTest) {
    std::cout << "Testing :: Parameter :: Test Defaults" << std::endl;
    TestParameter parameters{};

    IndexType eche_level = parameters.Get<IndexType>("echo_level");
    BOOST_CHECK_EQUAL(eche_level, 0UL);

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    BOOST_CHECK(embedding_flag);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    BOOST_CHECK_LT(std::abs(initial_triangle_edge_length-1.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<IndexType>("min_num_boundary_triangles");
    BOOST_CHECK_EQUAL(min_num_boundary_triangles, 1000L);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    BOOST_CHECK_LT( std::abs(1.0e-10-moment_fitting_residual)/1.0e-10, 1.0e-10);

    IndexType init_point_distribution_factor = parameters.Get<IndexType>("init_point_distribution_factor");
    BOOST_CHECK_LT( std::abs(1.0e-10-moment_fitting_residual)/1.0e-10, 1.0e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    BOOST_CHECK_EQUAL( polynomial_order[0], 2UL);
    BOOST_CHECK_EQUAL( polynomial_order[1], 2UL);
    BOOST_CHECK_EQUAL( polynomial_order[2], 2UL);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    BOOST_CHECK_EQUAL( integration_method, IntegrationMethod::Gauss);
}

BOOST_AUTO_TEST_CASE(ParameterCustomConstructorTest) {
    std::cout << "Testing :: Parameter :: Test Custom Constructor" << std::endl;

    TestParameter parameters( {Component("input_filename", std::string("date/test.stl")),
                               Component("postprocess_filename", std::string("date/test2.stl")),
                               Component("echo_level", 2UL),
                               Component("embedding_flag", false),
                               Component("lower_bound", PointType(-1.0, 0.0, 1.22)),
                               Component("upper_bound", PointType(1.1, 3.3, 4.4)),
                               Component("polynomial_order", Vector3i(3,2,4)),
                               Component("number_of_knot_spans", Vector3i(10,12,14)),
                               Component("initial_triangle_edge_length", 5.0),
                               Component("min_num_boundary_triangles", 2000UL),
                               Component("moment_fitting_residual", 0.5e-5),
                               Component("init_point_distribution_factor", 5UL),
                               Component("integration_method", IntegrationMethod::ReducedExact) });

    std::string input_filename = parameters.Get<std::string>("input_filename");
    BOOST_CHECK_EQUAL(input_filename, std::string("date/test.stl"));

    std::string postprocess_filename = parameters.Get<std::string>("postprocess_filename");
    BOOST_CHECK_EQUAL(postprocess_filename, std::string("date/test2.stl"));

    IndexType eche_level = parameters.Get<IndexType>("echo_level");
    BOOST_CHECK_EQUAL(eche_level, 2UL);

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    BOOST_CHECK(!embedding_flag);

    PointType lower_bound = parameters.Get<PointType>("lower_bound");
    BOOST_CHECK_LT(std::abs(lower_bound[0]+1.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[1]+0.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[2]-1.22), 1e-10);

    PointType upper_bound = parameters.Get<PointType>("upper_bound");
    BOOST_CHECK_LT(std::abs(upper_bound[0]-1.1), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[1]-3.3), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[2]-4.4), 1e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    BOOST_CHECK_EQUAL(polynomial_order[0], 3Ul);
    BOOST_CHECK_EQUAL(polynomial_order[1], 2Ul);
    BOOST_CHECK_EQUAL(polynomial_order[2], 4Ul);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    BOOST_CHECK_LT(std::abs(initial_triangle_edge_length-5.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<IndexType>("min_num_boundary_triangles");
    BOOST_CHECK_EQUAL(min_num_boundary_triangles, 2000UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    BOOST_CHECK_LT( std::abs(0.5e-5-moment_fitting_residual)/0.5e-5, 1.0e-10);

    IndexType init_point_distribution_factor = parameters.Get<IndexType>("init_point_distribution_factor");
    BOOST_CHECK_EQUAL(init_point_distribution_factor, 5UL);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    BOOST_CHECK_EQUAL( integration_method, IntegrationMethod::ReducedExact);
}

BOOST_AUTO_TEST_CASE(ParameterCustomSetTest) {
    std::cout << "Testing :: Parameter :: Test Custom Set" << std::endl;

    TestParameter parameters{};
    parameters.Set("input_filename", std::string("date/test.stl") );
    parameters.Set("postprocess_filename", std::string("date/test2.stl"));
    parameters.Set("echo_level", 2UL);
    parameters.Set("embedding_flag", false);
    parameters.Set("lower_bound", PointType(-1.0, 0.0, 1.22));
    parameters.Set("upper_bound", PointType(1.1, 3.3, 4.4));
    parameters.Set("polynomial_order", Vector3i(3,2,4));
    parameters.Set("number_of_knot_spans", Vector3i(10,12,14));
    parameters.Set("initial_triangle_edge_length", 5.0);
    parameters.Set("min_num_boundary_triangles", 2000UL);
    parameters.Set("moment_fitting_residual", 0.5e-5);
    parameters.Set("init_point_distribution_factor", 5UL);
    parameters.Set("integration_method", IntegrationMethod::ReducedExact);

    std::string input_filename = parameters.Get<std::string>("input_filename");
    BOOST_CHECK_EQUAL(input_filename, std::string("date/test.stl"));

    std::string postprocess_filename = parameters.Get<std::string>("postprocess_filename");
    BOOST_CHECK_EQUAL(postprocess_filename, std::string("date/test2.stl"));

    IndexType eche_level = parameters.Get<IndexType>("echo_level");
    BOOST_CHECK_EQUAL(eche_level, 2UL);

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    BOOST_CHECK(!embedding_flag);

    PointType lower_bound = parameters.Get<PointType>("lower_bound");
    BOOST_CHECK_LT(std::abs(lower_bound[0]+1.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[1]+0.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[2]-1.22), 1e-10);

    PointType upper_bound = parameters.Get<PointType>("upper_bound");
    BOOST_CHECK_LT(std::abs(upper_bound[0]-1.1), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[1]-3.3), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[2]-4.4), 1e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    BOOST_CHECK_EQUAL(polynomial_order[0], 3Ul);
    BOOST_CHECK_EQUAL(polynomial_order[1], 2Ul);
    BOOST_CHECK_EQUAL(polynomial_order[2], 4Ul);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    BOOST_CHECK_LT(std::abs(initial_triangle_edge_length-5.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<IndexType>("min_num_boundary_triangles");
    BOOST_CHECK_EQUAL(min_num_boundary_triangles, 2000UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    BOOST_CHECK_LT( std::abs(0.5e-5-moment_fitting_residual)/0.5e-5, 1.0e-10);

    IndexType init_point_distribution_factor = parameters.Get<IndexType>("init_point_distribution_factor");
    BOOST_CHECK_EQUAL(init_point_distribution_factor, 5UL);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    BOOST_CHECK_EQUAL( integration_method, IntegrationMethod::ReducedExact);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace tibra

