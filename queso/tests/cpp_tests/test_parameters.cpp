// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "utilities/parameters.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ParamtersTestSuite )

BOOST_AUTO_TEST_CASE(ParameterCheckComponentsTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Wrong Types" << std::endl;
    /// Constructor Tests
    // Check if wrong types are actually detected.
    BOOST_CHECK_THROW( Parameters parameters( {Component("lower_bound_xyz", Vector3i(1, 2, 3))} ), std::exception );
    BOOST_CHECK_THROW( Parameters parameters( {Component("lower_bound_xyz", false)} ), std::exception );
    BOOST_CHECK_THROW( Parameters parameters( {Component("lower_bound_xyz", 1.0)} ), std::exception );
    BOOST_CHECK_THROW( Parameters parameters( {Component("lower_bound_xyz", 1UL)} ), std::exception );
    BOOST_CHECK_THROW( Parameters parameters( {Component("lower_bound_xyz", IntegrationMethod::Gauss)} ), std::exception );
    BOOST_CHECK_THROW( Parameters parameters( {Component("lower_bound_xyz", std::string("Hallo"))} ), std::exception );

    BOOST_CHECK_THROW( Parameters parameters( {Component("integration_method", 1UL)} ), std::exception );
    BOOST_CHECK_THROW( Parameters parameters( {Component("integration_method", std::string("Hallo"))} ), std::exception );

    // Check if wrong Names (typos) are detected.
    BOOST_CHECK_THROW( Parameters parameters( {Component("integration_metod", IntegrationMethod::Gauss)} ), std::exception );

    /// Set function Tests
    // Check if wrong types are detected.
    Parameters parameters{};
    BOOST_CHECK_THROW( parameters.Set<bool>("moment_fitting_residual", true), std::exception);
    BOOST_CHECK_THROW( parameters.Set<unsigned long>("moment_fitting_residual", 1UL), std::exception);
    BOOST_CHECK_THROW( parameters.Set<IntegrationMethodType>("moment_fitting_residual", IntegrationMethod::Gauss), std::exception);
    BOOST_CHECK_THROW( parameters.Set<double>("echo_level", 1.0), std::exception);
    BOOST_CHECK_THROW( parameters.Set<IntegrationMethodType>("echo_level", IntegrationMethod::Gauss), std::exception);
    // Check if wrong Names (typos) are detected.
    BOOST_CHECK_THROW( parameters.Set<unsigned long>("echo_leve", 1UL), std::exception);

    /// Check if error for non-default types is detected.
    Parameters parameters2{};
    BOOST_CHECK_THROW( parameters2.Get<PointType>("upper_bound_xyz"), std::exception );
}

BOOST_AUTO_TEST_CASE(ParameterDefaultTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Defaults" << std::endl;

    Parameters parameters{};

    IndexType eche_level = parameters.Get<unsigned long>("echo_level");
    BOOST_CHECK_EQUAL(eche_level, 0UL);
    BOOST_CHECK_EQUAL(eche_level, parameters.EchoLevel());

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    BOOST_CHECK(embedding_flag);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    BOOST_CHECK_LT(std::abs(initial_triangle_edge_length-1.0),1e-10);
    BOOST_CHECK_LT(std::abs(parameters.InitialTriangleEdgeLength()-1.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<unsigned long>("min_num_boundary_triangles");
    BOOST_CHECK_EQUAL(min_num_boundary_triangles, 500UL);
    BOOST_CHECK_EQUAL(parameters.MinimumNumberOfTriangles(), 500UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    BOOST_CHECK_LT( std::abs(1.0e-10-moment_fitting_residual)/1.0e-10, 1.0e-10);
    BOOST_CHECK_LT( std::abs(1.0e-10-parameters.MomentFittingResidual())/1.0e-10, 1.0e-10);

    IndexType init_point_distribution_factor = parameters.Get<unsigned long>("init_point_distribution_factor");
    BOOST_CHECK_EQUAL( init_point_distribution_factor, 1UL);
    BOOST_CHECK_EQUAL( parameters.GetPointDistributionFactor(), 1UL);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    BOOST_CHECK_EQUAL( polynomial_order[0], 2UL);
    BOOST_CHECK_EQUAL( polynomial_order[1], 2UL);
    BOOST_CHECK_EQUAL( polynomial_order[2], 2UL);
    BOOST_CHECK_EQUAL( parameters.Order()[0], 2UL);
    BOOST_CHECK_EQUAL( parameters.Order()[1], 2UL);
    BOOST_CHECK_EQUAL( parameters.Order()[2], 2UL);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    BOOST_CHECK_EQUAL( integration_method, IntegrationMethod::Gauss);
    BOOST_CHECK_EQUAL( parameters.IntegrationMethod(), IntegrationMethod::Gauss);

    bool use_customized_trimmed_points = parameters.Get<bool>("use_customized_trimmed_points");
    BOOST_CHECK( !use_customized_trimmed_points );
}

BOOST_AUTO_TEST_CASE(ParameterCustomConstructorTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Custom Constructor" << std::endl;

    Parameters parameters( {Component("input_filename", std::string("date/test.stl")),
                               Component("postprocess_filename", std::string("date/test2.stl")),
                               Component("echo_level", 2UL),
                               Component("embedding_flag", false),
                               Component("lower_bound_xyz", PointType(-1.0, 0.0, 1.22)),
                               Component("upper_bound_xyz", PointType(1.1, 3.3, 4.4)),
                               Component("lower_bound_uvw", PointType(-1.0, 0.0, 1.22)),
                               Component("upper_bound_uvw", PointType(1.1, 3.3, 4.4)),
                               Component("polynomial_order", Vector3i(3,2,4)),
                               Component("number_of_elements", Vector3i(10,12,14)),
                               Component("initial_triangle_edge_length", 5.0),
                               Component("min_num_boundary_triangles", 2000UL),
                               Component("moment_fitting_residual", 0.5e-5),
                               Component("init_point_distribution_factor", 5UL),
                               Component("integration_method", IntegrationMethod::GGQ_Optimal) });

    std::string input_filename = parameters.Get<std::string>("input_filename");
    BOOST_CHECK_EQUAL(input_filename, std::string("date/test.stl"));

    std::string postprocess_filename = parameters.Get<std::string>("postprocess_filename");
    BOOST_CHECK_EQUAL(postprocess_filename, std::string("date/test2.stl"));

    IndexType eche_level = parameters.Get<unsigned long>("echo_level");
    BOOST_CHECK_EQUAL(eche_level, 2UL);

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    BOOST_CHECK(!embedding_flag);

    PointType lower_bound = parameters.Get<PointType>("lower_bound_xyz");
    BOOST_CHECK_LT(std::abs(lower_bound[0]+1.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[1]+0.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[2]-1.22), 1e-10);

    PointType upper_bound = parameters.Get<PointType>("upper_bound_xyz");
    BOOST_CHECK_LT(std::abs(upper_bound[0]-1.1), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[1]-3.3), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[2]-4.4), 1e-10);

    PointType lower_bound_uvw = parameters.Get<PointType>("lower_bound_uvw");
    BOOST_CHECK_LT(std::abs(lower_bound_uvw[0]+1.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound_uvw[1]+0.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound_uvw[2]-1.22), 1e-10);

    PointType upper_bound_uvw = parameters.Get<PointType>("upper_bound_uvw");
    BOOST_CHECK_LT(std::abs(upper_bound_uvw[0]-1.1), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound_uvw[1]-3.3), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound_uvw[2]-4.4), 1e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    BOOST_CHECK_EQUAL(polynomial_order[0], 3Ul);
    BOOST_CHECK_EQUAL(polynomial_order[1], 2Ul);
    BOOST_CHECK_EQUAL(polynomial_order[2], 4Ul);

    Vector3i number_elements = parameters.Get<Vector3i>("number_of_elements");
    BOOST_CHECK_EQUAL(number_elements[0], 10Ul);
    BOOST_CHECK_EQUAL(number_elements[1], 12Ul);
    BOOST_CHECK_EQUAL(number_elements[2], 14Ul);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    BOOST_CHECK_LT(std::abs(initial_triangle_edge_length-5.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<unsigned long>("min_num_boundary_triangles");
    BOOST_CHECK_EQUAL(min_num_boundary_triangles, 2000UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    BOOST_CHECK_LT( std::abs(0.5e-5-moment_fitting_residual)/0.5e-5, 1.0e-10);

    IndexType init_point_distribution_factor = parameters.Get<unsigned long>("init_point_distribution_factor");
    BOOST_CHECK_EQUAL(init_point_distribution_factor, 5UL);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    BOOST_CHECK_EQUAL( integration_method, IntegrationMethod::GGQ_Optimal);

    bool use_customized_trimmed_points = parameters.Get<bool>("use_customized_trimmed_points");
    BOOST_CHECK( !use_customized_trimmed_points );
}

BOOST_AUTO_TEST_CASE(ParameterCustomSetTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Custom Set" << std::endl;

    Parameters parameters{};
    parameters.Set("input_filename", std::string("date/test.stl") );
    parameters.Set("postprocess_filename", std::string("date/test2.stl"));
    parameters.Set("echo_level", 2UL);
    parameters.Set("embedding_flag", false);
    parameters.Set("lower_bound_xyz", PointType(-1.0, 0.0, 1.22));
    parameters.Set("upper_bound_xyz", PointType(1.1, 3.3, 4.4));
    parameters.Set("lower_bound_uvw", PointType(-2.0, 1.0, 2.22));
    parameters.Set("upper_bound_uvw", PointType(2.1, 6.3, 6.4));
    parameters.Set("polynomial_order", Vector3i(3,2,4));
    parameters.Set("number_of_elements", Vector3i(10,12,14));
    parameters.Set("initial_triangle_edge_length", 5.0);
    parameters.Set("min_num_boundary_triangles", 2000UL);
    parameters.Set("moment_fitting_residual", 0.5e-5);
    parameters.Set("init_point_distribution_factor", 5UL);
    parameters.Set("integration_method", IntegrationMethod::GGQ_Optimal);
    parameters.Set("use_customized_trimmed_points", true);

    std::string input_filename = parameters.Get<std::string>("input_filename");
    BOOST_CHECK_EQUAL(input_filename, std::string("date/test.stl"));

    std::string postprocess_filename = parameters.Get<std::string>("postprocess_filename");
    BOOST_CHECK_EQUAL(postprocess_filename, std::string("date/test2.stl"));

    IndexType eche_level = parameters.Get<unsigned long>("echo_level");
    BOOST_CHECK_EQUAL(eche_level, 2UL);

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    BOOST_CHECK(!embedding_flag);

    PointType lower_bound = parameters.Get<PointType>("lower_bound_xyz");
    BOOST_CHECK_LT(std::abs(lower_bound[0]+1.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[1]+0.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound[2]-1.22), 1e-10);

    PointType upper_bound = parameters.Get<PointType>("upper_bound_xyz");
    BOOST_CHECK_LT(std::abs(upper_bound[0]-1.1), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[1]-3.3), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound[2]-4.4), 1e-10);

    PointType lower_bound_uvw = parameters.LowerBoundUVW();
    BOOST_CHECK_LT(std::abs(lower_bound_uvw[0]+2.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound_uvw[1]-1.0), 1e-10);
    BOOST_CHECK_LT(std::abs(lower_bound_uvw[2]-2.22), 1e-10);

    PointType upper_bound_uvw = parameters.UpperBoundUVW();
    BOOST_CHECK_LT(std::abs(upper_bound_uvw[0]-2.1), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound_uvw[1]-6.3), 1e-10);
    BOOST_CHECK_LT(std::abs(upper_bound_uvw[2]-6.4), 1e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    BOOST_CHECK_EQUAL(polynomial_order[0], 3Ul);
    BOOST_CHECK_EQUAL(polynomial_order[1], 2Ul);
    BOOST_CHECK_EQUAL(polynomial_order[2], 4Ul);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    BOOST_CHECK_LT(std::abs(initial_triangle_edge_length-5.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<unsigned long>("min_num_boundary_triangles");
    BOOST_CHECK_EQUAL(min_num_boundary_triangles, 2000UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    BOOST_CHECK_LT( std::abs(0.5e-5-moment_fitting_residual)/0.5e-5, 1.0e-10);

    IndexType init_point_distribution_factor = parameters.Get<unsigned long>("init_point_distribution_factor");
    BOOST_CHECK_EQUAL(init_point_distribution_factor, 5UL);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    BOOST_CHECK_EQUAL( integration_method, IntegrationMethod::GGQ_Optimal);

    bool use_customized_trimmed_points = parameters.Get<bool>("use_customized_trimmed_points");
    BOOST_CHECK( use_customized_trimmed_points );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
