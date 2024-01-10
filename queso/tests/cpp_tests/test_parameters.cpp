// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "includes/checks.hpp"
#include "includes/parameters.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( ParamtersTestSuite )

BOOST_AUTO_TEST_CASE(ParameterCheckComponentsTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Wrong Types" << std::endl;
    /// Constructor Tests
    // Check if wrong types are actually detected.
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
    Parameters parameters2{};
    BOOST_CHECK_THROW( parameters2.Set<double>("echo_level", 1.0), std::exception);
    Parameters parameters3{};
    BOOST_CHECK_THROW( parameters3.Set<IntegrationMethodType>("moment_fitting_residual", IntegrationMethod::Gauss), std::exception);
    Parameters parameters4{};
    BOOST_CHECK_THROW( parameters4.Set<double>("echo_level", 1.0), std::exception);
    Parameters parameters5{};
    BOOST_CHECK_THROW( parameters5.Set<IntegrationMethodType>("echo_level", IntegrationMethod::Gauss), std::exception);
    Parameters parameters6{};
    BOOST_CHECK_THROW( parameters6.Set<Vector3d>("number_of_elements", Vector3d(1.0, 1.0, 1.0)), std::exception);

    // Check if wrong Names (typos) are detected.
    Parameters parameters7{};
    BOOST_CHECK_THROW( parameters7.Set<unsigned long>("echo_leve", 1UL), std::exception);
}

BOOST_AUTO_TEST_CASE(ParameterDefaultTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Defaults" << std::endl;

    Parameters parameters{};

    IndexType eche_level = parameters.Get<unsigned long>("echo_level");
    QuESo_CHECK_EQUAL(eche_level, 1UL);
    QuESo_CHECK_EQUAL(eche_level, parameters.EchoLevel());

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    QuESo_CHECK(embedding_flag);

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    QuESo_CHECK_NEAR(initial_triangle_edge_length, 1.0, 1e-10);
    QuESo_CHECK_NEAR(parameters.InitialTriangleEdgeLength(), 1.0,1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<unsigned long>("min_num_boundary_triangles");
    QuESo_CHECK_EQUAL(min_num_boundary_triangles, 100UL);
    QuESo_CHECK_EQUAL(parameters.MinimumNumberOfTriangles(), 100UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    QuESo_CHECK_RELATIVE_NEAR(moment_fitting_residual, 1e-10, 1e-10);
    QuESo_CHECK_RELATIVE_NEAR(parameters.MomentFittingResidual(), 1e-10, 1e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    QuESo_CHECK_Vector3i_EQUAL(polynomial_order, Vector3i(2, 2, 2));
    QuESo_CHECK_Vector3i_EQUAL(polynomial_order, parameters.Order());

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    QuESo_CHECK_EQUAL( integration_method, IntegrationMethod::Gauss);
    QuESo_CHECK_EQUAL( parameters.IntegrationMethod(), IntegrationMethod::Gauss);

    bool use_customized_trimmed_points = parameters.Get<bool>("use_customized_trimmed_points");
    QuESo_CHECK_IS_FALSE( use_customized_trimmed_points );
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
                            Component("integration_method", IntegrationMethod::GGQ_Optimal) });

    std::string input_filename = parameters.Get<std::string>("input_filename");
    QuESo_CHECK_EQUAL(input_filename, std::string("date/test.stl"));

    std::string postprocess_filename = parameters.Get<std::string>("postprocess_filename");
    QuESo_CHECK_EQUAL(postprocess_filename, std::string("date/test2.stl"));

    IndexType eche_level = parameters.Get<unsigned long>("echo_level");
    QuESo_CHECK_EQUAL(eche_level, 2UL);

    bool embedding_flag = parameters.Get<bool>("embedding_flag");
    QuESo_CHECK_IS_FALSE(embedding_flag);

    PointType lower_bound = parameters.Get<PointType>("lower_bound_xyz");
    QuESo_CHECK_POINT_NEAR(lower_bound, PointType(-1.0, 0.0, 1.22), 1e-10);

    PointType upper_bound = parameters.Get<PointType>("upper_bound_xyz");
    QuESo_CHECK_POINT_NEAR(upper_bound, PointType(1.1, 3.3, 4.4), 1e-10);

    PointType lower_bound_uvw = parameters.Get<PointType>("lower_bound_uvw");
    QuESo_CHECK_POINT_NEAR(lower_bound_uvw, PointType(-1.0, 0.0, 1.22), 1e-10);

    PointType upper_bound_uvw = parameters.Get<PointType>("upper_bound_uvw");
    QuESo_CHECK_POINT_NEAR(upper_bound_uvw, PointType(1.1, 3.3, 4.4), 1e-10);

    Vector3i polynomial_order = parameters.Get<Vector3i>("polynomial_order");
    QuESo_CHECK_Vector3i_EQUAL(polynomial_order, Vector3i(3, 2, 4));

    Vector3i number_elements = parameters.Get<Vector3i>("number_of_elements");
    QuESo_CHECK_Vector3i_EQUAL(number_elements, Vector3i(10, 12, 14));

    double initial_triangle_edge_length = parameters.Get<double>("initial_triangle_edge_length");
    QuESo_CHECK_LT(std::abs(initial_triangle_edge_length-5.0),1e-10);

    IndexType min_num_boundary_triangles = parameters.Get<unsigned long>("min_num_boundary_triangles");
    QuESo_CHECK_EQUAL(min_num_boundary_triangles, 2000UL);

    double moment_fitting_residual = parameters.Get<double>("moment_fitting_residual");
    QuESo_CHECK_LT( std::abs(0.5e-5-moment_fitting_residual)/0.5e-5, 1.0e-10);

    IntegrationMethod integration_method = parameters.Get<IntegrationMethod>("integration_method");
    QuESo_CHECK_EQUAL( integration_method, IntegrationMethod::GGQ_Optimal);

    bool use_customized_trimmed_points = parameters.Get<bool>("use_customized_trimmed_points");
    QuESo_CHECK_IS_FALSE( use_customized_trimmed_points );
}

void CheckValuesParameterCustomSetTest(const Parameters& rParameters ) {
    std::string input_filename = rParameters.Get<std::string>("input_filename");

    QuESo_CHECK_EQUAL(input_filename, std::string("date/test.stl"));

    std::string postprocess_filename = rParameters.Get<std::string>("postprocess_filename");
    QuESo_CHECK_EQUAL(postprocess_filename, std::string("date/test2.stl"));

    IndexType eche_level = rParameters.Get<unsigned long>("echo_level");
    QuESo_CHECK_EQUAL(eche_level, 2UL);

    bool embedding_flag = rParameters.Get<bool>("embedding_flag");
    QuESo_CHECK_IS_FALSE(embedding_flag);

    PointType lower_bound = rParameters.Get<PointType>("lower_bound_xyz");
    QuESo_CHECK_POINT_NEAR(lower_bound, PointType(-1.0, 0.0, 1.22), 1e-10);

    PointType upper_bound = rParameters.Get<PointType>("upper_bound_xyz");
    QuESo_CHECK_POINT_NEAR(upper_bound, PointType(1.1, 3.3, 4.4), 1e-10);

    PointType lower_bound_uvw = rParameters.Get<PointType>("lower_bound_uvw");
    QuESo_CHECK_POINT_NEAR(lower_bound_uvw, PointType(-2.0, 1.0, 2.22), 1e-10);

    PointType upper_bound_uvw = rParameters.Get<PointType>("upper_bound_uvw");
    QuESo_CHECK_POINT_NEAR(upper_bound_uvw, PointType(2.1, 6.3, 6.4), 1e-10);

    Vector3i polynomial_order = rParameters.Get<Vector3i>("polynomial_order");
    QuESo_CHECK_Vector3i_EQUAL(polynomial_order, Vector3i(3, 2, 4));

    double initial_triangle_edge_length = rParameters.Get<double>("initial_triangle_edge_length");
    QuESo_CHECK_LT(std::abs(initial_triangle_edge_length-5.0),1e-10);

    IndexType min_num_boundary_triangles = rParameters.Get<unsigned long>("min_num_boundary_triangles");
    QuESo_CHECK_EQUAL(min_num_boundary_triangles, 2000UL);

    double moment_fitting_residual = rParameters.Get<double>("moment_fitting_residual");
    QuESo_CHECK_LT( std::abs(0.5e-5-moment_fitting_residual)/0.5e-5, 1.0e-10);

    IntegrationMethod integration_method = rParameters.Get<IntegrationMethod>("integration_method");
    QuESo_CHECK_EQUAL( integration_method, IntegrationMethod::GGQ_Optimal);

    bool use_customized_trimmed_points = rParameters.Get<bool>("use_customized_trimmed_points");
    QuESo_CHECK( use_customized_trimmed_points );
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
    parameters.Set("integration_method", IntegrationMethod::GGQ_Optimal);
    parameters.Set("use_customized_trimmed_points", true);

    CheckValuesParameterCustomSetTest(parameters);
    Parameters parameters_copy_const(parameters);
    CheckValuesParameterCustomSetTest(parameters_copy_const);
    Parameters parameters_copy_assign = parameters;
    CheckValuesParameterCustomSetTest(parameters_copy_assign);
}

BOOST_AUTO_TEST_CASE(ParameterCastWrongTypesTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Cast Wrong Types" << std::endl;

    Parameters parameters{};
    // Check if IndexTypes are properly converted to doubles
    parameters.Set<double>("moment_fitting_residual", 3.0);

    parameters.Set<unsigned long>("moment_fitting_residual", 4UL);
    BOOST_CHECK_THROW( parameters.Get<unsigned long>("moment_fitting_residual"), std::exception );
    double value = parameters.Get<double>("moment_fitting_residual");
    QuESo_CHECK_NEAR(value, 4.0, 1e-14);

    parameters.Set<unsigned long>("moment_fitting_residual", 5UL);
    BOOST_CHECK_THROW( parameters.Get<unsigned long>("moment_fitting_residual"), std::exception );
    value = parameters.Get<double>("moment_fitting_residual");
    QuESo_CHECK_NEAR(value, 5.0, 1e-14);

    parameters.Set<double>("moment_fitting_residual", 5UL);
    BOOST_CHECK_THROW( parameters.Get<unsigned long>("moment_fitting_residual"), std::exception );
    value = parameters.Get<double>("moment_fitting_residual");
    QuESo_CHECK_NEAR(value, 5.0, 1e-14);

    // Check if Vector3i are properly converted to Vector3d
    parameters.Set<Vector3d>("lower_bound_xyz", Vector3d(1.0, 1.0, 1.0));

    parameters.Set<Vector3i>("lower_bound_xyz", Vector3i(2, 2, 2));
    BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
    auto vector = parameters.Get<Vector3d>("lower_bound_xyz");
    QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d(2.0, 2.0, 2.0), vector, 1e-14);

    parameters.Set<Vector3i>("lower_bound_xyz", Vector3i(54, 5, 5));
    BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
    vector = parameters.Get<Vector3d>("lower_bound_xyz");
    QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d(54, 5, 5), vector, 1e-14);

    parameters.Set<Vector3d>("lower_bound_xyz", Vector3d(3.0, 4.0, 5.0));
    BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
    vector = parameters.Get<Vector3d>("lower_bound_xyz");
    QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d(3.0, 4.0, 5.0), vector, 1e-14);

    parameters.Set<Vector3d>("lower_bound_xyz", Vector3d(3.0, 4.0, 5.0));
    BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
    vector = parameters.Get<Vector3d>("lower_bound_xyz");
    QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d(3.0, 4.0, 5.0), vector, 1e-14);
}


void CheckValuesConditionParameteters(const Parameters& rParameters){

    for( const auto& r_cond : rParameters.GetConditionsSettingsVector() ){
        if( r_cond.Get<std::string>("type") == "PenaltySupportCondition" ){
            QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
            QuESo_CHECK_RELATIVE_NEAR(r_cond.Get<double>("penalty_factor"), 5.0, 1e-14 );
            QuESo_CHECK_POINT_RELATIVE_NEAR(r_cond.Get<PointType>("value"), PointType(5.0, 3.0, 4.0), 1e-14 );
        }
        if( r_cond.Get<std::string>("type") == "LagrangeSupportCondition" ){
            QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
            QuESo_CHECK_POINT_RELATIVE_NEAR(r_cond.Get<PointType>("value"), PointType(5.0, 3.0, 4.0), 1e-14 );
        }
        if( r_cond.Get<std::string>("type") == "SurfaceLoadCondition" ){
            QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
            QuESo_CHECK_RELATIVE_NEAR(r_cond.Get<double>("modulus"), 10.0, 1e-14 );
            QuESo_CHECK_POINT_RELATIVE_NEAR(r_cond.Get<PointType>("direction"), PointType(5.0, 3.0, 4.0), 1e-14 );
        }
        if( r_cond.Get<std::string>("type") == "PressureLoadCondition" ){
            QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
            QuESo_CHECK_RELATIVE_NEAR(r_cond.Get<double>("modulus"), 15.0, 1e-14 );
        }
    }
}

BOOST_AUTO_TEST_CASE(ConditionParametetersCustomSetTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Condition Parameters Custom Set" << std::endl;

    Parameters parameters{};
    ConditionParameters penalty_support_params("PenaltySupportCondition");
    penalty_support_params.Set("input_filename", std::string("test.stl"));
    penalty_support_params.Set("penalty_factor", 5.0);
    penalty_support_params.Set("value", PointType(5.0, 3.0, 4.0));
    parameters.AddConditionSettings(penalty_support_params);

    ConditionParameters langrange_support_params("LagrangeSupportCondition");
    langrange_support_params.Set("input_filename", std::string("test.stl"));
    langrange_support_params.Set("value", PointType(5.0, 3.0, 4.0));
    parameters.AddConditionSettings(langrange_support_params);

    ConditionParameters surface_load_params("SurfaceLoadCondition");
    surface_load_params.Set("input_filename", std::string("test.stl"));
    surface_load_params.Set("modulus", 10.0);
    surface_load_params.Set("direction", PointType(5.0, 3.0, 4.0));
    parameters.AddConditionSettings(surface_load_params);

    ConditionParameters surface_pressure_params("PressureLoadCondition");
    surface_pressure_params.Set("input_filename", std::string("test.stl"));
    surface_pressure_params.Set("modulus", 15.0);
    parameters.AddConditionSettings(surface_pressure_params);

    CheckValuesConditionParameteters(parameters);
    Parameters parameters_copy_const(parameters);
    CheckValuesConditionParameteters(parameters_copy_const);
    Parameters parameters_copy_assign = parameters;
    CheckValuesConditionParameteters(parameters_copy_assign);
}

BOOST_AUTO_TEST_CASE(ConditionParametetersCustomConstructorTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Condition Parameters Custom Constructor" << std::endl;

    Parameters parameters{};
    ConditionParameters penalty_support_params({Component("type", std::string("PenaltySupportCondition")),
                                                Component("input_filename", std::string("test.stl")),
                                                Component("penalty_factor", 5.0),
                                                Component("value", PointType(5.0, 3.0, 4.0)) });
    parameters.AddConditionSettings(penalty_support_params);

    ConditionParameters langrange_support_params({Component("type", std::string("LagrangeSupportCondition")),
                                                  Component("input_filename", std::string("test.stl")),
                                                  Component("value", PointType(5.0, 3.0, 4.0)) });
    parameters.AddConditionSettings(langrange_support_params);

    ConditionParameters surface_load_params({Component("type", std::string("SurfaceLoadCondition")),
                                             Component("input_filename", std::string("test.stl")),
                                             Component("direction", PointType(5.0, 3.0, 4.0)),
                                             Component("modulus", 10.0) });
    parameters.AddConditionSettings(surface_load_params);

    ConditionParameters surface_pressure_params({Component("type", std::string("PressureLoadCondition")),
                                                 Component("input_filename", std::string("test.stl")),
                                                 Component("modulus", 15.0) });
    parameters.AddConditionSettings(surface_pressure_params);

    CheckValuesConditionParameteters(parameters);
    Parameters parameters_copy_const(parameters);
    CheckValuesConditionParameteters(parameters_copy_const);
    Parameters parameters_copy_assign = parameters;
    CheckValuesConditionParameteters(parameters_copy_assign);
}

BOOST_AUTO_TEST_CASE(ConditionParametetersWrongTypesTest) {
    QuESo_INFO << "Testing :: Test Parameter :: Test Condition Parameters Wrong Types" << std::endl;

    ConditionParameters penalty_support_params("PenaltySupportCondition");
    BOOST_CHECK_THROW( penalty_support_params.Set("penalty_factor", PointType(5.0, 3.0, 4.0)), std::exception );

    ConditionParameters penalty_support_params_2("PenaltySupportCondition");
    BOOST_CHECK_THROW( penalty_support_params_2.Set("modulus", 5.0), std::exception );

    ConditionParameters penalty_support_params_3("PenaltySupportCondition");
    BOOST_CHECK_THROW( penalty_support_params_3.Set("penalty_facto", 5.0), std::exception );
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

